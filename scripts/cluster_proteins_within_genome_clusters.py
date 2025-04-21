#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
import csv
import tempfile
import duckdb
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple
import numba

def set_memory_limit(limit_in_gb):
    current_os = platform.system()
    if current_os == "Linux":
        limit_in_bytes = limit_in_gb * 1024**3
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))

log_level = logging.DEBUG if snakemake.params.debug else logging.INFO
log_file = snakemake.params.log
logging.basicConfig(
    level=log_level,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()
logging.getLogger("numba").setLevel(logging.INFO)

print("========================================================================\n"
      "          Step 14/22: Cluster proteins within genome clusters           \n"
      "========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n"
              "          Step 14/22: Cluster proteins within genome clusters           \n"
              "========================================================================\n")

def map_genomes_proteins_to_ids(conn: duckdb.DuckDBPyConnection):
    # Check if protein_alignments_mapped already exists; if so, skip the function
    exists = conn.execute(
        "SELECT 1 FROM information_schema.tables WHERE table_name = 'protein_alignments_mapped'"
    ).fetchone()
    if exists:
        logger.debug("protein_alignments_mapped already exists, skipping mapping.")
        return
    
    logger.debug("Mapping genome_ids into genome_clusters using existing columns from gene_map")
    conn.execute("""
        CREATE TEMP TABLE genome_ids AS
        SELECT DISTINCT genome_id, genome
        FROM gene_map
    """)
    conn.execute("DROP TABLE IF EXISTS genome_clusters_mapped")
    conn.execute("""
        CREATE TABLE genome_clusters_mapped AS
        SELECT
            gc.*,
            gi.genome
        FROM genome_clusters gc
        LEFT JOIN genome_ids gi ON gc.genome_id = gi.genome_id
    """)
    conn.execute("DROP TABLE genome_clusters")
    conn.execute("ALTER TABLE genome_clusters_mapped RENAME TO genome_clusters")

    logger.debug("Replacing query_protein_id and target_protein_id with actual protein names")
    conn.execute("DROP TABLE IF EXISTS protein_alignments_mapped")
    conn.execute("""
        CREATE TABLE protein_alignments_mapped AS
        SELECT
            pa.query_protein_id,
            gm_q.protein AS query,
            gm_q.genome AS query_genome,
            pa.target_protein_id,
            gm_t.protein AS target,
            gm_t.genome AS target_genome,
            pa.fident,
            pa.alnlen,
            pa.bits,
            pa.qlen,
            pa.tlen
        FROM protein_alignments pa
        JOIN gene_map gm_q ON pa.query_protein_id = gm_q.protein_id
        JOIN gene_map gm_t ON pa.target_protein_id = gm_t.protein_id
    """)
    conn.execute("DROP TABLE protein_alignments")
    conn.execute("ALTER TABLE protein_alignments_mapped RENAME TO protein_alignments")

def connect_and_prepare(db_path: str, min_seq_id: float, cov_fraction: float, max_threads: int, mem_limit: int) -> duckdb.DuckDBPyConnection:
    logger.info(f"Connecting to existing DuckDB at {db_path}")
    conn = duckdb.connect(db_path, read_only=False)

    # Heuristic from DuckDB docs: ~5-10 GB memory per thread
    connection_mem_limit = min(round(int(max_threads * 5)), int(mem_limit * 0.9))
    logger.debug(f"Setting DuckDB threads={max_threads}, memory_limit={connection_mem_limit}GB")
    conn.execute(f"SET threads={max_threads}")
    conn.execute(f"SET memory_limit='{connection_mem_limit}GB'")

    required_tables = ['gene_map', 'genome_clusters', 'protein_alignments']
    for table in required_tables:
        exists = conn.execute(f"""
            SELECT COUNT(*) FROM information_schema.tables 
            WHERE table_name = '{table}'
        """).fetchone()[0]
        if not exists:
            logger.error(f"Required table '{table}' does not exist in DuckDB.")
            sys.exit(1)

    map_genomes_proteins_to_ids(conn)

    logger.info(f"Creating filtered protein_alignments with fident >= {min_seq_id} "
                f"and alignment_coverage >= {cov_fraction}")
    conn.execute(f"""
        CREATE VIEW IF NOT EXISTS protein_alignments_filtered AS
        SELECT
            query, target, fident, alnlen,
            bits, qlen, tlen,
            (alnlen * 1.0) / LEAST(qlen, tlen) AS alignment_coverage
        FROM protein_alignments
        WHERE fident >= {min_seq_id:.6f} 
        AND ((alnlen * 1.0) / LEAST(qlen, tlen)) >= {cov_fraction:.6f}
    """)

    rcount = conn.execute("SELECT COUNT(*) FROM protein_alignments_filtered").fetchone()[0]
    logger.info(f"Filtered alignments has {rcount:,} rows.")

    return conn

@numba.njit
def uf_find(x, parent):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x

@numba.njit
def uf_union(x, y, parent, rank):
    rx = uf_find(x, parent)
    ry = uf_find(y, parent)
    if rx != ry:
        if rank[rx] < rank[ry]:
            parent[rx] = ry
        elif rank[rx] > rank[ry]:
            parent[ry] = rx
        else:
            parent[ry] = rx
            rank[rx] += 1

def create_protein_partition_table(conn: duckdb.DuckDBPyConnection, block_size: int,
                                   cluster_ranks: List[str]):
    """
    Creates protein_partition_map for all proteins in gene_map, even those not in alignments.
    Left joins genome_clusters, assigning each protein a [rank]_cluster if present, or a
    fallback (rank_unassigned) if not. This ensures singletons or missing mappings still
    receive a valid cluster ID.
    """
    logger.debug("Creating 'protein_partition_map' in DuckDB with a left join for all proteins.")
    valid_ranks = []
    for r in cluster_ranks:
        cex = conn.execute(f"""
            SELECT COUNT(*) 
            FROM pragma_table_info('genome_clusters')
            WHERE name='{r}_cluster'
        """).fetchone()[0]
        if cex:
            valid_ranks.append(r)
        else:
            logger.warning(f"Skipping rank={r}, not found in genome_clusters.")

    # Build COALESCE expressions for each rank
    rank_cols = []
    for r in valid_ranks:
        rank_cols.append(f"COALESCE(c.{r}_cluster, '{r}_unassigned') AS {r}_cluster")
    if rank_cols:
        rank_select = ",\n       ".join(rank_cols)
    else:
        rank_select = "NULL AS dummy_col"

    conn.execute("DROP TABLE IF EXISTS protein_partition_map")

    # Step1: create tmp_partition with a LEFT JOIN
    tmp_sql = f"""
        CREATE TEMP TABLE tmp_partition AS
        SELECT DISTINCT
           g.protein
           {(","+rank_select) if rank_cols else ""}
        FROM gene_map g
        LEFT JOIN genome_clusters c ON g.genome = c.genome
    """
    conn.execute("DROP TABLE IF EXISTS tmp_partition")
    conn.execute(tmp_sql)

    # Step2: enumerate with ROW_NUMBER
    enumer_sql = f"""
        CREATE TEMP TABLE tmp_enumer AS
        SELECT
            protein,
            ROW_NUMBER() OVER () - 1 AS row_id
          {(","+ ",".join([r+"_cluster" for r in valid_ranks])) if valid_ranks else ""}
        FROM tmp_partition
        ORDER BY protein
    """
    conn.execute("DROP TABLE IF EXISTS tmp_enumer")
    conn.execute(enumer_sql)

    # Step3: final table => block_id, local_idx
    final_cols = []
    final_cols.append("protein")
    final_cols.append("row_id::BIGINT AS row_id")
    final_cols.append(f"(row_id/{block_size})::BIGINT AS block_id")
    final_cols.append(f"(row_id%{block_size})::BIGINT AS local_idx")
    for r in valid_ranks:
        final_cols.append(f"{r}_cluster")

    final_sql = f"""
        CREATE TABLE protein_partition_map AS
        SELECT {", ".join(final_cols)}
        FROM tmp_enumer
        ORDER BY row_id
    """
    conn.execute(final_sql)

    # Cleanup
    conn.execute("DROP TABLE tmp_partition")
    conn.execute("DROP TABLE tmp_enumer")
    
def load_partition_data(
    conn: duckdb.DuckDBPyConnection,
    cluster_ranks: List[str]
):
    logger.debug("Reading 'protein_partition_map' into Python structures.")
    colinfo = conn.execute("""
        SELECT name FROM pragma_table_info('protein_partition_map')
    """).fetchall()
    colnames = [r[0] for r in colinfo]

    # figure out which ranks exist
    valid_ranks = []
    for r in cluster_ranks:
        if f"{r}_cluster" in colnames:
            valid_ranks.append(r)

    rows = conn.execute("SELECT * FROM protein_partition_map").fetchall()
    idx_protein = colnames.index("protein")
    idx_row_id  = colnames.index("row_id")
    idx_block   = colnames.index("block_id")
    idx_local   = colnames.index("local_idx")

    rank_idx_map={}
    for r in valid_ranks:
        rank_idx_map[r] = colnames.index(f"{r}_cluster")

    n=len(rows)
    protein_id_map={}
    block_of_protein_list = [0]*n
    index_in_block_list   = [0]*n
    rank_cluster_data = {r:{} for r in valid_ranks}
    max_block=0

    for rline in rows:
        protein = rline[idx_protein]
        rid     = rline[idx_row_id]
        blk     = rline[idx_block]
        loc     = rline[idx_local]
        protein_id_map[protein]=rid
        block_of_protein_list[rid]=blk
        index_in_block_list[rid]=loc
        if blk>max_block: max_block=blk
        for r in valid_ranks:
            cval = rline[ rank_idx_map[r] ]
            rank_cluster_data[r][(blk,loc)] = cval

    block_count=max_block+1
    logger.debug(f"Loaded partition data: {n:,} proteins, block_count={block_count:,} valid_ranks={valid_ranks}")
    return (protein_id_map, block_of_protein_list, index_in_block_list, block_count, rank_cluster_data, valid_ranks, n)

def unify_edges_chunkwise(
    conn: duckdb.DuckDBPyConnection,
    protein_id_map: Dict[str,int],
    block_of_protein_list: List[int],
    index_in_block_list: List[int],
    rank_cluster_data: Dict[str,Dict[Tuple[int,int],str]],
    valid_ranks: List[str],
    all_parents: Dict[str,np.ndarray],
    all_ranks: Dict[str,np.ndarray],
    chunk_size: int
):
    logger.debug("Reading edges chunkwise, immediate unify in union-find for each rank.")
    alignment_cur = conn.execute("SELECT query, target FROM protein_alignments_filtered")
    total_edges=0
    chunk_count=0

    while True:
        batch=alignment_cur.fetchmany(chunk_size)
        if not batch:
            break
        chunk_count+=1
        for (q,t) in batch:
            total_edges+=1
            qid=protein_id_map.get(q,None)
            tid=protein_id_map.get(t,None)
            if qid is None or tid is None:
                continue
            bq=block_of_protein_list[qid]
            iq=index_in_block_list[qid]
            bt=block_of_protein_list[tid]
            it=index_in_block_list[tid]
            # check each rank
            for r in valid_ranks:
                cval_q = rank_cluster_data[r].get((bq,iq), None)
                cval_t = rank_cluster_data[r].get((bt,it), None)
                if cval_q is not None and cval_t is not None and cval_q==cval_t:
                    uf_union(qid, tid, all_parents[r], all_ranks[r])
        logger.debug(f"Processed chunk {chunk_count}, edges in chunk={len(batch):,}, total={total_edges:,} so far.")
    logger.info(f"Finished reading edges. {total_edges:,} total edges have been processed.")

def store_final_clusters(conn: duckdb.DuckDBPyConnection,
                         valid_ranks: List[str],
                         all_parents: Dict[str, np.ndarray],
                         proteins_list: List[str],
                         block_of_protein_list: List[int],
                         index_in_block_list: List[int],
                         rank_cluster_data: Dict[str,Dict[Tuple[int,int],str]],
                         chunk_size: int = 2_000_000):
    """
    Writes final <rank>_protein_clusters, but uses the corresponding genome cluster name
    plus an underscore plus a unique integer ID. Singletons each get their own ID as well.
    """
    n = len(proteins_list)
    logger.info("Storing final subclusters to CSV and loading into DuckDB for each rank.")

    # Build subcluster names for each rank
    subcluster_names = {}
    for r in valid_ranks:
        parent = all_parents[r]
        subcluster_names[r] = [None]*n
        # (genomeClusterName, rootIndex) -> subclusterID
        cluster_map = {}
        cluster_counts = {}
        for i in range(n):
            # rank_cluster_data[r] always has an entry since the rank is valid
            gc_name = rank_cluster_data[r][(block_of_protein_list[i], index_in_block_list[i])]
            root_i = uf_find(i, parent)
            # get or assign a new subcluster ID within gc_name
            if (gc_name, root_i) not in cluster_map:
                cluster_counts[gc_name] = cluster_counts.get(gc_name, 0) + 1
                cluster_map[(gc_name, root_i)] = cluster_counts[gc_name]
            sc_id = cluster_map[(gc_name, root_i)]
            # final label is genomeClusterName_subclusterID
            subcluster_names[r][i] = f"{gc_name}_{sc_id}"

    # Write these subclusters to DuckDB tables
    for r in valid_ranks:
        tbl = f"{r}_protein_clusters"
        col = f"{r}_protein_cluster"
        conn.execute(f"DROP TABLE IF EXISTS {tbl}")
        conn.execute(f"CREATE TABLE {tbl} (protein TEXT, {col} TEXT)")
        tmp_file = tempfile.NamedTemporaryFile(delete=False, mode='w', newline='', suffix=".csv")
        writer = csv.writer(tmp_file)
        start = 0
        batch_count = 0
        while start < n:
            batch_count += 1
            end = min(start + chunk_size, n)
            for i in range(start, end):
                writer.writerow([proteins_list[i], subcluster_names[r][i]])
            start = end
            logger.debug(f"Wrote chunk {batch_count} (size {end - start}) for rank {r} to CSV.")
        tmp_file.close()
        conn.execute(f"COPY {tbl} FROM '{tmp_file.name}' (FORMAT 'csv')")
        os.remove(tmp_file.name)
        logger.info(f"Finished loading CSV for rank {r}; total {n:,} rows in {tbl}.")
        
def finalize_protein_clusters(conn: duckdb.DuckDBPyConnection, cluster_ranks: List[str], output_tsv: str):
    """
    Merges <rank>_protein_clusters with genome_clusters, writes a TSV.
    Singletons also have valid cluster names of the form <genomeCluster>_<id>.
    """
    logger.info("Finalizing protein clusters...")
    conn.execute("DROP TABLE IF EXISTS gene_map_extended")
    conn.execute("""
        CREATE TEMP TABLE gene_map_extended AS
        SELECT protein, genome, REGEXP_REPLACE(protein, '_(\\d+)$', '') AS contig
        FROM gene_map
    """)

    base_cols = ["g.protein", "g.genome", "g.contig"]
    join_expr = "FROM gene_map_extended g JOIN genome_clusters gc ON g.genome = gc.genome"
    for rank in cluster_ranks:
        rcol = f"{rank}_cluster"
        base_cols.append(f"gc.{rcol} AS {rank}_genome_cluster")

    base_query = f"SELECT {', '.join(base_cols)} {join_expr}"
    conn.execute("DROP TABLE IF EXISTS final_merged_genome_clusters")
    conn.execute(f"CREATE TEMP TABLE final_merged_genome_clusters AS {base_query}")

    merged_q = "SELECT fmgc.protein, fmgc.genome, fmgc.contig"
    for rank in cluster_ranks:
        merged_q += f", fmgc.{rank}_genome_cluster AS {rank}_genome_cluster"
    merged_q += " FROM final_merged_genome_clusters fmgc"

    for rank in cluster_ranks:
        rtab = f"{rank}_protein_clusters"
        texists = conn.execute(f"""
            SELECT COUNT(*) 
            FROM information_schema.tables 
            WHERE table_name = '{rtab}'
        """).fetchone()[0]
        if texists:
            merged_q += f"\n LEFT JOIN {rtab} AS {rank} ON fmgc.protein = {rank}.protein"
            merged_q = merged_q.replace(
                f"fmgc.{rank}_genome_cluster AS {rank}_genome_cluster",
                f"fmgc.{rank}_genome_cluster, {rank}.{rank}_protein_cluster"
            )

    conn.execute("DROP TABLE IF EXISTS final_merged_allranks")
    conn.execute(f"CREATE TABLE final_merged_allranks AS {merged_q}")

    final_cols = ["genome", "contig", "protein"]
    for r in cluster_ranks:
        final_cols += [f"{r}_genome_cluster"]
        final_cols += [f"{r}_protein_cluster"]
        
    col_expr = []
    for c in final_cols:
        cexists = conn.execute(f"""
            SELECT COUNT(*) 
            FROM pragma_table_info('final_merged_allranks') 
            WHERE name = '{c}'
        """).fetchone()[0]
        if cexists:
            col_expr.append(c)
        else:
            # still must select something for output
            col_expr.append(f"'{c}_MISSING' AS {c}")

    final_sel = f"SELECT {', '.join(col_expr)} FROM final_merged_allranks ORDER BY genome, protein"
    conn.execute(f"COPY ({final_sel}) TO '{output_tsv}' (HEADER, DELIMITER '\\t')")
    rc = conn.execute("SELECT COUNT(*) FROM final_merged_allranks").fetchone()[0]
    logger.info(f"Wrote final clusters ({rc:,} rows) to {output_tsv}")

def select_optimal_sizes(conn: duckdb.DuckDBPyConnection, mem_limit: int, num_workers: int) -> Tuple[int, int]:
    """
    Dynamically select optimal block size and chunk size based on input data characteristics,
    available memory, and number of workers.
    
    This function queries DuckDB for:
      - total number of distinct proteins (from gene_map)
      - total number of alignment edges (from protein_alignments_filtered)
    
    Then, it heuristically sets:
      - block_size: number of proteins per partition block.
      - chunk_size: number of alignment edges to process in one batch.
    
    Returns a tuple (block_size, chunk_size).
    """
    total_proteins = conn.execute("SELECT COUNT(DISTINCT protein) FROM gene_map").fetchone()[0]
    total_edges = conn.execute("SELECT COUNT(*) FROM protein_alignments_filtered").fetchone()[0]
    
    # For this heuristic, assume:
    # - block_size: aim for a block to hold roughly between 10e6 and 100e6 proteins.
    block_size = int(min(max(total_proteins // 100, 10_000_000), 100_000_000))
    
    # - chunk_size: choose a chunk size between 10e6 and 100e6 edges.
    chunk_size = int(min(max(total_edges // 100, 10_000_000), 100_000_000))
    
    logger.debug(f"Selected block_size: {block_size:,}, chunk_size: {chunk_size:,} based on "
                f"total_proteins: {total_proteins:,}, total_edges: {total_edges:,}, "
                f"mem_limit: {mem_limit:,} GB, num_workers: {num_workers}.")
    return block_size, chunk_size

def main():
    output_tsv = snakemake.params.protein_clusters_tsv
    min_seq_id = snakemake.params.min_seq_id
    cov_fraction = snakemake.params.cov_fraction
    cluster_ranks = snakemake.params.cluster_taxa_levels
    db_path = snakemake.params.duckdb
    mem_limit = snakemake.resources.mem
    num_workers = snakemake.threads
    os.environ["NUMEXPR_MAX_THREADS"] = str(num_workers)

    # block_size = 100_000_000
    # # block_size = 5_000 # Debugging
    # chunk_size = 100_000_000
    # # chunk_size = 1_000 # Debugging    
    # Larger chunk size is faster, but requires more memory
    # Smaller chunk size is slower, but requires less memory

    logger.debug(f"Setting memory limit to {mem_limit} GB")
    set_memory_limit(mem_limit)

    logger.info("Protein clustering within genome clusters starting...")

    # 1) Connect
    conn = connect_and_prepare(db_path, min_seq_id, cov_fraction, num_workers, mem_limit)
    block_size, chunk_size = select_optimal_sizes(conn, mem_limit, num_workers)

    # 2) Create partition table => protein_partition_map
    create_protein_partition_table(conn, block_size, cluster_ranks)

    # 3) Load partition data => arrays/dicts
    (protein_id_map,
     block_of_protein_list,
     index_in_block_list,
     block_count,
     rank_cluster_data,
     valid_ranks,
     n_proteins) = load_partition_data(conn, cluster_ranks)

    # 4) Build one union-find per rank
    logger.debug(f"Building union-find arrays for {len(valid_ranks)} ranks, n_proteins={n_proteins:,}")
    all_parents={}
    all_ranks={}
    for r in valid_ranks:
        par = np.arange(n_proteins, dtype=np.int64)
        rak = np.zeros(n_proteins, dtype=np.int64)
        all_parents[r] = par
        all_ranks[r]   = rak

    # 5) unify edges chunk by chunk
    unify_edges_chunkwise(
        conn=conn,
        protein_id_map=protein_id_map,
        block_of_protein_list=block_of_protein_list,
        index_in_block_list=index_in_block_list,
        rank_cluster_data=rank_cluster_data,
        valid_ranks=valid_ranks,
        all_parents=all_parents,
        all_ranks=all_ranks,
        chunk_size=chunk_size
    )

    # 6) produce final <rank>_protein_clusters
    # invert protein_id_map => i->protein
    proteins_list = [None] * n_proteins
    for pstr, pid in protein_id_map.items():
        proteins_list[pid] = pstr

    # Use the new function to assign final labels and store them.
    store_final_clusters(conn, valid_ranks, all_parents, proteins_list,
                        block_of_protein_list, index_in_block_list, rank_cluster_data,
                        chunk_size)

    # 7) Finalize and export
    finalize_protein_clusters(conn, cluster_ranks, output_tsv)
    conn.close()

    logger.info("Protein clustering completed.")

if __name__ == "__main__":
    main()
