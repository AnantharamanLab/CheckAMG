#!/usr/bin/env python3

import os
import sys
import tempfile
import shutil
import logging
from typing import Dict, List
import subprocess
import resource
import platform
import duckdb
import rustworkx as rx
from rustworkx import PyGraph, connected_components
from pathlib import Path
from typing import Literal

FilePath = str | Path
ClusterTaxaLevel = Literal["class", "family", "genus", "species"]

def set_memory_limit(limit_in_gb):
    current_os = platform.system()
    if current_os == "Linux":
        limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))

log_level = logging.DEBUG if snakemake.params.debug else logging.INFO
log_file = snakemake.params.log
logging.basicConfig(
    level=log_level,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(log_file, mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()

print("========================================================================\n         Step 13/22: Cluster genomes using amino-acid identities        \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n         Step 13/22: Cluster genomes using amino-acid identities        \n========================================================================\n")

def compute_sparsity(
    db_conn: duckdb.DuckDBPyConnection,
    table_name: str,
    genome_list: list[str]
) -> float:
    """
    Returns fraction of possible edges (undirected) that exist among a set of genomes.
    The table should have (query_genome_id, target_genome_id, score).
    """
    num_nodes = len(genome_list)
    if num_nodes <= 1:
        return 0.0
    max_edges = (num_nodes * (num_nodes - 1)) // 2
    if max_edges == 0:
        return 0.0

    db_conn.execute("CREATE TEMPORARY TABLE tmp_spar (genome_id INT)")
    for g in genome_list:
        db_conn.execute("INSERT INTO tmp_spar VALUES (?)", [g])

    edge_count = db_conn.execute(f"""
        SELECT COUNT(*) 
        FROM {table_name} f
        JOIN tmp_spar g1 ON f.query_genome_id = g1.genome_id
        JOIN tmp_spar g2 ON f.target_genome_id = g2.genome_id
        WHERE f.query_genome_id < f.target_genome_id
    """).fetchone()[0]

    db_conn.execute("DROP TABLE tmp_spar")

    return min(edge_count / max_edges, 1.0)

def run_mcl_file_based(
    db_conn: duckdb.DuckDBPyConnection,
    select_sql: str,
    inflation: float,
    chunk_size: int,
    out_prefix: str,
    rank: str,
    rank_cluster_counter: Dict[str, int],
    threads: int = 1
) -> Dict[str, str]:
    """
    1) Writes edges to an .abc file from DuckDB (in chunks).
    2) Runs mcxload -> .mci + .tab
    3) Runs mcl
    4) Parses clusters => genome_id->cluster
    """
    abc_path = out_prefix + ".abc"
    mci_path = out_prefix + ".mci"
    tab_path = out_prefix + ".tab"
    clusters_path = out_prefix + ".clusters"

    edge_cursor = db_conn.execute(select_sql)
    with open(abc_path, "w") as f_abc:
        batch = edge_cursor.fetchmany(chunk_size)
        while batch:
            for (qg, tg, score) in batch:
                f_abc.write(f"{qg}\t{tg}\t{score}\n")
            batch = edge_cursor.fetchmany(chunk_size)

    mcxload_cmd = [
        "mcxload",
        "-abc", abc_path,
        "-o", mci_path,
        "-write-tab", tab_path
    ]
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".log", delete=False) as tmp_mcx:
        tmp_mcx_name = tmp_mcx.name
    try:
        subprocess.run(
            mcxload_cmd,
            check=True,
            stdout=open(tmp_mcx_name, "w"),
            stderr=open(tmp_mcx_name, "a")
        )
    finally:
        if os.path.exists(tmp_mcx_name):
            os.remove(tmp_mcx_name)

    mcl_cmd = [
        "mcl", mci_path,
        "-te", str(threads),
        "-I", str(inflation),
        "-o", clusters_path,
        "-use-tab", tab_path
    ]
    with tempfile.NamedTemporaryFile(mode="w", suffix=".log", delete=False) as tmp_mcl:
        tmp_mcl_name = tmp_mcl.name
    try:
        subprocess.run(
            mcl_cmd,
            check=True,
            stdout=open(tmp_mcl_name, "w"),
            stderr=open(tmp_mcl_name, "a")
        )
    finally:
        if os.path.exists(tmp_mcl_name):
            os.remove(tmp_mcl_name)

    cluster_assignments = {}
    with open(clusters_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            rank_cluster_counter[rank] += 1
            c_id = rank_cluster_counter[rank]
            cname = f"{rank}_{c_id}"
            members = line.split()
            for genome in members:
                cluster_assignments[genome] = cname

    os.remove(abc_path)
    os.remove(mci_path)
    os.remove(tab_path)
    os.remove(clusters_path)
    return cluster_assignments

def process_subgraph_rustworkx(
    genome_list: List[str],
    edges: List[tuple]
) -> Dict[str, str]:
    """
    Builds a rustworkx graph for these nodes/edges and returns {node: subcomp_label}.
    Node strings must match edge endpoints. This uses connected_components.
    """
    import rustworkx as rx
    graph = rx.PyGraph()
    idx_map = {}
    for g in genome_list:
        nid = graph.add_node(g)
        idx_map[g] = nid

    for (qg, tg, _) in edges:
        qg_s = str(qg)
        tg_s = str(tg)
        if qg_s in idx_map and tg_s in idx_map:
            graph.add_edge(idx_map[qg_s], idx_map[tg_s], None)

    comps = rx.connected_components(graph)
    assignments = {}
    cid = 0
    for cset in comps:
        cid += 1
        lbl = f"tempC_{cid}"
        for nid in cset:
            node_name = graph[nid]
            assignments[node_name] = lbl
    return assignments

def partial_supernode_merge(
    db_conn: duckdb.DuckDBPyConnection,
    label_map: Dict[str, List[str]],
    node_label: Dict[str, str],
    persistent_table: str,
    chunk_size: int,
    max_super_block: int
) -> bool:
    """
    Merges super-nodes in smaller partial slices. Returns True if merges happened, False if stable.
    label_map: super_lbl -> list_of_nodes
    node_label: node -> super_lbl
    max_super_block: how many super-nodes to handle in one pass
    """
    import rustworkx as rx
    super_list = list(label_map.keys())
    merges = False

    # Partition super_list into small blocks
    subblocks = []
    start_idx = 0
    while start_idx < len(super_list):
        end_idx = min(start_idx + max_super_block, len(super_list))
        subblocks.append(super_list[start_idx:end_idx])
        start_idx = end_idx

    # We'll do a pass for each subblock, joining with the rest
    new_label_map: Dict[str, List[str]] = {}
    for sblA in subblocks:
        # Build a PyGraph for these super-nodes + cross edges
        # Next pass merges them. Then we re-update node_label.
        sgraph = rx.PyGraph()
        sidx_map = {}
        for lbl in sblA:
            idx = sgraph.add_node(lbl)
            sidx_map[lbl] = idx

        # get all nodes for sblA
        db_conn.execute("CREATE TEMPORARY TABLE tmp_superA (genome_id TEXT)")
        for lbl in sblA:
            for node in label_map[lbl]:
                db_conn.execute("INSERT INTO tmp_superA VALUES (?)", [node])

        # We'll unify with all other super-nodes not in sblA
        # but not building a giant adjacency at once. We'll do another loop
        other_labels = [x for x in super_list if x not in sblA]
        subblocksB = []
        idxB = 0
        while idxB < len(other_labels):
            endB = min(idxB + max_super_block, len(other_labels))
            subblocksB.append(other_labels[idxB:endB])
            idxB = endB

        for sblB in subblocksB:
            # Add nodes for sblB if not in sgraph
            tmpB_idx_map = {}
            for lblb in sblB:
                idx2 = sgraph.add_node(lblb)
                tmpB_idx_map[lblb] = idx2

            db_conn.execute("CREATE TEMPORARY TABLE tmp_superB (genome_id TEXT)")
            for lblb in sblB:
                for node in label_map[lblb]:
                    db_conn.execute("INSERT INTO tmp_superB VALUES (?)", [node])

            cross_cur = db_conn.execute(f"""
                SELECT f.query_genome_id, f.target_genome_id
                FROM {persistent_table} f
                JOIN tmp_superA sA ON f.query_genome_id = sA.genome_id
                JOIN tmp_superB sB ON f.target_genome_id = sB.genome_id
            """)
            while True:
                part = cross_cur.fetchmany(chunk_size)
                if not part:
                    break
                for (qq, tt) in part:
                    lq = node_label[qq]
                    lt = node_label[tt]
                    if lq != lt:
                        # add edge in supergraph
                        # lq in sidx_map or tmpB_idx_map
                        if lq in sidx_map:
                            n1 = sidx_map[lq]
                        else:
                            n1 = tmpB_idx_map[lq]
                        if lt in sidx_map:
                            n2 = sidx_map[lt]
                        else:
                            n2 = tmpB_idx_map[lt]
                        sgraph.add_edge(n1, n2, None)
            db_conn.execute("DROP TABLE tmp_superB")

        db_conn.execute("DROP TABLE tmp_superA")

        # run connected_components
        comps = rx.connected_components(sgraph)
        if len(comps) < len(sgraph.nodes()):
            merges = True

        # apply merges
        for comp in comps:
            # pick a single new super-label for the entire component
            # merges all old super-nodes in that comp
            rep_lbl = None
            for node_id in comp:
                oldlbl = sgraph[node_id]
                if rep_lbl is None:
                    rep_lbl = oldlbl
                else:
                    # unify oldlbl into rep_lbl
                    if oldlbl != rep_lbl:
                        # move nodes from oldlbl -> rep_lbl
                        for nd in label_map[oldlbl]:
                            node_label[nd] = rep_lbl
                        label_map[rep_lbl].extend(label_map[oldlbl])
                        label_map[oldlbl].clear()
                        label_map.pop(oldlbl, None)

    return merges

def iterative_merge_supernodes(
    db_conn: duckdb.DuckDBPyConnection,
    label_map: Dict[str, List[str]],
    node_label: Dict[str, str],
    persistent_table: str,
    chunk_size: int
):
    """
    Repeated partial merges of super-nodes in small slices until no merges occur.
    """
    stable = False
    while not stable:
        merges_happened = partial_supernode_merge(
            db_conn=db_conn,
            label_map=label_map,
            node_label=node_label,
            persistent_table=persistent_table,
            chunk_size=chunk_size,
            max_super_block=50000 # 50k super-nodes at a time, lower if needed to reduce memory usage
        )
        if not merges_happened:
            stable = True

def multi_pass_greedy(
    db_conn: duckdb.DuckDBPyConnection,
    persistent_table: str,
    current_genomes: List[int],
    rank: str,
    ephemeral_id: int,
    chunk_size: int,
    block_size: int
) -> Dict[str, str]:
    """
    Multi-pass approach with multiple levels of partial merges. Each cross-block step is further
    subdivided to avoid building giant adjacency at once.
    """
    g_list = [str(g) for g in current_genomes]
    if len(g_list) <= 1:
        if len(g_list) == 1:
            return {g_list[0]: f"{rank}_temp_{ephemeral_id}_1"}
        return []

    # 1) Intra-block subgraphs
    blocks = []
    start_idx = 0
    n = len(g_list)
    while start_idx < n:
        end_idx = min(start_idx + block_size, n)
        blocks.append(g_list[start_idx:end_idx])
        start_idx = end_idx

    block_subassign = []
    block_num = 0
    for bnodes in blocks:
        block_num += 1
        block_edges = []
        db_conn.execute("CREATE TEMPORARY TABLE tmp_blk (genome_id TEXT)")
        for gm in bnodes:
            db_conn.execute("INSERT INTO tmp_blk VALUES (?)", [gm])

        qsql = f"""
            SELECT f.query_genome_id, f.target_genome_id, f.score
            FROM {persistent_table} f
            JOIN tmp_blk b1 ON f.query_genome_id = b1.genome_id
            JOIN tmp_blk b2 ON f.target_genome_id = b2.genome_id
        """
        cur = db_conn.execute(qsql)
        while True:
            part = cur.fetchmany(chunk_size)
            if not part:
                break
            for (qg, tg, sc) in part:
                block_edges.append((str(qg), str(tg), sc))

        db_conn.execute("DROP TABLE tmp_blk")

        subc_map = process_subgraph_rustworkx(bnodes, block_edges)
        mapped = {}
        for g in bnodes:
            mapped[g] = f"Blk{block_num}_{subc_map[g]}"
        block_subassign.append(mapped)

    # combine into node->super_lbl
    node_label = {}
    for d in block_subassign:
        node_label.update(d)

    label_map: Dict[str, List[str]] = {}
    for node, lbl in node_label.items():
        label_map.setdefault(lbl, []).append(node)

    # 2) Iterative merges of super-nodes in smaller partial slices
    iterative_merge_supernodes(
        db_conn=db_conn,
        label_map=label_map,
        node_label=node_label,
        persistent_table=persistent_table,
        chunk_size=chunk_size
    )

    return node_label

def worker_greedy_subgraph(args):
    """
    Parallel worker for a single subgraph with multi-level partial merges.
    """
    (db_path, persistent_table, current_genomes, rank, ephemeral_id, threads, mem_limit) = args
    local_conn = duckdb.connect(db_path, read_only=False)
    conn_mem_limit = max(1, round(mem_limit / threads))
    local_conn.execute(f"SET memory_limit = '{conn_mem_limit}GB'")

    gcount = len(current_genomes)
    threshold_nodes = 100_000
    if gcount <= 1:
        if gcount == 1:
            return {str(current_genomes[0]): f"{rank}_temp_{ephemeral_id}_1"}
        return {}
    if gcount <= threshold_nodes:
        # single pass
        local_conn.execute("CREATE TEMPORARY TABLE tmp_g (genome_id TEXT)")
        for g in current_genomes:
            local_conn.execute("INSERT INTO tmp_g VALUES (?)", [str(g)])
        sub_select = f"""
            SELECT f.query_genome_id, f.target_genome_id, f.score
            FROM {persistent_table} f
            JOIN tmp_g g1 ON f.query_genome_id=g1.genome_id
            JOIN tmp_g g2 ON f.target_genome_id=g2.genome_id
        """
        csize = 100_000_000
        edges = []
        cur = local_conn.execute(sub_select)
        while True:
            part = cur.fetchmany(csize)
            if not part:
                break
            edges.extend(part)
        local_conn.execute("DROP TABLE tmp_g")
        local_conn.close()

        G = PyGraph()
        idx_map = {}
        for gm in current_genomes:
            nid = G.add_node(str(gm))
            idx_map[gm] = nid
        for (qg, tg, _) in edges:
            G.add_edge(idx_map[qg], idx_map[tg], None)
        comps = connected_components(G)
        results = {}
        cid = 0
        for cset in comps:
            cid += 1
            lbl = f"{rank}_temp_{ephemeral_id}_{cid}"
            for nd in cset:
                results[G[nd]] = lbl
        return results
    else:
        # multi-level partial merges
        block_size = 300_000
        out_map = multi_pass_greedy(
            db_conn=local_conn,
            persistent_table=persistent_table,
            current_genomes=current_genomes,
            rank=rank,
            ephemeral_id=ephemeral_id,
            chunk_size=100_000_000,
            block_size=block_size
        )
        local_conn.close()
        return out_map

def hierarchical_clustering(
    db_conn: duckdb.DuckDBPyConnection,
    db_path: str,
    wdir: Path,
    aai_table_name: str,
    all_genomes: list[str],
    cluster_ranks: List[str],
    out_table_name: str,
    sparsity_threshold: float = 0.01,
    chunk_size: int = 10_000_000,
    threads: int = 1,
    mem_limit: int = 100
):
    logger.info("Creating cluster_assignments table.")
    cols = ", ".join(f"{r}_cluster TEXT" for r in cluster_ranks)
    db_conn.execute("DROP TABLE IF EXISTS cluster_assignments")
    db_conn.execute(f"""
        CREATE TABLE cluster_assignments (
            genome_id INT PRIMARY KEY,
            {cols}
        )
    """)

    db_conn.execute("BEGIN")
    for g in all_genomes:
        db_conn.execute("INSERT INTO cluster_assignments (genome_id) VALUES (?)", [g])
    db_conn.execute("COMMIT")

    rank_cluster_counter = {r: 0 for r in cluster_ranks}

    inflation_params = {
        "class": 0.5,
        "family": 1.5,
        "genus": 2.0,
        "species": 4.0
    }

    def get_thresholds(level):
        if level == "class":
            return (1, 0.05, 0.10)
        elif level == "family":
            return (4, 0.10, 0.20)
        elif level == "genus":
            return (8, 0.20, 0.80)
        elif level == "species":
            return (16, 0.50, 0.95)
        else:
            raise ValueError(f"Unknown level: {level}")

    from multiprocessing import Pool

    for i, level in enumerate(cluster_ranks):
        prev_level = cluster_ranks[i-1] if i > 0 else None
        logger.info(f"Clustering genomes at the {level}-level...")

        min_genes, min_shared, min_aai = get_thresholds(level)
        filtered_table = f"filtered_{level}"

        db_conn.execute(f"DROP TABLE IF EXISTS {filtered_table}")
        db_conn.execute(f"""
            CREATE TABLE {filtered_table} AS
            SELECT
              query_genome_id,
              target_genome_id,
              (aai * LEAST(query_shared, target_shared)) AS score
            FROM {aai_table_name}
            WHERE (
                (shared_genes >= {min_genes})
                OR (query_shared >= {min_shared} AND target_shared >= {min_shared})
            )
            AND aai >= {min_aai}
        """)
        logger.debug(f"Created persistent table {filtered_table} for {level} filtering.")

        if not prev_level:
            parent_clusters = [None]
            force_top_greedy = True
        else:
            rows = db_conn.execute(f"""
                SELECT {prev_level}_cluster
                FROM cluster_assignments
                GROUP BY {prev_level}_cluster
                HAVING COUNT(*) > 1
            """).fetchall()
            parent_clusters = [r[0] for r in rows if r[0] is not None]
            force_top_greedy = False

        mcl_dir = os.path.join(wdir, "mcl_clusters", level)
        os.makedirs(mcl_dir, exist_ok=True)

        mcl_count = 0
        greedy_count = 0
        tasks = []
        for parent_clust in parent_clusters:
            if parent_clust is None:
                current_genomes = all_genomes
            else:
                g_rows = db_conn.execute(f"""
                    SELECT genome_id
                    FROM cluster_assignments
                    WHERE {prev_level}_cluster = ?
                """, [parent_clust]).fetchall()
                current_genomes = [r[0] for r in g_rows]

            if len(current_genomes) <= 1:
                if len(current_genomes) == 1:
                    rank_cluster_counter[level] += 1
                    c_id = rank_cluster_counter[level]
                    c_name = f"{level}_{c_id}"
                    db_conn.execute(f"""
                        UPDATE cluster_assignments
                        SET {level}_cluster = ?
                        WHERE genome_id = ?
                    """, [c_name, current_genomes[0]])
                continue

            sp = compute_sparsity(db_conn, filtered_table, [str(g) for g in current_genomes])

            if sp == 0:
                for gm in current_genomes:
                    rank_cluster_counter[level] += 1
                    c_id = rank_cluster_counter[level]
                    c_name = f"{level}_{c_id}"
                    db_conn.execute(f"""
                        UPDATE cluster_assignments
                        SET {level}_cluster = ?
                        WHERE genome_id = ?
                    """, [c_name, gm])
                continue

            if force_top_greedy:
                greedy_count += 1
                ephemeral_id = rank_cluster_counter[level]
                task_args = (
                    str(db_path),
                    filtered_table,
                    current_genomes,
                    level,
                    ephemeral_id,
                    threads,
                    mem_limit
                )
                tasks.append(task_args)
            else:
                if sp <= sparsity_threshold:
                    mcl_count += 1
                    out_prefix = os.path.join(
                        mcl_dir,
                        f"{level}_{parent_clust if parent_clust else 'top'}"
                    )
                    db_conn.execute("CREATE TEMPORARY TABLE tmp_parent (genome_id INT)")
                    for g in current_genomes:
                        db_conn.execute("INSERT INTO tmp_parent VALUES (?)", [g])
                    sub_select = f"""
                        SELECT f.query_genome_id, f.target_genome_id, f.score
                        FROM {filtered_table} f
                        JOIN tmp_parent g1 ON f.query_genome_id = g1.genome_id
                        JOIN tmp_parent g2 ON f.target_genome_id = g2.genome_id
                    """
                    assignments = run_mcl_file_based(
                        db_conn=db_conn,
                        select_sql=sub_select,
                        inflation=inflation_params[level],
                        chunk_size=chunk_size,
                        out_prefix=out_prefix,
                        rank=level,
                        rank_cluster_counter=rank_cluster_counter,
                        threads=threads
                    )
                    db_conn.execute("DROP TABLE tmp_parent")

                    for g, cname in assignments.items():
                        db_conn.execute(f"""
                            UPDATE cluster_assignments
                            SET {level}_cluster = ?
                            WHERE genome_id = ?
                        """, [cname, g])
                else:
                    greedy_count += 1
                    ephemeral_id = rank_cluster_counter[level]
                    task_args = (
                        str(db_path),
                        filtered_table,
                        current_genomes,
                        level,
                        ephemeral_id,
                        threads,
                        mem_limit
                    )
                    tasks.append(task_args)

        logger.info(f"For {level}-level clustering, used MCL on {mcl_count:,} parent clusters, greedy on {greedy_count:,}.")

        if tasks:
            with Pool(processes=threads) as p:
                res_list = p.map(worker_greedy_subgraph, tasks)
            for sub_assignments in res_list:
                # sub_assignments: {node: ephemeral_label}
                lbl_map = {}
                for gid, ephemeral_label in sub_assignments.items():
                    lbl_map.setdefault(ephemeral_label, []).append(gid)
                for ephemeral_lbl, glist in lbl_map.items():
                    rank_cluster_counter[level] += 1
                    c_id = rank_cluster_counter[level]
                    c_name = f"{level}_{c_id}"
                    for gm in glist:
                        db_conn.execute(f"""
                            UPDATE cluster_assignments
                            SET {level}_cluster = ?
                            WHERE genome_id = ?
                        """, [c_name, gm])

        db_conn.execute(f"""
            SELECT genome_id FROM cluster_assignments
            WHERE {level}_cluster IS NULL
        """)
        leftover = db_conn.fetchall()
        for (g,) in leftover:
            rank_cluster_counter[level] += 1
            c_id = rank_cluster_counter[level]
            c_name = f"{level}_{c_id}"
            db_conn.execute(f"""
                UPDATE cluster_assignments
                SET {level}_cluster = ?
                WHERE genome_id = ?
            """, [c_name, g])

        n_clusters = db_conn.execute(f"""
            SELECT COUNT(DISTINCT {level}_cluster) 
            FROM cluster_assignments
        """).fetchone()[0]
        multi_query = f"""
            WITH cluster_sizes AS (
              SELECT {level}_cluster AS c, COUNT(*) AS size
              FROM cluster_assignments
              GROUP BY {level}_cluster
            )
            SELECT COUNT(*) 
            FROM cluster_sizes 
            WHERE size > 1
        """
        multi_count = db_conn.execute(multi_query).fetchone()[0]
        logger.info(f"Multi-member clusters at the {level}-level (no singletons): {multi_count:,}.")
        logger.info(f"Total clusters at the {level}-level (including-singletons): {n_clusters:,}.")

        db_conn.execute(f"DROP TABLE {filtered_table}")

    concatenated_expr = " || '-' || ".join(f"{r}_cluster" for r in cluster_ranks)
    db_conn.execute("ALTER TABLE cluster_assignments ADD COLUMN final_cluster TEXT")
    db_conn.execute(f"UPDATE cluster_assignments SET final_cluster = {concatenated_expr}")

    db_conn.execute(f"DROP TABLE IF EXISTS {out_table_name}")
    db_conn.execute(f"""
        CREATE TABLE {out_table_name} AS
        SELECT 
            genome_id, 
            final_cluster, 
            {', '.join(r + '_cluster' for r in cluster_ranks)}
        FROM cluster_assignments
        ORDER BY genome_id
    """)
    logger.info(f"Saved final genome clusters to DuckDB table '{out_table_name}'.")

def main():
    gene_map_path = snakemake.params.gene_map
    output_dir = Path(snakemake.params.wdir)
    db = Path(snakemake.params.duckdb)
    cluster_ranks = snakemake.params.cluster_taxa_levels
    num_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    os.environ["NUMEXPR_MAX_THREADS"] = str(num_cpus)

    set_memory_limit(mem_limit)
    logger.info("Starting hierarchical clustering of genomes by AAI...")
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Using DuckDB database at {db}")
    conn = duckdb.connect(db)
    conn.execute(f"SET threads={num_cpus}")
    conn.execute(f"SET memory_limit = '{max(1, round((mem_limit / num_cpus)*0.9))}GB'")

    logger.info("Loading/verifying AAI table in DuckDB...")
    conn.execute("SELECT COUNT(*) FROM duckdb_tables() WHERE table_name='final_aai'")

    logger.info("Loading/verifying genome map in DuckDB...")
    conn.execute("SELECT COUNT(*) FROM duckdb_tables() WHERE table_name='gene_map'")

    all_g = conn.execute("SELECT DISTINCT genome_id FROM gene_map").fetchall()
    all_genomes = [r[0] for r in all_g]
    logger.info(f"Found {len(all_genomes):,} distinct genomes.")

    hierarchical_clustering(
        db_conn=conn,
        db_path=db,
        wdir=output_dir,
        aai_table_name='final_aai',
        all_genomes=all_genomes,
        cluster_ranks=cluster_ranks,
        out_table_name="genome_clusters",
        sparsity_threshold=0.2,
        chunk_size=100_000_000,
        threads=num_cpus,
        mem_limit=mem_limit
    )

    conn.execute("DROP TABLE IF EXISTS cluster_assignments")
    conn.execute("DROP TABLE IF EXISTS genome_info")
    conn.close()
    mcl_temp_dir = os.path.join(output_dir, "mcl_clusters")
    shutil.rmtree(mcl_temp_dir, ignore_errors=True)
    logger.info("Hierarchical clustering completed.")

if __name__ == "__main__":
    main()