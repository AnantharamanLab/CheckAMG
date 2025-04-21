#!/usr/bin/env python3

import os
from pathlib import Path
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
pl.enable_string_cache()
import shutil
import duckdb
import gc

def set_memory_limit(limit_in_gb):
    """Sets the maximum address space (virtual memory) for the process on Linux."""
    current_os = platform.system()
    if current_os == "Linux":
        limit_in_bytes = limit_in_gb * 1024**3
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

print("========================================================================\n    Step 12/22: Compute the amino-acid identity (AAI) of all proteins   \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n    Step 12/22: Compute the amino-acid identity (AAI) of all proteins   \n========================================================================\n")

def connect_to_duckdb(db_path: str, tmp_spill_path: Path, num_cpus: int, mem_limit: int) -> duckdb.DuckDBPyConnection:
    """
    Connects to a DuckDB database at the given path and configures settings.
    """
    logger.debug(f"Using DuckDB database at: {db_path}")
    con = duckdb.connect(db_path)
    con.execute("SET preserve_insertion_order = false")
    con.execute(f"SET threads={num_cpus}")
    con.execute(f"SET memory_limit = '{round(int(num_cpus * 5))}GB'") # From DuckDB docs: Aim for 5-10 GB memory per thread
    os.makedirs(tmp_spill_path, exist_ok=True)
    con.execute(f"SET temp_directory='{tmp_spill_path}'")
    con.execute(f"SET max_temp_directory_size='{mem_limit*4}GB'")
    if snakemake.params.debug:
        con.execute("PRAGMA enable_progress_bar;")
        
    return con
        
def table_exists(con: duckdb.DuckDBPyConnection, table_name: str) -> bool:
    """
    Returns True if the given table_name exists in the 'main' schema,
    otherwise False.
    """
    sql = f"""
        SELECT table_name
        FROM information_schema.tables
        WHERE table_schema = 'main'
          AND table_name = '{table_name}'
    """
    result = con.execute(sql).fetchone()
    return (result is not None)

def table_has_n_rows(con: duckdb.DuckDBPyConnection, table_name: str, expected_rows: int) -> bool:
    """
    Checks whether the given table_name has exactly `expected_rows` rows.
    """
    sql = f"SELECT COUNT(*) FROM {table_name}"
    row_count = con.execute(sql).fetchone()[0]
    return (row_count == expected_rows)

def read_gene_map(file) -> pl.LazyFrame:
    """
    Reads the gene-to-genome mapping file and logs the counts.
    """
    cols = {
        "protein": pl.String,
        "genome": pl.String,
        "protein_id": pl.UInt64,
        "genome_id": pl.UInt64
    }
    # Only need ID cols for this script, but we'll need the others for downstream steps
    gene_map = pl.read_csv(file, separator="\t", schema_overrides=cols)
    if "genome_id" not in gene_map.columns or "protein_id" not in gene_map.columns:
        logger.error("Gene map file must contain 'genome_id' and 'protein_id' columns.")
        raise ValueError("Gene map file must contain 'genome_id' and 'protein_id' columns.")
        
    logger.info(f"There are {gene_map.shape[0]:,} proteins from {gene_map.select('genome_id').n_unique():,} genomes.")
    
    gc.collect()
    return gene_map

def check_gene_map_table(con: duckdb.DuckDBPyConnection, gene_map_df: pl.DataFrame):
    """
    Checks if the gene_map table exists in DuckDB and if the row count matches the input dataframe.
    If the table does not exist, it is created.
    """
    if table_exists(con, "gene_map"):
        row_count_gmap = gene_map_df.shape[0] # how many rows to expect
        if table_has_n_rows(con, "gene_map", row_count_gmap):
            logger.debug("gene_map table exists and has the correct row count; skipping re-creation.")
        else:
            logger.debug("gene_map table exists but row count mismatch; dropping and re-creating.")
            con.execute("DROP TABLE IF EXISTS gene_map")
            con.execute("""
            CREATE TABLE gene_map (
                protein STRING,
                protein_id UINTEGER,
                genome STRING,
                genome_id UINTEGER
            )
            """)
            arrow_gmap = gene_map_df.to_arrow()
            con.register("tmp_gene_map", arrow_gmap)
            con.execute("""
                INSERT INTO gene_map
                SELECT protein, protein_id, genome, genome_id
                FROM tmp_gene_map
            """)
            con.unregister("tmp_gene_map")
    else:
        logger.debug("gene_map table does not exist; creating it now.")
        con.execute("""
        CREATE TABLE gene_map (
            protein STRING,
            protein_id UINTEGER,
            genome STRING,
            genome_id UINTEGER
        )
        """)
        arrow_gmap = gene_map_df.to_arrow()
        con.register("tmp_gene_map", arrow_gmap)
        con.execute("""
            INSERT INTO gene_map
            SELECT protein, protein_id, genome, genome_id
            FROM tmp_gene_map
        """)
        con.unregister("tmp_gene_map")
        
def count_lines(file_path: str) -> int:
    """
    Counts lines in a file in a memory-efficient manner by processing fixed-size chunks.
    """
    count = 0
    chunk_size = 1024**3 # 1 GB per chunk
    with open(file_path, 'rb') as f:
        while chunk := f.read(chunk_size):
            count += chunk.count(b'\n')
    return count

def load_alignments_in_duckdb(alignments_file: str, con: duckdb.DuckDBPyConnection, gene_map_table: str, min_cov: float):
    """
    Loads a large alignment file (> 1 billion rows) into DuckDB, performs filtering, joins, and
    computes reciprocal best hits. The final results are placed in a table named 'protein_alignments'.
    """
    logger.info("Loading alignments into DuckDB and finding reciprocal best-hits.")
    logger.debug(f"Source alignment file: {alignments_file}, coverage filter: {min_cov}")

    # 1) raw_alignments
    if table_exists(con, "raw_alignments"):
        # Compare row count to the file's line count
        expected_num_rows = count_lines(alignments_file)
        if table_has_n_rows(con, "raw_alignments", expected_num_rows):
            logger.debug("Raw alignments table already exists with matching row count; skipping re-creation.")
        else:
            logger.debug("Raw alignments table exists but row count differs; dropping and re-creating.")
            con.execute("DROP TABLE IF EXISTS raw_alignments")
            con.execute("""
            CREATE TABLE raw_alignments (
                query_protein_id UINTEGER,
                target_protein_id UINTEGER,
                fident FLOAT,
                alnlen USMALLINT,
                qlen USMALLINT,
                tlen USMALLINT,
                qcov FLOAT,
                tcov FLOAT,
                bits FLOAT
            )
            """)
            copy_sql = f"""
            COPY raw_alignments FROM '{alignments_file}'
            (DELIM '\t', HEADER FALSE);
            """
            con.execute(copy_sql)
            logger.debug("Raw alignments loaded into DuckDB (table: raw_alignments).")
    else:
        logger.debug("Raw alignments table does not exist; creating it now.")
        con.execute("""
        CREATE TABLE raw_alignments (
            query_protein_id UINTEGER,
            target_protein_id UINTEGER,
            fident FLOAT,
            alnlen USMALLINT,
            qlen USMALLINT,
            tlen USMALLINT,
            qcov FLOAT,
            tcov FLOAT,
            bits FLOAT
        )
        """)
        copy_sql = f"""
        COPY raw_alignments FROM '{alignments_file}'
        (DELIM '\t', HEADER FALSE);
        """
        con.execute(copy_sql)
        logger.debug("Raw alignments loaded into DuckDB (table: raw_alignments).")

    # 2) alignments_filtered
    if table_exists(con, "alignments_filtered"):
        logger.debug("alignments_filtered table already exists; skipping creation.")
    else:
        logger.debug("Filtered alignments table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS alignments_filtered")
        con.execute(f"""
        CREATE TABLE alignments_filtered AS
        SELECT *
        FROM raw_alignments
        WHERE qcov >= {min_cov}
          AND tcov >= {min_cov};
        """)
        logger.debug("Alignment table filtered by coverage (table: alignments_filtered).")

    # 3) alignments_joined
    if table_exists(con, "alignments_joined"):
        logger.debug("alignments_joined table already exists; skipping creation.")
    else:
        logger.debug("Joined alignments table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS alignments_joined")
        con.execute(f"""
        CREATE TABLE alignments_joined AS
        SELECT
            a.query_protein_id,
            a.target_protein_id,
            a.fident,
            a.alnlen,
            a.qlen,
            a.tlen,
            a.qcov,
            a.tcov,
            a.bits,
            gm1.genome_id AS query_genome_id,
            gm2.genome_id AS target_genome_id
        FROM alignments_filtered a
        JOIN {gene_map_table} gm1 ON a.query_protein_id = gm1.protein_id
        JOIN {gene_map_table} gm2 ON a.target_protein_id = gm2.protein_id
        WHERE gm1.genome_id <> gm2.genome_id
        """)
        logger.debug("Joined alignments with gene_map on query_protein_id and target_protein_id (table: alignments_joined).")

    # 4) alignments_grouped
    if table_exists(con, "alignments_grouped"):
        logger.debug("alignments_grouped table already exists; skipping creation.")
    else:
        logger.debug("Grouped alignments table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS alignments_grouped")
        con.execute("""
        CREATE TABLE alignments_grouped AS
        SELECT
            LEAST(query_genome_id, target_genome_id)   AS left_genome_id,
            GREATEST(query_genome_id, target_genome_id) AS right_genome_id,
            LEAST(query_protein_id, target_protein_id)  AS left_protein_id,
            GREATEST(query_protein_id, target_protein_id) AS right_protein_id,
            AVG(fident)         AS fident,
            CAST(AVG(alnlen) AS INT) AS alnlen,
            FIRST(qlen)         AS qlen,
            FIRST(tlen)         AS tlen,
            FIRST(qcov)         AS qcov,
            FIRST(tcov)         AS tcov,
            AVG(bits)           AS bits
        FROM alignments_joined
        GROUP BY 1,2,3,4
        """)
        logger.debug("Grouped alignments to merge duplicates for A<->B vs B<->A (table: alignments_grouped).")

    # 5) GET RECIPROCAL BEST HITS
    # for each pair of genomes, only want to keep the most significant matches per proteins
    # ie suppose we have genomes A and B with 3 ptns each
    # we could have ptn A1 with matches to B1 AND B2
    # but suppose that A1-B2 is the more significant match -- that's the one we want to keep
    # basically AAI is only calculated for the most significant RECIPROCAL matches
    # ie ptn A1 is most related to B2 out of all 3 ptns in genome B AND
    # ptn B2 is most related to A1 out of all 3 ptns in genome A
    
    # 5.1) best_by_left
    if table_exists(con, "best_by_left"):
        logger.debug("best_by_left table already exists; skipping creation.")
    else:
        logger.debug("Left best hits table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS best_by_left")
        con.execute("""
        CREATE TABLE best_by_left AS
        WITH cte AS (
          SELECT
            *,
            ROW_NUMBER() OVER (
              PARTITION BY left_genome_id, right_genome_id, left_protein_id
              ORDER BY bits DESC
            ) AS rn
          FROM alignments_grouped
        )
        SELECT *
        FROM cte
        WHERE rn = 1
        """)
        logger.debug("Chose the best alignment (highest bits) per left_protein_id (table: best_by_left).")

    # 5.2) best_by_right
    if table_exists(con, "best_by_right"):
        logger.debug("best_by_right table already exists; skipping creation.")
    else:
        logger.debug("Right best hits table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS best_by_right")
        con.execute("""
        CREATE TABLE best_by_right AS
        WITH cte AS (
          SELECT
            *,
            ROW_NUMBER() OVER (
              PARTITION BY left_genome_id, right_genome_id, right_protein_id
              ORDER BY bits DESC
            ) AS rn
          FROM alignments_grouped
        )
        SELECT *
        FROM cte
        WHERE rn = 1
        """)
        logger.debug("Chose the best alignment (highest bits) per right_protein_id (table: best_by_right).")

    # 5.3) reciprocal_best
    if table_exists(con, "reciprocal_best"):
        logger.debug("reciprocal_best table already exists; skipping creation.")
    else:
        logger.debug("Reciprocal best hits table does not exist; creating it now.")
        con.execute("DROP TABLE IF EXISTS reciprocal_best")
        con.execute("""
        CREATE TABLE reciprocal_best AS
        SELECT bleft.*
        FROM best_by_left bleft
        JOIN best_by_right bright
          ON bleft.left_genome_id   = bright.left_genome_id
         AND bleft.right_genome_id  = bright.right_genome_id
         AND bleft.left_protein_id  = bright.left_protein_id
         AND bleft.right_protein_id = bright.right_protein_id
        """)
        logger.debug("Retained only reciprocal best hits (table: reciprocal_best).")

    # 6) protein_alignments
    ## Check if 'reciprocal_best' exists & row-count matches
    if table_exists(con, "reciprocal_best"):
        reciprocal_best_count = con.execute("SELECT COUNT(*) FROM reciprocal_best").fetchone()[0]
    else:
        reciprocal_best_count = None # force creation anyway

    if table_exists(con, "protein_alignments"):
        if reciprocal_best_count is not None and table_has_n_rows(con, "protein_alignments", reciprocal_best_count):
            logger.debug("protein_alignments table exists and row count matches reciprocal_best; skipping creation.")
        else:
            logger.debug("protein_alignments table exists but row count mismatch or no reciprocal_best; re-creating.")
            con.execute("DROP TABLE IF EXISTS protein_alignments")
            con.execute("""
            CREATE TABLE protein_alignments AS
            SELECT
                left_genome_id   AS query_genome_id,
                right_genome_id  AS target_genome_id,
                left_protein_id  AS query_protein_id,
                right_protein_id AS target_protein_id,
                fident,
                alnlen,
                qlen,
                tlen,
                qcov,
                tcov,
                bits
            FROM reciprocal_best
            """)
            logger.debug("Created final table 'protein_alignments'.")
    else:
        logger.debug("Final protein alignments table does not exist; creating it now.")
        con.execute("""
        CREATE TABLE protein_alignments AS
        SELECT
            left_genome_id   AS query_genome_id,
            right_genome_id  AS target_genome_id,
            left_protein_id  AS query_protein_id,
            right_protein_id AS target_protein_id,
            fident,
            alnlen,
            qlen,
            tlen,
            qcov,
            tcov,
            bits
        FROM reciprocal_best
        """)
        logger.debug("Created final table in DuckDB 'protein_alignments'.")

def compute_aai_in_duckdb(con: duckdb.DuckDBPyConnection):
    """
    Compute AAI inside DuckDB without ever loading the entire protein_alignments
    table into Polars. 
    This creates a new table final_aai in DuckDB.
    """
    logger.info("Computing AAI in DuckDB from the protein_alignments table.")
    
    # 1) Create a small helper table: genome_info
    logger.debug("Creating genome_info table.")
    con.execute("DROP TABLE IF EXISTS genome_info")
    con.execute("""
    CREATE TABLE genome_info AS
    SELECT
        genome_id,
        COUNT(DISTINCT protein_id) AS num_proteins
    FROM gene_map
    GROUP BY genome_id
    """)

    # 2) Create a final table with all aggregated results:
    logger.debug("Creating final_aai table.")
    con.execute("DROP TABLE IF EXISTS final_aai")
    con.execute("""
    CREATE TABLE final_aai AS
    SELECT 
        pa.query_genome_id                       AS query_genome_id,
        pa.target_genome_id                      AS target_genome_id,
        AVG(pa.fident)                           AS aai,
        COUNT(*)                                 AS shared_genes,
        MAX(gi1.num_proteins)                    AS query_n_ptns,
        MAX(gi2.num_proteins)                    AS target_n_ptns,
        (COUNT(*) * 1.0 / MAX(gi1.num_proteins)) AS query_shared,
        (COUNT(*) * 1.0 / MAX(gi2.num_proteins)) AS target_shared
    FROM protein_alignments pa
    JOIN genome_info gi1
    ON pa.query_genome_id = gi1.genome_id
    JOIN genome_info gi2
    ON pa.target_genome_id = gi2.genome_id
    GROUP BY pa.query_genome_id, pa.target_genome_id
    ORDER BY pa.query_genome_id, pa.target_genome_id
    """)
    
    logger.info("AAI computation completed in DuckDB.")

def check_unaccounted_genomes(con: duckdb.DuckDBPyConnection, gene_map_df: pl.DataFrame, db_path: str):
    """
    Identifies any genome_ids present in gene_map_df but missing from final_aai,
    i.e., genome_id does not appear in final_aai.query_genome_id or .target_genome_id.
    Writes them to a file if any are found. Avoids correlated subqueries for better performance.
    """
    logger.info("Checking for genomes not present in the final AAI table")

    total_genomes = gene_map_df["genome_id"].n_unique()
    logger.debug(f"Total input genomes from gene_map: {total_genomes:,}")

    # 1) Create a table 'used_genomes' that is the set of genome_id from BOTH query_genome_id and target_genome_id
    con.execute("DROP TABLE IF EXISTS used_genomes")
    con.execute("""
    CREATE TEMP TABLE used_genomes AS
    SELECT query_genome_id AS genome_id
    FROM final_aai
    UNION
    SELECT target_genome_id AS genome_id
    FROM final_aai
    """)

    # 2) Create a table 'all_genomes' from gene_map (only the distinct genome_id)
    con.execute("DROP TABLE IF EXISTS all_genomes")
    con.execute("""
    CREATE TEMP TABLE all_genomes AS
    SELECT DISTINCT genome_id
    FROM gene_map
    """)

    # 3) Now do a LEFT JOIN to find which genome_ids from all_genomes are NOT in used_genomes
    con.execute("DROP TABLE IF EXISTS unaccounted_genomes")
    con.execute("""
    CREATE TABLE unaccounted_genomes AS
    SELECT a.genome_id
    FROM all_genomes a
    LEFT JOIN used_genomes u
         ON a.genome_id = u.genome_id
    WHERE u.genome_id IS NULL
    """)

    # 4) Count how many were missing
    unaccounted_count = con.execute("SELECT COUNT(*) FROM unaccounted_genomes").fetchone()[0]
    accounted_count = total_genomes - unaccounted_count
    logger.debug(f"accounted: {accounted_count:,}, unaccounted: {unaccounted_count:,}")

    if unaccounted_count > 0:
        percent_unaccounted = round(unaccounted_count / total_genomes * 100, 2)
        logger.warning(
            f"There are {unaccounted_count:,} genomes present in the gene-to-genome "
            f"mapping file that are not in the AAI results ({percent_unaccounted}%). "
            "This is likely due to the previous mmseqs search sensitivity, the "
            "minimum coverage parameter for mmseqs search, and/or a large number of "
            "proteins present in only one input genome. This should be fine, but "
            "the unaccounted genomes will end up as unclustered singletons in the "
            "genome clustering step."
        )
        unaccounted_file = os.path.join(os.path.dirname(db_path), "genomes_not_in_aai_results.txt")
        logger.warning(f"Writing IDs of unaccounted genomes to {unaccounted_file}")

        # 5) Use DuckDB's COPY to write all unaccounted IDs to disk without loading them into Python memory
        con.execute(f"COPY unaccounted_genomes TO '{unaccounted_file}' (DELIM '\t', HEADER false)")
        logger.debug("Wrote unaccounted genome IDs to file via DuckDB's COPY command.")

        # Drop unaccounted_genomes
        con.execute("DROP TABLE unaccounted_genomes")
    else:
        logger.info("No unaccounted genomes found. All input genomes appear in final_aai.")

    # Cleanup temp tables
    con.execute("DROP TABLE all_genomes")
    con.execute("DROP TABLE used_genomes")

def cleanup(con: duckdb.DuckDBPyConnection, tmp_spill_path: Path, alignments_file: str) -> None:
    """
    Drosp DuckDB temporary tables that are no longer needed to save disk space.
    This should only occur if the key tables gene_map, protein_alignments, and final_aai exist 
    and are not empty. If any of these tables are missing or empty, this should not drop the 
    temporary tables, nor proceed with removing any other files, because something went wrong.
    """
    key_tables = ["gene_map", "protein_alignments", "final_aai"]
    drop_temp_tables = True
    for tbl in key_tables:
        if not table_exists(con, tbl):
            drop_temp_tables = False
            logger.error(f"The needed final table {tbl} does not exist after AAI computation; exiting.")
            raise ValueError(f"The needed final table {tbl} does not exist after AAI computation.")
        else:
            count = con.execute(f"SELECT COUNT(*) FROM {tbl}").fetchone()[0]
            if count == 0:
                drop_temp_tables = False
                logger.error(f"The needed final table {tbl} is empty after AAI computation; exiting.")
                raise ValueError(f"The needed final table {tbl} is empty after AAI computation.")
    if drop_temp_tables:
        con.execute("DROP TABLE IF EXISTS raw_alignments")
        con.execute("DROP TABLE IF EXISTS alignments_filtered")
        con.execute("DROP TABLE IF EXISTS alignments_joined")
        con.execute("DROP TABLE IF EXISTS alignments_grouped")
        con.execute("DROP TABLE IF EXISTS best_by_left")
        con.execute("DROP TABLE IF EXISTS best_by_right")
        con.execute("DROP TABLE IF EXISTS reciprocal_best")
        con.execute("DROP TABLE IF EXISTS genome_info")
        logger.debug("Dropped temporary tables: raw_alignments, alignments_filtered, alignments_joined, alignments_grouped, best_by_left, best_by_right, reciprocal_best, and genome_info.")
    else:
        logger.debug("Not dropping temporary tables as one or more key tables are missing or empty.")

    # Remove the DuckDB temporary spill directory to save space
    try:
        if os.path.exists(tmp_spill_path):
            shutil.rmtree(tmp_spill_path)
            logger.debug(f"Removed temporary spill directory: {tmp_spill_path}")
    except Exception as e:
        logger.warning(f"Failed to remove spill directory: {e}")
    
    # Remove the alignment file directory to save space
    try:
        alignments_dir = os.path.dirname(alignments_file)
        if os.path.exists(alignments_dir):
            shutil.rmtree(alignments_dir)
            logger.debug(f"Removed temporary alignments directory: {alignments_dir}")
    except Exception as e:
        logger.warning(f"Failed to remove alignments directory: {e}")
        
def main():
    alignments_file = snakemake.params.search_output
    gene_map_file = snakemake.params.gene_map_file
    db_path = os.path.join(snakemake.params.wdir, "alignments_and_clusters.duckdb")
    min_cov = snakemake.params.cov_fraction
    mem_limit = snakemake.resources.mem
    num_cpus = snakemake.threads
    os.environ["NUMEXPR_MAX_THREADS"] = str(num_cpus)

    logger.info("Computing pairwise amino-acid identities of all genomes...")
    logger.debug(f"Memory limit set to {mem_limit} GB")
    set_memory_limit(mem_limit)

    # Connect to DuckDB
    tmp_spill_path = Path(str(db_path) + "_spill_tmp")
    con = connect_to_duckdb(db_path, tmp_spill_path, num_cpus, mem_limit)

    # Read gene map in Polars
    gene_map_df = read_gene_map(gene_map_file)
    
    # Check if if the final_aai table already exists
    ## This is in case the pipeline is re-run after the final_aai table
    ## was already created and the temporary tables were dropped, but 
    ## something went wrong after that.
    if table_exists(con, "final_aai") and con.execute(f"SELECT COUNT(*) FROM final_aai").fetchone()[0] > 0:
        logger.info("final_aai table exists and has rows; skipping AAI computation.")
    else:
        logger.debug("No AAI results found in DuckDB; proceeding with workflow.")

        # Check if gene_map table exists in DuckDB and if row count matches
        check_gene_map_table(con, gene_map_df)

        # Now read alignments into DuckDB, filter by coverage, and compute reciprocal best hits
        load_alignments_in_duckdb(alignments_file, con, "gene_map", min_cov)
        logger.info(f"Alignments table successfully read, filtered, and processed with DuckDB.")
        
        # Compute AAI in DuckDB
        compute_aai_in_duckdb(con)

    # Check for genomes in the gene map that are not present in the final_aai table
    check_unaccounted_genomes(con, gene_map_df, db_path)
    
    # Clean up and remove temporary files
    cleanup(con, tmp_spill_path, alignments_file)

    logger.info("AAI computation and processing of results completed.")

if __name__ == "__main__":
    main()
