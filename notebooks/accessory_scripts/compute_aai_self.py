import os
import psutil
import gc
import logging
import argparse
import resource
from pathlib import Path
from typing import Union

import duckdb
from pyfastatools import Parser
import polars as pl

FilePath = Union[str, Path]

logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
logger = logging.getLogger(__name__)


def set_memory_limit(limit_in_gb: int) -> None:
    if not limit_in_gb or limit_in_gb <= 0:
        return
    limit_in_bytes = int(limit_in_gb) * 1024**3
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
        logger.info(f"Memory limit set to {limit_in_gb} GB")
    except (ValueError, OSError, AttributeError):
        pass


def genome_sizes(fasta: FilePath) -> pl.DataFrame:
    g: dict[str, int] = {}
    for n, rec in enumerate(Parser(str(fasta)), 1):
        genome = "_".join(rec.header.name.split("_")[:-1])
        g[genome] = g.get(genome, 0) + 1
        if n % 5_000_000 == 0:
            logger.info(f"Scanned {n:,} proteins")
    return pl.DataFrame({"genome": list(g.keys()), "n_ptns": list(g.values())})


def write_aai(
    con: duckdb.DuckDBPyConnection,
    alignments_file: FilePath,
    sizes: pl.DataFrame,
    out_path: FilePath,
) -> None:
    # Register genome sizes as a DuckDB view from the Polars DataFrame
    con.register("genome_sizes", sizes.to_arrow())

    aln_path = str(alignments_file)
    out_path = str(out_path)

    # Detect TSV vs Parquet
    suffix = Path(alignments_file).suffix.lower()
    if suffix == ".parquet":
        aln_source = f"read_parquet('{aln_path}')"
        # Parquet already has named columns; project only what is needed
        col_select = """
            query::VARCHAR                                        AS query,
            target::VARCHAR                                       AS target,
            fident::FLOAT                                         AS fident,
            bits::FLOAT                                           AS bits,
            regexp_replace(query,  '_[^_]+$', '')::VARCHAR        AS query_genome,
            regexp_replace(target, '_[^_]+$', '')::VARCHAR        AS target_genome
        """
    else:
        # TSV with no header, name columns positionally
        aln_source = f"""read_csv('{aln_path}', delim='\t', header=false,
            columns={{
                'query':'VARCHAR','target':'VARCHAR','fident':'FLOAT','alnlen':'INTEGER',
                'mismatch':'INTEGER','gapopen':'INTEGER','qstart':'INTEGER','qend':'INTEGER',
                'tstart':'INTEGER','tend':'INTEGER','evalue':'DOUBLE','bits':'FLOAT',
                'qlen':'INTEGER','tlen':'INTEGER','qcov':'FLOAT','tcov':'FLOAT'
            }})"""
        col_select = """
            query, target, fident, bits,
            regexp_replace(query,  '_[^_]+$', '') AS query_genome,
            regexp_replace(target, '_[^_]+$', '') AS target_genome
        """

    sql = f"""
    COPY (
        WITH

        -- 1. Strip self-genome hits and project only needed columns
        alns AS (
            SELECT {col_select}
            FROM {aln_source}
            WHERE regexp_replace(query,  '_[^_]+$', '')
               != regexp_replace(target, '_[^_]+$', '')
        ),

        -- 2. Best hit per (query protein -> target genome)
        ranked_q AS (
            SELECT *,
                   RANK() OVER (PARTITION BY query, target_genome ORDER BY bits DESC) AS rank_q
            FROM alns
        ),
        best_q AS (SELECT * FROM ranked_q WHERE rank_q = 1),

        -- 3. Reciprocal: best hit per (target protein -> query genome)
        ranked_t AS (
            SELECT *,
                   RANK() OVER (PARTITION BY target, query_genome ORDER BY bits DESC) AS rank_t
            FROM best_q
        ),
        rbh AS (
            SELECT query_genome, target_genome, fident
            FROM ranked_t
            WHERE rank_t = 1
        ),

        -- 4. Aggregate RBH pairs
        aai AS (
            SELECT
                query_genome,
                target_genome,
                AVG(fident)::FLOAT   AS aai,
                COUNT(*)::INTEGER    AS shared_genes
            FROM rbh
            GROUP BY query_genome, target_genome
        )

        -- 5. Attach genome sizes and compute shared fractions
        SELECT
            a.query_genome,
            a.target_genome,
            a.aai,
            a.shared_genes,
            (a.shared_genes / q.n_ptns)::FLOAT  AS query_shared,
            (a.shared_genes / t.n_ptns)::FLOAT  AS target_shared
        FROM aai a
        JOIN genome_sizes q ON a.query_genome = q.genome
        JOIN genome_sizes t ON a.target_genome = t.genome

    ) TO '{out_path}' (FORMAT PARQUET);
    """

    logger.info("Running AAI query...")
    con.execute(sql)
    logger.info(f"AAI written to {out_path}")


def main(ptn_fa: FilePath, alignments_file: FilePath, out_dir: FilePath, threads: int, mem_gb: int) -> None:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    sizes = genome_sizes(ptn_fa)
    logger.info(f"Genome sizes computed: {len(sizes):,} genomes")

    con = duckdb.connect()
    con.execute(f"PRAGMA threads={threads}")
    if mem_gb > 0:
        con.execute(f"PRAGMA memory_limit='{mem_gb*0.67}GB'")
        # Spill directory, DuckDB will write temp files here if it exceeds the limit
        con.execute(f"PRAGMA temp_directory='{out_dir}/duckdb_tmp'")

    con.execute("SET enable_progress_bar = true;")
    con.execute("SET enable_progress_bar_print = true;")
    write_aai(con, alignments_file, sizes, out_dir / "aai.parquet")
    con.close()
    gc.collect()


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--ptn_path",        required=True)
    p.add_argument("--alignments_file", required=True)
    p.add_argument("--output_dir",      required=True)
    p.add_argument("--threads",  type=int, default=1)
    p.add_argument("--mem_limit", type=int, default=0, help="RAM budget in GB; DuckDB spills beyond this")
    args = p.parse_args()

    if args.mem_limit == 0:
        # If no memory limit specified, use 80% of available RAM
        args.mem_limit = int(psutil.virtual_memory().available * 0.8 / 1024**3)
    set_memory_limit(args.mem_limit)
    os.environ["POLARS_MAX_THREADS"] = str(args.threads)
    os.environ["NUMEXPR_MAX_THREADS"] = str(args.threads)
    import polars as pl
    import duckdb
    main(args.ptn_path, args.alignments_file, args.output_dir, args.threads, args.mem_limit)