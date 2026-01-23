#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import math
import re
from pathlib import Path

def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    except (ValueError, OSError, AttributeError) as e:
        logger.warning(f"Unable to set memory limit. Error: {e}")

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

print("========================================================================\n   Step 9/11: Curate the predicted functions based on genomic context   \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n   Step 9/11: Curate the predicted functions based on genomic context   \n========================================================================\n")

# Global caches for thresholds
KEGG_THRESHOLDS = {}
FOAM_THRESHOLDS = {}
CAMPER_THRESHOLDS = {}

def _nan_to_none(v):
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    return float(v)

# Load KEGG thresholds
KEGG_THRESHOLDS = {}
KEGG_THRESHOLDS_PATH = snakemake.params.kegg_cutoff_file
if Path(KEGG_THRESHOLDS_PATH).exists():
    kegg_df = pl.read_csv(
        KEGG_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "threshold": pl.Float64},
        separator="\t",
    )
    # Treat null/NaN thresholds as missing to fall back to bypass_min_bitscore downstream
    ids = kegg_df["id"].to_list()
    thrs = [_nan_to_none(v) for v in kegg_df["threshold"].to_list()]
    KEGG_THRESHOLDS = {i: t for i, t in zip(ids, thrs) if t is not None}
else:
    logger.warning(f"KEGG thresholds file not found at {KEGG_THRESHOLDS_PATH}. KEGG thresholds will not be used to filter HMMsearch results!")

# Load FOAM thresholds
FOAM_THRESHOLDS = {}
FOAM_THRESHOLDS_PATH = snakemake.params.foam_cutoff_file
if Path(FOAM_THRESHOLDS_PATH).exists():
    foam_df = pl.read_csv(
        FOAM_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "cutoff_full": pl.Float64, "cutoff_domain": pl.Float64},
        separator="\t",
    )
    ids = foam_df["id"].to_list()
    fulls = foam_df["cutoff_full"].to_list()
    doms = foam_df["cutoff_domain"].to_list()
    FOAM_THRESHOLDS = {i: (_nan_to_none(f), _nan_to_none(d)) for i, f, d in zip(ids, fulls, doms)}
else:
    logger.warning(f"FOAM thresholds file not found at {FOAM_THRESHOLDS_PATH}. FOAM thresholds will not be used to filter HMMsearch results!")

# Load CAMPER thresholds
CAMPER_THRESHOLDS = {}
CAMPER_THRESHOLDS_PATH = snakemake.params.camper_cutoff_file
if Path(CAMPER_THRESHOLDS_PATH).exists():
    camper_df = pl.read_csv(
        CAMPER_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "cutoff_full": pl.Float64, "cutoff_domain": pl.Float64},
        separator="\t",
    )
    ids = camper_df["id"].to_list()
    fulls = camper_df["cutoff_full"].to_list()
    doms = camper_df["cutoff_domain"].to_list()
    CAMPER_THRESHOLDS = {i: (_nan_to_none(f), _nan_to_none(d)) for i, f, d in zip(ids, fulls, doms)}
else:
    logger.warning(f"CAMPER thresholds file not found at {CAMPER_THRESHOLDS_PATH}. CAMPER thresholds will not be used to filter HMMsearch results!")

KEGG_THRESHOLDS_IDS = set(KEGG_THRESHOLDS.keys())
FOAM_THRESHOLDS_IDS = set(FOAM_THRESHOLDS.keys())
CAMPER_THRESHOLDS_IDS = set(CAMPER_THRESHOLDS.keys())

def dedup_desc_by_id(df, id_col="id", db=None):
    logger.debug(f"Deduplicating descriptions for {db} HMMs by id column: {id_col}")
    logger.debug(f"Before deduplication: {df}")
    # Deterministic 1 row per id: prefer longer non-null name, tie-break by name lexicographically
    df_dedup = (
        df.with_columns(
            pl.col("name").cast(pl.Utf8).fill_null("").alias("name"),
            pl.col("name").cast(pl.Utf8).fill_null("").str.len_bytes().alias("_name_len"),
        )
        .sort([id_col, "_name_len", "name"], descending=[False, True, False])
        .unique(subset=[id_col], keep="first")
        .drop("_name_len")
    )
    logger.debug(f"After deduplication: {df_dedup}")
    return df_dedup

def dedup_desc_by_id_norm(df, norm_col="id_norm", db=None):
    logger.debug(f"Deduplicating descriptions for {db} HMMs by normalized id column: {norm_col}")
    logger.debug(f"Before deduplication: {df}")
    df_dedup = (
        df.with_columns(
            pl.col("name").cast(pl.Utf8).fill_null("").alias("name"),
            pl.col("name").cast(pl.Utf8).fill_null("").str.len_bytes().alias("_name_len"),
        )
        .sort([norm_col, "_name_len", "name"], descending=[False, True, False])
        .unique(subset=[norm_col], keep="first")
        .drop("_name_len")
    )
    logger.debug(f"After deduplication: {df_dedup}")
    return df_dedup

def summarize_annot_table(table, hmm_descriptions):
    """
    Summarize per-gene annotations and attach reference descriptions.

    Inputs:
    - table (pl.DataFrame): gene context table containing per-database HMM IDs/scores/coverage and
    flanking-distance columns (some may be missing)
    - hmm_descriptions (pl.DataFrame): reference HMM descriptions table with columns ["id", "db", "name"]

    Behavior:
    - Ensures required annotation columns exist (fills missing with nulls)
    - Computes min left/right distances to viral and MGE flanking genes
    - Joins in per-database descriptions (with simple ID normalization for Pfam/dbCAN)
    - Computes top_hit_hmm_id / top_hit_description / top_hit_db from the max score across databases
    - De-duplicates and sorts by genomic position

    Returns:
    - pl.DataFrame: curated annotation table (one row per gene) with standardized columns, joined
    descriptions, top-hit fields, and flanking-distance context.
    """
    logger.debug(f"Initial table columns: {table.columns}")

    # Normalize coverage fields: map *_coverage to the *_coverage_hmm (preferred) or *_coverage_sequence columns if present.
    coverage_map = {
        "KEGG": ("KEGG_coverage_hmm", "KEGG_coverage_sequence"),
        "FOAM": ("FOAM_coverage_hmm", "FOAM_coverage_sequence"),
        "Pfam": ("Pfam_coverage_hmm", "Pfam_coverage_sequence"),
        "dbCAN": ("dbCAN_coverage_hmm", "dbCAN_coverage_sequence"),
        "METABOLIC": ("METABOLIC_coverage_hmm", "METABOLIC_coverage_sequence"),
        "CAMPER": ("CAMPER_coverage_hmm", "CAMPER_coverage_sequence"),
        "PHROG": ("PHROG_coverage_hmm", "PHROG_coverage_sequence"),
    }
    for src, (hmm_col, seq_col) in coverage_map.items():
        out_col = f"{src}_coverage"
        if out_col not in table.columns:
            if hmm_col in table.columns:
                table = table.with_columns(pl.col(hmm_col).cast(pl.Float64).alias(out_col))
            elif seq_col in table.columns:
                table = table.with_columns(pl.col(seq_col).cast(pl.Float64).alias(out_col))
            else:
                table = table.with_columns(pl.lit(None, dtype=pl.Float64).alias(out_col))

    window_map = {
        "window_avg_KEGG_VL-score_viral": "window_avg_KEGG_VL-score",
        "window_avg_Pfam_VL-score_viral": "window_avg_Pfam_VL-score",
        "window_avg_PHROG_VL-score_viral": "window_avg_PHROG_VL-score",
    }
    for out_col, in_col in window_map.items():
        if out_col not in table.columns:
            if in_col in table.columns:
                table = table.with_columns(pl.col(in_col).cast(pl.Float64).alias(out_col))
            else:
                table = table.with_columns(pl.lit(None, dtype=pl.Float64).alias(out_col))

    table = table.with_columns(
        [
            pl.min_horizontal(
                [pl.col("KEGG_viral_left_dist"), pl.col("Pfam_viral_left_dist"), pl.col("PHROG_viral_left_dist")]
            ).alias("Viral_Flanking_Genes_Left_Dist"),
            pl.min_horizontal(
                [pl.col("KEGG_viral_right_dist"), pl.col("Pfam_viral_right_dist"), pl.col("PHROG_viral_right_dist")]
            ).alias("Viral_Flanking_Genes_Right_Dist"),
            pl.min_horizontal(
                [pl.col("KEGG_MGE_left_dist"), pl.col("Pfam_MGE_left_dist"), pl.col("PHROG_MGE_left_dist")]
            ).alias("MGE_Flanking_Genes_Left_Dist"),
            pl.min_horizontal(
                [pl.col("KEGG_MGE_right_dist"), pl.col("Pfam_MGE_right_dist"), pl.col("PHROG_MGE_right_dist")]
            ).alias("MGE_Flanking_Genes_Right_Dist"),
        ]
    )

    required_cols = [
        "protein",
        "contig",
        "circular_contig",
        "genome",
        "gene_number",
        "KEGG_hmm_id",
        "FOAM_hmm_id",
        "Pfam_hmm_id",
        "dbCAN_hmm_id",
        "METABOLIC_hmm_id",
        "CAMPER_hmm_id",
        "PHROG_hmm_id",
        "KEGG_score",
        "FOAM_score",
        "Pfam_score",
        "dbCAN_score",
        "METABOLIC_score",
        "CAMPER_score",
        "PHROG_score",
        "KEGG_coverage",
        "FOAM_coverage",
        "Pfam_coverage",
        "dbCAN_coverage",
        "METABOLIC_coverage",
        "CAMPER_coverage",
        "PHROG_coverage",
        "KEGG_V-score",
        "Pfam_V-score",
        "PHROG_V-score",
        "window_avg_KEGG_VL-score_viral",
        "window_avg_Pfam_VL-score_viral",
        "window_avg_PHROG_VL-score_viral",
        "Viral_Flanking_Genes_Left_Dist",
        "Viral_Flanking_Genes_Right_Dist",
        "MGE_Flanking_Genes_Left_Dist",
        "MGE_Flanking_Genes_Right_Dist",
        "Viral_Origin_Confidence",
    ]

    for col in required_cols:
        if col not in table.columns:
            if col.endswith("_id") or col.endswith("_Description") or col in ["protein", "contig", "genome", "circular_contig", "Viral_Origin_Confidence"]:
                dtype = pl.Utf8
            elif col.endswith("_score") or col.endswith("_coverage") or col.endswith("_V-score") or col.endswith("_VL-score"):
                dtype = pl.Float64
            else:
                dtype = pl.Utf8
            table = table.with_columns(pl.lit(None, dtype=dtype).alias(col))

    table = table.select(required_cols).rename({"protein": "Protein"})
    logger.debug(f"Table columns after ensuring required columns: {table.columns}")
    logger.debug(f"Number of rows before description joins: {table.height:,}")
    logger.debug(f"Number of unique proteins before description joins: {table.select(pl.col('Protein').unique()).height:,}")
    logger.debug(f"Number of unique contigs before description joins: {table.select(pl.col('contig').unique()).height:,}")
    logger.debug(f"Number of unique genomes before description joins: {table.select(pl.col('genome').unique()).height:,}")

    # Build per-db, de-duplicated description tables
    hmm_kegg = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "KEGG").select(["id", "name"]), db="KEGG")
    hmm_foam = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "FOAM").select(["id", "name"]), db="FOAM")
    hmm_dbcan = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "dbCAN").select(["id", "name"]), db="dbCAN")
    hmm_metabolic = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "METABOLIC").select(["id", "name"]), db="METABOLIC")
    hmm_camper = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "CAMPER").select(["id", "name"]), db="CAMPER")
    hmm_phrog = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "PHROG").select(["id", "name"]), db="PHROG")

    # Pfam needs id normalization; dedupe on id_norm
    hmm_pfam = (
        hmm_descriptions.filter(pl.col("db") == "Pfam")
        .with_columns(pl.col("id").str.replace(r"\.\d+$", "", literal=False).alias("id_norm"))
        .select(["id_norm", "name"])
    )
    hmm_pfam = dedup_desc_by_id_norm(hmm_pfam, norm_col="id_norm", db="Pfam")

    def _join_and_rename(table, hmm_df, left_on, right_on, db_name):
        logger.debug(f"Joining {db_name} descriptions")
        logger.debug(f"Table columns before {db_name} join: {table.columns}")
        logger.debug(f"hmm_{db_name.lower()} columns: {hmm_df.columns}")
        table = table.join(hmm_df, left_on=left_on, right_on=right_on, how="left").rename({"name": f"{db_name}_Description"})
        if right_on in table.columns:
            table = table.drop(right_on)
        logger.debug(f"Table columns after {db_name} join: {table.columns}")
        return table

    # KEGG join (db-filtered, deduped)
    table = _join_and_rename(table, hmm_kegg, "KEGG_hmm_id", "id", "KEGG")

    # FOAM join (db-filtered, deduped)
    table = _join_and_rename(table, hmm_foam, "FOAM_hmm_id", "id", "FOAM")

    # Pfam join (normalized, db-filtered, deduped)
    table = table.with_columns(pl.col("Pfam_hmm_id").str.replace(r"\.\d+$", "", literal=False).alias("Pfam_hmm_id_norm"))
    table = _join_and_rename(table, hmm_pfam, "Pfam_hmm_id_norm", "id_norm", "Pfam")
    if "Pfam_hmm_id_norm" in table.columns:
        table = table.drop("Pfam_hmm_id_norm")
    logger.debug(f"Table columns after dropping columns: {table.columns}")

    # dbCAN join (underscore normalization + db-filtered, deduped)
    table = table.with_columns(pl.col("dbCAN_hmm_id").str.replace(r"_(.*)", "", literal=False).alias("dbCAN_hmm_id_no_underscore"))
    table = _join_and_rename(table, hmm_dbcan, "dbCAN_hmm_id_no_underscore", "id", "dbCAN")
    if "dbCAN_hmm_id_no_underscore" in table.columns:
        table = table.drop("dbCAN_hmm_id_no_underscore")
    logger.debug(f"Table columns after dropping columns: {table.columns}")

    # METABOLIC join (db-filtered, deduped)
    table = _join_and_rename(table, hmm_metabolic, "METABOLIC_hmm_id", "id", "METABOLIC")

    # CAMPER join (db-filtered, deduped)
    table = _join_and_rename(table, hmm_camper, "CAMPER_hmm_id", "id", "CAMPER")

    # PHROG join (db-filtered, deduped)
    table = _join_and_rename(table, hmm_phrog, "PHROG_hmm_id", "id", "PHROG")

    # Top-hit selection
    logger.debug("Computing top hits across databases")
    score_cols = ["KEGG_score", "FOAM_score", "Pfam_score", "dbCAN_score", "METABOLIC_score", "CAMPER_score", "PHROG_score"]
    table = table.with_columns([pl.col(c).cast(pl.Float64).fill_null(float("-inf")).alias(c) for c in score_cols])

    table = table.with_columns(pl.max_horizontal(score_cols).alias("max_score"))
    table = table.with_columns(
        pl.when(pl.col("max_score").is_null())
        .then(None)
        .otherwise(
            pl.struct(score_cols).map_elements(
                lambda row: list(row.values()).index(max(row.values())),
                return_dtype=pl.Int64,
            )
        )
        .alias("best_idx")
    ).drop("max_score")

    table = table.with_columns(
        [
            pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_hmm_id"))
            .when(pl.col("best_idx") == 1).then(pl.col("FOAM_hmm_id"))
            .when(pl.col("best_idx") == 2).then(pl.col("Pfam_hmm_id"))
            .when(pl.col("best_idx") == 3).then(pl.col("dbCAN_hmm_id"))
            .when(pl.col("best_idx") == 4).then(pl.col("METABOLIC_hmm_id"))
            .when(pl.col("best_idx") == 5).then(pl.col("CAMPER_hmm_id"))
            .when(pl.col("best_idx") == 6).then(pl.col("PHROG_hmm_id"))
            .otherwise(pl.lit(None))
            .alias("top_hit_hmm_id"),
            pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_Description"))
            .when(pl.col("best_idx") == 1).then(pl.col("FOAM_Description"))
            .when(pl.col("best_idx") == 2).then(pl.col("Pfam_Description"))
            .when(pl.col("best_idx") == 3).then(pl.col("dbCAN_Description"))
            .when(pl.col("best_idx") == 4).then(pl.col("METABOLIC_Description"))
            .when(pl.col("best_idx") == 5).then(pl.col("CAMPER_Description"))
            .when(pl.col("best_idx") == 6).then(pl.col("PHROG_Description"))
            .otherwise(pl.lit(None))
            .alias("top_hit_description"),
            pl.when(pl.col("best_idx") == 0).then(pl.lit("KEGG"))
            .when(pl.col("best_idx") == 1).then(pl.lit("FOAM"))
            .when(pl.col("best_idx") == 2).then(pl.lit("Pfam"))
            .when(pl.col("best_idx") == 3).then(pl.lit("dbCAN"))
            .when(pl.col("best_idx") == 4).then(pl.lit("METABOLIC"))
            .when(pl.col("best_idx") == 5).then(pl.lit("CAMPER"))
            .when(pl.col("best_idx") == 6).then(pl.lit("PHROG"))
            .otherwise(pl.lit(None))
            .alias("top_hit_db"),
        ]
    ).drop("best_idx")

    table = table.select(
        [
            "Protein",
            "contig",
            "genome",
            "gene_number",
            "KEGG_V-score",
            "Pfam_V-score",
            "PHROG_V-score",
            "KEGG_hmm_id",
            "KEGG_Description",
            "KEGG_score",
            "KEGG_coverage",
            "FOAM_hmm_id",
            "FOAM_Description",
            "FOAM_score",
            "FOAM_coverage",
            "Pfam_hmm_id",
            "Pfam_Description",
            "Pfam_score",
            "Pfam_coverage",
            "dbCAN_hmm_id",
            "dbCAN_Description",
            "dbCAN_score",
            "dbCAN_coverage",
            "METABOLIC_hmm_id",
            "METABOLIC_Description",
            "METABOLIC_score",
            "METABOLIC_coverage",
            "CAMPER_hmm_id",
            "CAMPER_Description",
            "CAMPER_score",
            "CAMPER_coverage",
            "PHROG_hmm_id",
            "PHROG_Description",
            "PHROG_score",
            "PHROG_coverage",
            "top_hit_hmm_id",
            "top_hit_description",
            "top_hit_db",
            "circular_contig",
            "Viral_Origin_Confidence",
            "Viral_Flanking_Genes_Left_Dist",
            "Viral_Flanking_Genes_Right_Dist",
            "MGE_Flanking_Genes_Left_Dist",
            "MGE_Flanking_Genes_Right_Dist",
        ]
    ).rename({"contig": "Contig", "genome": "Genome", "circular_contig": "Circular_Contig"})
    logger.debug(f"Table columns before deduplication and sorting: {table.columns}")
    logger.debug(f"Number of rows before deduplication and sorting: {table.height:,}")
    logger.debug(f"Number of unique proteins before deduplication and sorting: {table.select(pl.col('Protein').unique()).height:,}")
    logger.debug(f"Number of unique contigs before deduplication and sorting: {table.select(pl.col('Contig').unique()).height:,}")
    logger.debug(f"Number of unique genomes before deduplication and sortings: {table.select(pl.col('Genome').unique()).height:,}")

    table = (
        table
        .unique(subset=["Protein", "Genome", "Contig", "gene_number"], keep="first")
        .sort(["Genome", "Contig", "gene_number"])
    )

    logger.debug(f"Final table columns after deduplication and sorting: {table.columns}")
    logger.debug(f"Number of rows after deduplication and sorting: {table.height:,}")
    logger.debug(f"Number of unique proteins after deduplication and sorting: {table.select(pl.col('Protein').unique()).height:,}")
    logger.debug(f"Number of unique contigs after deduplication and sorting: {table.select(pl.col('Contig').unique()).height:,}")
    logger.debug(f"Number of unique genomes before deduplication and sorting: {table.select(pl.col('Genome').unique()).height:,}")
    return table

# Helpers for Function column construction
def _clean_text_expr(expr: pl.Expr) -> pl.Expr:
    # Convert empty/whitespace-only strings to null for consistent coalesce behavior
    expr = expr.cast(pl.Utf8)
    return (
        pl.when(expr.is_null())
        .then(pl.lit(None, dtype=pl.Utf8))
        .otherwise(
            pl.when(expr.str.strip_chars() == "")
            .then(pl.lit(None, dtype=pl.Utf8))
            .otherwise(expr)
        )
    )

def _best_desc_for_ref_ids(df: pl.DataFrame, ref_ids: list[str]) -> pl.Expr:
    # Pick the highest-score description among HMM hits whose IDs appear in ref_ids
    def _pfam_id_expr():
        if "Pfam_hmm_id_clean" in df.columns:
            return pl.col("Pfam_hmm_id_clean").cast(pl.Utf8)
        if "Pfam_hmm_id" in df.columns:
            return pl.col("Pfam_hmm_id").cast(pl.Utf8).str.replace(r"\.\d+$", "", literal=False)
        return pl.lit(None, dtype=pl.Utf8)

    sources = [
        ("KEGG_hmm_id", "KEGG_score", "KEGG_Description", pl.col("KEGG_hmm_id").cast(pl.Utf8) if "KEGG_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("FOAM_hmm_id", "FOAM_score", "FOAM_Description", pl.col("FOAM_hmm_id").cast(pl.Utf8) if "FOAM_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("Pfam_hmm_id", "Pfam_score", "Pfam_Description", _pfam_id_expr()),
        ("dbCAN_hmm_id", "dbCAN_score", "dbCAN_Description", pl.col("dbCAN_hmm_id").cast(pl.Utf8) if "dbCAN_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("METABOLIC_hmm_id", "METABOLIC_score", "METABOLIC_Description", pl.col("METABOLIC_hmm_id").cast(pl.Utf8) if "METABOLIC_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("CAMPER_hmm_id", "CAMPER_score", "CAMPER_Description", pl.col("CAMPER_hmm_id").cast(pl.Utf8) if "CAMPER_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("PHROG_hmm_id", "PHROG_score", "PHROG_Description", pl.col("PHROG_hmm_id").cast(pl.Utf8) if "PHROG_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
    ]

    cand_exprs = []
    desc_cols = []
    for _, score_col, desc_col, id_expr in sources:
        if score_col not in df.columns or desc_col not in df.columns:
            continue
        cond = id_expr.is_in(ref_ids).fill_null(False)
        score = pl.col(score_col).cast(pl.Float64).fill_null(float("-inf"))
        cand_exprs.append(pl.when(cond).then(score).otherwise(pl.lit(float("-inf"))))
        desc_cols.append(desc_col)

    if not cand_exprs:
        return pl.lit(None, dtype=pl.Utf8)

    max_expr = pl.max_horizontal(cand_exprs)
    out = pl.when(max_expr == float("-inf")).then(pl.lit(None, dtype=pl.Utf8))
    for cand_expr, desc_col in zip(cand_exprs, desc_cols):
        out = out.when(cand_expr == max_expr).then(_clean_text_expr(pl.col(desc_col)))
    return out.otherwise(pl.lit(None, dtype=pl.Utf8))

def add_function_column_for_category(df: pl.DataFrame, ref_ids: list[str]) -> pl.DataFrame:
    # Add Function from best-scoring hit restricted to reference IDs; fallback to top hit annotation
    top_desc_expr = _clean_text_expr(pl.col("top_hit_description")) if "top_hit_description" in df.columns else pl.lit(None, dtype=pl.Utf8)
    restricted_best_expr = _best_desc_for_ref_ids(df, ref_ids)
    df = df.with_columns(pl.coalesce([restricted_best_expr, top_desc_expr]).alias("Function"))

    # Place Function immediately after the Genome column
    after_candidates = ["Genome", "genome"]
    after_col = next((c for c in after_candidates if c in df.columns), None)
    if after_col is not None:
        cols = [c for c in df.columns if c != "Function"]
        idx = cols.index(after_col) + 1
        cols.insert(idx, "Function")
        df = df.select(cols)

    return df

def add_function_to_annot_table(
    annot_df: pl.DataFrame,
    metab_df: pl.DataFrame,
    phys_df: pl.DataFrame,
    reg_df: pl.DataFrame,
) -> pl.DataFrame:
    # Populate Function in annot_table from category tables; remaining entries fallback to top hit annotation
    key_cols = ["Protein", "Genome", "Contig", "gene_number"]
    for c in key_cols:
        if c not in annot_df.columns:
            raise ValueError(f"annot_table missing required key column: {c}")

    def _map(df: pl.DataFrame) -> pl.DataFrame:
        for c in key_cols + ["Function"]:
            if c not in df.columns:
                raise ValueError(f"category table missing required column: {c}")
        return df.select(key_cols + ["Function"])

    # Priority order for overlaps: metabolism, physiology, regulation
    function_map = pl.concat([_map(metab_df), _map(phys_df), _map(reg_df)], how="vertical")
    function_map = function_map.unique(subset=key_cols, keep="first")

    annot_df = annot_df.join(function_map, on=key_cols, how="left")

    top_desc_expr = _clean_text_expr(pl.col("top_hit_description")) if "top_hit_description" in annot_df.columns else pl.lit(None, dtype=pl.Utf8)
    annot_df = annot_df.with_columns(pl.coalesce([_clean_text_expr(pl.col("Function")), top_desc_expr]).alias("Function"))

    # Place Function immediately after the Genome column
    after_candidates = ["Genome", "genome"]
    after_col = next((c for c in after_candidates if c in annot_df.columns), None)
    if after_col is not None:
        cols = [c for c in annot_df.columns if c != "Function"]
        idx = cols.index(after_col) + 1
        cols.insert(idx, "Function")
        annot_df = annot_df.select(cols)

    return annot_df

def filter_flagged_annots(
    table: pl.DataFrame,
    flagged_ids: pl.DataFrame,
    valid_hmm_ids,
    filter_presets,
):
    """
    flagged_ids must have:
    - id (string)
    - columns named filter_<EXCEPTION>

    table is the HMMsearch results table with per-db columns like:
    <SRC>_hmm_id, <SRC>_score, <SRC>_coverage, <SRC>_Description (Descriptions not used here)

    - hard matches remove unless exception is enabled in filter_presets
    - no_filter disables all filtering
    - audit provides removed/kept + remove_reason/keep_reason
    """

    # normalize presets (exceptions are case-insensitive)
    if isinstance(filter_presets, str):
        presets = {p.strip().lower() for p in filter_presets.split(",") if p.strip()}
    else:
        presets = {str(p).strip().lower() for p in (filter_presets or [])}

    sources = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "CAMPER", "PHROG"]

    # discover exceptions from the flagged_ids columns
    cols = flagged_ids.columns
    filter_cols = sorted([c for c in cols if c.startswith("filter_")])
    exceptions = sorted({c[len("filter_"):] for c in filter_cols})

    # ensure all filter_ exception columns exist (missing -> False)
    for exc in exceptions:
        fc = f"filter_{exc}"
        if fc not in flagged_ids.columns:
            flagged_ids = flagged_ids.with_columns(pl.lit(False).alias(fc))

    # build ID sets per exception, split into "active" (not allowed) vs "allowed" (in presets)
    flagged_ids = flagged_ids.with_columns(pl.col("id").cast(pl.Utf8))
    filter_ids_active_by_exc = {}
    filter_ids_allowed_by_exc = {}

    for exc in exceptions:
        fc = f"filter_{exc}"

        ids = set(
            flagged_ids
            .filter(pl.col(fc).fill_null(False))
            .select("id")
            .to_series()
            .to_list()
        )

        if exc.strip().lower() in presets:
            filter_ids_allowed_by_exc[exc] = ids
            filter_ids_active_by_exc[exc] = set()
        else:
            filter_ids_active_by_exc[exc] = ids
            filter_ids_allowed_by_exc[exc] = set()

    # helpers: does any <SRC>_hmm_id match any IDs in ids_by_exc?
    def any_flag_match(ids_by_exc):
        exprs = []
        for src in sources:
            id_col = f"{src}_hmm_id"
            if id_col not in table.columns:
                continue
            for ids in ids_by_exc.values():
                if not ids:
                    continue
                exprs.append(pl.col(id_col).cast(pl.Utf8).is_in(ids).fill_null(False))
        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    # helper: first match token for audit, deterministic order
    # token format: "<etype>|<exception>|<src>|<id>"
    def first_flag_token(etype: str, ids_by_exc):
        pieces = []
        for exc in exceptions:
            ids = ids_by_exc.get(exc, set())
            if not ids:
                continue
            for src in sources:
                id_col = f"{src}_hmm_id"
                if id_col not in table.columns:
                    continue
                m = pl.col(id_col).cast(pl.Utf8).is_in(ids).fill_null(False)
                tok = pl.concat_str([
                    pl.lit(f"{etype}|"),
                    pl.lit(f"{exc}|"),
                    pl.lit(f"{src}|"),
                    pl.col(id_col).cast(pl.Utf8),
                ])
                pieces.append(pl.when(m).then(tok).otherwise(pl.lit(None, dtype=pl.Utf8)))
        return pl.coalesce(pieces) if pieces else pl.lit(None, dtype=pl.Utf8)

    # active flags for actual removal
    filter_match_active = any_flag_match(filter_ids_active_by_exc).alias("FILTER_MATCH")

    # all flags for audit (includes allowed exceptions too)
    filter_all_by_exc = {exc: (filter_ids_active_by_exc.get(exc, set()) | filter_ids_allowed_by_exc.get(exc, set())) for exc in exceptions}
    filter_match_all = any_flag_match(filter_all_by_exc).alias("FILTER_MATCH_ALL")
    filter_first_tok_all = first_flag_token("hard", filter_all_by_exc).alias("FILTER_REMOVE_TOKEN_ALL")
    filter_exc_match_all = any_flag_match(filter_ids_allowed_by_exc).alias("FILTER_EXCEPTION_ALL")
    filter_exc_tok_all = first_flag_token("hard", filter_ids_allowed_by_exc).alias("FILTER_EXCEPTION_TOKEN_ALL")

    table_with_flags = table.with_columns([
        filter_match_active,
        filter_match_all,
        filter_first_tok_all,
        filter_exc_match_all,
        filter_exc_tok_all,
    ])

    # filtering
    if "no_filter" in presets:
        table_filtered = table
        removed_expr = pl.lit(False)
    else:
        table_filtered = (
            table_with_flags
            .filter(~pl.col("FILTER_MATCH"))
            .drop("FILTER_MATCH")
        )
        removed_expr = pl.col("FILTER_MATCH")

    kept_expr = ~removed_expr

    # remove_reason: active cause; if kept, show hypothetical cause
    remove_reason_expr = (
        pl.when(removed_expr & pl.col("FILTER_MATCH"))
        .then(pl.col("FILTER_REMOVE_TOKEN_ALL"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )

    kept_side_remove_reason = (
        pl.when(pl.col("FILTER_MATCH_ALL"))
        .then(pl.col("FILTER_REMOVE_TOKEN_ALL"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )
    remove_reason = pl.when(kept_expr).then(kept_side_remove_reason).otherwise(remove_reason_expr)

    # keep_reason with strict precedence, applied only when kept
    # 1) no_filter
    # 2) exception:<token> for hard (must have hard match)
    def token_to_exception(expr: pl.Expr) -> pl.Expr:
        # token is "etype|exception|src|id"
        return expr.str.split("|").list.get(1)

    keep_no_filter = pl.when(kept_expr & pl.lit("no_filter" in presets)).then(pl.lit("no_filter")).otherwise(pl.lit(None, dtype=pl.Utf8))

    keep_exc_filter = pl.when(kept_expr & pl.col("FILTER_EXCEPTION_ALL") & pl.col("FILTER_MATCH_ALL")).then(
        pl.concat_str([pl.lit("exception:"), token_to_exception(pl.col("FILTER_EXCEPTION_TOKEN_ALL"))])
    ).otherwise(pl.lit(None, dtype=pl.Utf8))

    keep_reason_expr = (
        pl.when(keep_no_filter.is_not_null()).then(keep_no_filter)
        .when(keep_exc_filter.is_not_null()).then(keep_exc_filter)
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )
    keep_reason_expr = pl.when(removed_expr).then(pl.lit(None, dtype=pl.Utf8)).otherwise(keep_reason_expr)

    audit = (
        table_with_flags
        .with_columns([
            removed_expr.alias("removed"),
            kept_expr.alias("kept"),
            remove_reason.alias("remove_reason"),
            keep_reason_expr.alias("keep_reason"),
        ])
        .drop([
            "FILTER_MATCH",
            "FILTER_MATCH_ALL",
            "FILTER_REMOVE_TOKEN_ALL",
            "FILTER_EXCEPTION_ALL",
            "FILTER_EXCEPTION_TOKEN_ALL",
        ])
    )

    return table_filtered, audit

def filter_metabolism_annots(table, metabolism_table, flagged_metabolism_table, filter_presets):
    """
    Identify metabolism-related genes based on input metabolism table.
    by checking any of the HMM ID columns for membership in metabolism_table["id"].
    Also, apply id filtering to remove non-metabolic genes.
    """
    metab_ids = metabolism_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(metab_ids) |
        pl.col("FOAM_hmm_id").is_in(metab_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(metab_ids) |
        pl.col("dbCAN_hmm_id").is_in(metab_ids) |
        pl.col("METABOLIC_hmm_id").is_in(metab_ids) |
        pl.col("CAMPER_hmm_id").is_in(metab_ids) |
        pl.col("PHROG_hmm_id").is_in(metab_ids)
    )
    table = table.filter(condition)

    # Apply idfiltering
    table, audit = filter_flagged_annots(table, flagged_metabolism_table, metab_ids, filter_presets)

    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()

    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])

def filter_physiology_annots(table, physiology_table, flagged_phys_table, filter_presets):
    """
    Identify physiology-related genes based on input physiology table.
    by checking any of the HMM ID columns for membership in physiology_table["id"].
    Also, apply id filtering to remove non-physiological genes.
    """
    phys_ids = physiology_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(phys_ids) |
        pl.col("FOAM_hmm_id").is_in(phys_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(phys_ids) |
        pl.col("dbCAN_hmm_id").is_in(phys_ids) |
        pl.col("METABOLIC_hmm_id").is_in(phys_ids) |
        pl.col("CAMPER_hmm_id").is_in(phys_ids) |
        pl.col("PHROG_hmm_id").is_in(phys_ids)
    )
    table = table.filter(condition)

    # Apply id filtering
    table, audit = filter_flagged_annots(table, flagged_phys_table, phys_ids, filter_presets)

    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()

    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])

def filter_regulation_annots(table, regulation_table, flagged_reg_table, filter_presets):
    """
    Identify regulation-related genes based on input regulation table.
    by checking any of the HMM ID columns for membership in regulation_table["id"].
    Also, apply id filtering to remove non-regulatory genes.
    """
    reg_ids = regulation_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(reg_ids) |
        pl.col("FOAM_hmm_id").is_in(reg_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(reg_ids) |
        pl.col("dbCAN_hmm_id").is_in(reg_ids) |
        pl.col("METABOLIC_hmm_id").is_in(reg_ids) |
        pl.col("CAMPER_hmm_id").is_in(reg_ids) |
        pl.col("PHROG_hmm_id").is_in(reg_ids)
    )
    table = table.filter(condition)

    # Apply id filtering
    table, audit = filter_flagged_annots(table, flagged_reg_table, reg_ids, filter_presets)

    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()

    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])

def main():
    input_table  = snakemake.params.context_table
    hmm_ref = snakemake.params.hmm_ref
    metabolism_ref = snakemake.params.metabolism_table
    physiology_ref = snakemake.params.physiology_table
    regulation_ref = snakemake.params.regulation_table
    flagged_metab_table_path = snakemake.params.flagged_amgs
    flagged_phys_table_path = snakemake.params.flagged_apgs
    flagged_reg_table_path = snakemake.params.flagged_aregs
    filter_presets = list(snakemake.params.filter_presets.split(","))

    out_metabolism_table = snakemake.params.metabolism_table_out
    out_metabolism_audit = snakemake.params.metabolism_table_audit
    out_physiology_table = snakemake.params.physiology_table_out
    out_physiology_audit = snakemake.params.physiology_table_audit
    out_regulation_table = snakemake.params.regulation_table_out
    out_regulation_audit = snakemake.params.regulation_table_audit
    all_annot_out_table = snakemake.params.all_annot_out_table
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Starting the curation of annotations for metabolism, physiology, and regulation...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")

    table = pl.read_csv(input_table, separator="\t")
    pl.Config.set_tbl_cols(-1)
    pl.Config.set_tbl_rows(20)
    pl.Config.set_fmt_str_lengths(200)

    hmm_descriptions = pl.read_csv(hmm_ref, schema={"id": pl.Utf8, "db": pl.Utf8, "name": pl.Utf8}, separator="\t")
    hmm_descriptions = hmm_descriptions.select(["id", "db", "name"])

    # Add a normalized ID column for all Pfam entries in hmm_descriptions (strip .number suffix)
    hmm_descriptions = hmm_descriptions.with_columns([
        pl.when(pl.col("db") == "Pfam")
        .then(pl.col("id").str.replace(r"\.\d+$", "", literal=False))
        .otherwise(pl.col("id")).alias("id_norm")
    ])

    metabolism_table = pl.read_csv(metabolism_ref, separator="\t", schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8})
    physiology_table = pl.read_csv(physiology_ref, separator="\t", schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8})
    regulation_table = pl.read_csv(regulation_ref, separator="\t", schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8})

    flagged_metab_table = pl.read_csv(flagged_metab_table_path, separator="\t")
    flagged_phys_table = pl.read_csv(flagged_phys_table_path, separator="\t")
    flagged_reg_table = pl.read_csv(flagged_reg_table_path, separator="\t")

    annot_table = summarize_annot_table(table, hmm_descriptions)

    # Remove .X or .XX suffixes from top_hit_hmm_id for proper matching of Pfam hits
    annot_table = annot_table.with_columns(
        pl.col("top_hit_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("top_hit_hmm_id_clean"),
        pl.col("Pfam_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("Pfam_hmm_id_clean"),
    )

    drop_cols = ["gene_number", "window_avg_KEGG_VL-score_viral", "window_avg_Pfam_VL-score_viral", "window_avg_PHROG_VL-score_viral", "top_hit_hmm_id_clean", "Pfam_hmm_id_clean"]

    metabolism_table_out, metabolism_filter_audit = filter_metabolism_annots(annot_table, metabolism_table, flagged_metab_table, filter_presets)
    physiology_table_out, physiology_filter_audit = filter_physiology_annots(annot_table, physiology_table, flagged_phys_table, filter_presets)
    regulation_table_out, regulation_filter_audit = filter_regulation_annots(annot_table, regulation_table, flagged_reg_table, filter_presets)

    metab_ids = metabolism_table["id"].unique().to_list()
    phys_ids = physiology_table["id"].unique().to_list()
    reg_ids = regulation_table["id"].unique().to_list()

    # Add Function to category tables using category membership and best-scoring reference hits
    metabolism_table_out = add_function_column_for_category(metabolism_table_out, metab_ids)
    physiology_table_out = add_function_column_for_category(physiology_table_out, phys_ids)
    regulation_table_out = add_function_column_for_category(regulation_table_out, reg_ids)

    # Add Function to audit tables using category membership and best-scoring reference hits
    metabolism_filter_audit = add_function_column_for_category(metabolism_filter_audit, metab_ids)
    physiology_filter_audit = add_function_column_for_category(physiology_filter_audit, phys_ids)
    regulation_filter_audit = add_function_column_for_category(regulation_filter_audit, reg_ids)

    # Add Function to annot_table using category tables first; remaining entries fallback to top hit annotation
    annot_table = add_function_to_annot_table(
        annot_table,
        metabolism_table_out,
        physiology_table_out,
        regulation_table_out,
    )

    out_dfs = {
        "annot_table": annot_table,
        "metabolism_table_out": metabolism_table_out,
        "metabolism_filter_audit": metabolism_filter_audit,
        "physiology_table_out": physiology_table_out,
        "physiology_filter_audit": physiology_filter_audit,
        "regulation_table_out": regulation_table_out,
        "regulation_filter_audit": regulation_filter_audit,
    }

    for table_name in out_dfs.keys():
        df = out_dfs[table_name]
        df = df.drop([col for col in df.columns if col in drop_cols])
        replacements = []
        extra_drop_cols = []
        for col in df.columns:
            if "audit" in table_name:
                if not col.endswith("_Description") and not col.endswith("_hmm_id") and not col.endswith("_coverage") and not col.endswith("_score") and not col in ["Protein", "Contig", "Genome", "removed", "kept", "remove_reason", "keep_reason", "Function"]:
                    extra_drop_cols.append(col)
            if col.endswith("_score"):
                replacements.append(
                    pl.when(pl.col(col) == -float("inf"))
                    .then(None)
                    .otherwise(pl.col(col))
                    .alias(col)
                )
        if replacements:
            df = df.with_columns(replacements)
        if extra_drop_cols:
            df = df.drop(extra_drop_cols)
        out_dfs[table_name] = df

    annot_table, \
        metabolism_table_out, metabolism_filter_audit, \
            physiology_table_out, physiology_filter_audit, \
                regulation_table_out, regulation_filter_audit \
                    = out_dfs["annot_table"], \
                        out_dfs["metabolism_table_out"], out_dfs["metabolism_filter_audit"], \
                            out_dfs["physiology_table_out"], out_dfs["physiology_filter_audit"], \
                                out_dfs["regulation_table_out"], out_dfs["regulation_filter_audit"]

    annot_table.write_csv(all_annot_out_table, separator="\t")
    metabolism_table_out.write_csv(out_metabolism_table, separator="\t")
    metabolism_filter_audit.write_csv(out_metabolism_audit, separator="\t")
    physiology_table_out.write_csv(out_physiology_table, separator="\t")
    physiology_filter_audit.write_csv(out_physiology_audit, separator="\t")
    regulation_table_out.write_csv(out_regulation_table, separator="\t")
    regulation_filter_audit.write_csv(out_regulation_audit, separator="\t")

    logger.info("Curation of annotations completed.")
    logger.info(f"Total number of genes analyzed: {annot_table.shape[0]:,}")
    logger.info(f"Number of curated metabolic genes: {metabolism_table_out.shape[0]:,}")
    logger.info(f"Number of curated physiology genes: {physiology_table_out.shape[0]:,}")
    logger.info(f"Number of curated regulatory genes: {regulation_table_out.shape[0]:,}")

if __name__ == "__main__":
    main()
