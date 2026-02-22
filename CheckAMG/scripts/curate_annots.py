#!/usr/bin/env python3

import os
import sys
import logging
import resource
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np

os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
pl.enable_string_cache()


def set_memory_limit(limit_in_gb):
    limit_in_bytes = int(limit_in_gb) * 1024 * 1024 * 1024
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    except (ValueError, OSError, AttributeError) as e:
        logger.warning(f"Unable to set memory limit. Error: {e}")


log_level = logging.DEBUG if snakemake.params.debug else logging.INFO
log_file = snakemake.params.log
logging.basicConfig(
    level=log_level,
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode="a"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger()

print(
    "========================================================================\n"
    "   Step 9/11: Curate the predicted functions based on genomic context   \n"
    "========================================================================"
)
with open(log_file, "a") as log:
    log.write(
        "========================================================================\n"
        "   Step 9/11: Curate the predicted functions based on genomic context   \n"
        "========================================================================\n"
    )


def _as_parquet_path_if_enabled(p, save_to_parquet: bool):
    p = Path(p)
    if save_to_parquet and p.suffix.lower() == ".tsv":
        return p.with_suffix(".parquet")
    return p


def _nan_to_none(v):
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    return float(v)


def _normalize_any_pfam_like_id_expr(id_expr: pl.Expr) -> pl.Expr:
    id_expr = id_expr.cast(pl.Utf8, strict=False)
    id_expr = id_expr.str.strip_chars()
    return (
        pl.when(id_expr.is_null() | (id_expr == ""))
        .then(pl.lit(None, dtype=pl.Utf8))
        .when(id_expr.str.starts_with("PF"))
        .then(id_expr.str.replace(r"\.\d+$", "", literal=False))
        .otherwise(id_expr)
    )


def normalize_reference_table_pfam_ids(df: pl.DataFrame, table_name: str) -> pl.DataFrame:
    if df is None or df.height == 0:
        return df
    if "id" not in df.columns:
        raise ValueError(f"[{table_name}] expected an 'id' column")

    before_pf = int(df.select(pl.col("id").cast(pl.Utf8).str.starts_with("PF").fill_null(False).sum()).item() or 0)
    df = df.with_columns(_normalize_any_pfam_like_id_expr(pl.col("id")).alias("id"))
    after_pf = int(df.select(pl.col("id").cast(pl.Utf8).str.starts_with("PF").fill_null(False).sum()).item() or 0)

    logger.debug(f"[{table_name}] normalized PF* ids (before PF* count={before_pf:,}, after PF* count={after_pf:,}).")
    return df


def dedup_desc_by_id(df, id_col="id", db=None):
    logger.debug(f"Deduplicating descriptions for {db} HMMs by id column: {id_col}")
    df_dedup = (
        df.with_columns(
            pl.col("name").cast(pl.Utf8).fill_null("").alias("name"),
            pl.col("name").cast(pl.Utf8).fill_null("").str.len_bytes().alias("_name_len"),
        )
        .sort([id_col, "_name_len", "name"], descending=[False, True, False])
        .unique(subset=[id_col], keep="first")
        .drop("_name_len")
    )
    return df_dedup


def dedup_desc_by_id_norm(df, norm_col="id_norm", db=None):
    logger.debug(f"Deduplicating descriptions for {db} HMMs by normalized id column: {norm_col}")
    df_dedup = (
        df.with_columns(
            pl.col("name").cast(pl.Utf8).fill_null("").alias("name"),
            pl.col("name").cast(pl.Utf8).fill_null("").str.len_bytes().alias("_name_len"),
        )
        .sort([norm_col, "_name_len", "name"], descending=[False, True, False])
        .unique(subset=[norm_col], keep="first")
        .drop("_name_len")
    )
    return df_dedup


EXC_LABEL_MAP = {
    "no_exception": "other non-auxiliary function",
    "allow_glycosyl": "GH/GT",
    "allow_nucleotide": "nucleotide metabolism",
    "allow_methyl": "methyltransferase",
    "allow_lipid": "lipid metabolism",
}


def validate_audit_table(audit_df: pl.DataFrame, name: str) -> None:
    if audit_df is None or audit_df.height == 0:
        return

    need = {"removed", "kept"}
    missing = need - set(audit_df.columns)
    if missing:
        raise ValueError(f"[{name}] audit is missing required columns: {sorted(missing)}")

    if "remove_reason" not in audit_df.columns:
        removed_n = int(
            audit_df.select(pl.col("removed").cast(pl.Boolean, strict=False).fill_null(False).sum()).item() or 0
        )
        if removed_n > 0:
            raise RuntimeError(f"[{name}] has {removed_n:,} removed rows but no remove_reason column exists.")
        return

    removed_mask = pl.col("removed").cast(pl.Boolean, strict=False).fill_null(False)
    rr = pl.col("remove_reason").cast(pl.Utf8)
    missing_rr = rr.is_null() | (rr.str.strip_chars() == "")

    bad_df = audit_df.filter(removed_mask & missing_rr)
    bad_n = int(bad_df.height)
    if bad_n > 0:
        cols = [c for c in ["Protein", "Genome", "Contig", "gene_number", "remove_reason", "keep_reason"] if c in audit_df.columns]
        sample = bad_df.select(cols).head(50)
        logger.error(f"[{name}] invariant violation: {bad_n:,} rows have removed=True but remove_reason is NULL/blank.")
        logger.error(f"[{name}] sample rows with missing remove_reason (up to 50):\n{sample}")
        raise RuntimeError(f"[{name}] removed rows missing remove_reason (n={bad_n:,}).")


def _pretty_filter_reason_token(tok: str, avg_array_len_limit: int | None, type_name: str) -> str:
    if tok is None:
        return "Unknown"
    tok = str(tok).strip()
    if tok == "":
        return "Unknown"

    if tok.startswith("hard|"):
        parts = tok.split("|")
        exc = parts[1].strip().lower() if len(parts) > 1 else ""
        return EXC_LABEL_MAP.get(exc, exc.replace("_", " ").title() if exc else "other")

    if tok == "context|nonviral_region":
        return "non-viral contig region"

    if tok.startswith("context|avg_array"):
        if avg_array_len_limit is not None:
            return f"{int(avg_array_len_limit)}+ {type_name}s in a row"
        m = re.search(r"avg_array_(?:gt|ge)_(\d+)$", tok)
        if m:
            return f"{int(m.group(1))}+ {type_name}s in a row"
        return f"{type_name} arrays"

    return tok


def log_filter_summary(audit_df: pl.DataFrame, type_name: str, avg_array_len_limit: int | None) -> None:
    if audit_df is None or audit_df.height == 0:
        logger.info(f"Filtered 0 {type_name}s from 0.")
        return

    validate_audit_table(audit_df, f"{type_name}_audit")

    total = int(audit_df.height)
    removed_mask = pl.col("removed").cast(pl.Boolean, strict=False).fill_null(False)

    removed_n = int(audit_df.select(removed_mask.sum()).item() or 0)
    logger.info(f"Filtered {removed_n:,} {type_name}s from {total:,}.")
    if removed_n == 0:
        return

    removed_reasons = (
        audit_df.filter(removed_mask)
        .select(pl.col("remove_reason").cast(pl.Utf8))
        .to_series()
        .to_list()
    )

    counts = Counter()
    for r in removed_reasons:
        toks = [t.strip() for t in str(r).split(",") if t.strip() != ""]
        for tok in toks:
            counts[_pretty_filter_reason_token(tok, avg_array_len_limit, type_name)] += 1

    for reason, n in sorted(counts.items(), key=lambda x: (-x[1], x[0])):
        logger.info(f"Filtered {n:,} {type_name}s for reason: {reason}")


def summarize_annot_table(table, hmm_descriptions, viral_true_col):
    logger.debug(f"Initial table columns: {table.columns}")

    coverage_map = {
        "KEGG": ("KEGG_coverage_hmm", "KEGG_coverage_sequence"),
        "FOAM": ("FOAM_coverage_hmm", "FOAM_coverage_sequence"),
        "Pfam": ("Pfam_coverage_hmm", "Pfam_coverage_sequence"),
        "dbCAN": ("dbCAN_coverage_hmm", "dbCAN_coverage_sequence"),
        "METABOLIC": ("METABOLIC_coverage_hmm", "METABOLIC_coverage_sequence"),
        "CAMPER": ("CAMPER_coverage_hmm", "CAMPER_coverage_sequence"),
        "PHROG": ("PHROG_coverage_hmm", "PHROG_coverage_sequence"),
        # "VOG": ("VOG_coverage_hmm", "VOG_coverage_sequence"),
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

    # These should already exist from the upstream genomic context step, but be defensive
    dist_cols = [
        "KEGG_viral_left_dist", "Pfam_viral_left_dist", "PHROG_viral_left_dist",
        "KEGG_viral_right_dist", "Pfam_viral_right_dist", "PHROG_viral_right_dist",
        "KEGG_MGE_left_dist", "Pfam_MGE_left_dist", "PHROG_MGE_left_dist",
        "KEGG_MGE_right_dist", "Pfam_MGE_right_dist", "PHROG_MGE_right_dist",
    ]
    for c in dist_cols:
        if c not in table.columns:
            table = table.with_columns(pl.lit(None, dtype=pl.Int64).alias(c))

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
        "contig_left_end_gene_dist",
        "contig_right_end_gene_dist",
        "KEGG_hmm_id",
        "FOAM_hmm_id",
        "Pfam_hmm_id",
        "dbCAN_hmm_id",
        "METABOLIC_hmm_id",
        "CAMPER_hmm_id",
        "PHROG_hmm_id",
        # "VOG_hmm_id",
        "KEGG_score",
        "FOAM_score",
        "Pfam_score",
        "dbCAN_score",
        "METABOLIC_score",
        "CAMPER_score",
        "PHROG_score",
        # "VOG_score",
        "KEGG_coverage",
        "FOAM_coverage",
        "Pfam_coverage",
        "dbCAN_coverage",
        "METABOLIC_coverage",
        "CAMPER_coverage",
        "PHROG_coverage",
        # "VOG_coverage",
        "KEGG_V-score",
        "Pfam_V-score",
        "PHROG_V-score",
        # "VOG_V-score",
        "window_avg_KEGG_VL-score_viral",
        "window_avg_Pfam_VL-score_viral",
        "window_avg_PHROG_VL-score_viral",
        "Viral_Flanking_Genes_Left_Dist",
        "Viral_Flanking_Genes_Right_Dist",
        "MGE_Flanking_Genes_Left_Dist",
        "MGE_Flanking_Genes_Right_Dist",
        viral_true_col,
        "nonviral_region_pos",
        "Viral_Origin_Confidence",
    ]

    for col in required_cols:
        if col not in table.columns:
            if col.endswith("_id") or col.endswith("_Description") or col in ["protein", "contig", "genome", "circular_contig", "Viral_Origin_Confidence"]:
                dtype = pl.Utf8
            elif col in ["gene_number", "contig_left_end_gene_dist", "contig_right_end_gene_dist"]:
                dtype = pl.Int64
            elif col in [viral_true_col]:
                dtype = pl.Boolean
            elif col.endswith("_score") or col.endswith("_coverage") or col.endswith("_V-score") or col.endswith("_VL-score"):
                dtype = pl.Float64
            else:
                dtype = pl.Utf8
            table = table.with_columns(pl.lit(None, dtype=dtype).alias(col))

    table = table.with_columns(
        [
            pl.col(viral_true_col).cast(pl.Boolean, strict=False).fill_null(False).alias(viral_true_col),
            pl.col("gene_number").cast(pl.Int64, strict=False).alias("gene_number"),
            pl.col("contig_left_end_gene_dist").cast(pl.Int64, strict=False).alias("contig_left_end_gene_dist"),
            pl.col("contig_right_end_gene_dist").cast(pl.Int64, strict=False).alias("contig_right_end_gene_dist"),
        ]
    )

    table = table.select(required_cols).rename(
        {"protein": "Protein", "contig": "Contig", "genome": "Genome", "circular_contig": "Circular_Contig"}
    )

    hmm_kegg = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "KEGG").select(["id", "name"]), db="KEGG")
    hmm_foam = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "FOAM").select(["id", "name"]), db="FOAM")
    hmm_dbcan = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "dbCAN").select(["id", "name"]), db="dbCAN")
    hmm_metabolic = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "METABOLIC").select(["id", "name"]), db="METABOLIC")
    hmm_camper = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "CAMPER").select(["id", "name"]), db="CAMPER")
    hmm_phrog = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "PHROG").select(["id", "name"]), db="PHROG")
    # hmm_vog = dedup_desc_by_id(hmm_descriptions.filter(pl.col("db") == "VOG").select(["id", "name"]), db="VOG")

    hmm_pfam = (
        hmm_descriptions.filter(pl.col("db") == "Pfam")
        .with_columns(pl.col("id").str.replace(r"\.\d+$", "", literal=False).alias("id_norm"))
        .select(["id_norm", "name"])
    )
    hmm_pfam = dedup_desc_by_id_norm(hmm_pfam, norm_col="id_norm", db="Pfam")

    def _join_and_rename(tbl: pl.DataFrame, hmm_df: pl.DataFrame, left_on: str, right_on: str, db_name: str) -> pl.DataFrame:
        if left_on in tbl.columns:
            tbl = tbl.with_columns(pl.col(left_on).cast(pl.Utf8, strict=False).alias(left_on))
        else:
            tbl = tbl.with_columns(pl.lit(None, dtype=pl.Utf8).alias(left_on))

        if right_on in hmm_df.columns:
            hmm_df = hmm_df.with_columns(pl.col(right_on).cast(pl.Utf8, strict=False).alias(right_on))
        else:
            raise ValueError(f"[{db_name}] hmm_df is missing join key column: {right_on}")

        tbl = tbl.join(hmm_df, left_on=left_on, right_on=right_on, how="left")

        if "name" in tbl.columns:
            tbl = tbl.rename({"name": f"{db_name}_Description"})
        if right_on in tbl.columns:
            tbl = tbl.drop(right_on)

        return tbl

    tbl = table
    tbl = _join_and_rename(tbl, hmm_kegg, "KEGG_hmm_id", "id", "KEGG")
    tbl = _join_and_rename(tbl, hmm_foam, "FOAM_hmm_id", "id", "FOAM")

    tbl = tbl.with_columns(
        pl.col("Pfam_hmm_id")
        .cast(pl.Utf8, strict=False)
        .str.strip_chars()
        .str.replace(r"\.\d+$", "", literal=False)
        .alias("Pfam_hmm_id_clean")
    )
    tbl = _join_and_rename(tbl, hmm_pfam, "Pfam_hmm_id_clean", "id_norm", "Pfam")

    tbl = tbl.with_columns(
        pl.col("dbCAN_hmm_id").cast(pl.Utf8, strict=False).str.replace(r"_(.*)", "", literal=False).alias("dbCAN_hmm_id_clean")
    )
    tbl = _join_and_rename(tbl, hmm_dbcan, "dbCAN_hmm_id_clean", "id", "dbCAN")
    tbl = _join_and_rename(tbl, hmm_metabolic, "METABOLIC_hmm_id", "id", "METABOLIC")
    tbl = _join_and_rename(tbl, hmm_camper, "CAMPER_hmm_id", "id", "CAMPER")
    tbl = _join_and_rename(tbl, hmm_phrog, "PHROG_hmm_id", "id", "PHROG")
    # tbl = _join_and_rename(tbl, hmm_vog, "VOG_hmm_id", "id", "VOG")

    score_cols = [
        "KEGG_score",
        "FOAM_score",
        "Pfam_score",
        "dbCAN_score",
        "METABOLIC_score",
        "CAMPER_score",
        "PHROG_score",
        # "VOG_score",
    ]
    tbl = tbl.with_columns([pl.col(c).cast(pl.Float64).fill_null(float("-inf")).alias(c) for c in score_cols])

    tbl = tbl.with_columns(pl.max_horizontal(score_cols).alias("_max_score"))
    tbl = tbl.with_columns(
        pl.when(pl.col("_max_score") == float("-inf"))
        .then(pl.lit(None, dtype=pl.Int64))
        .otherwise(pl.concat_list(score_cols).list.arg_max())
        .alias("_best_idx")
    ).drop("_max_score")

    tbl = tbl.with_columns(
        [
            pl.when(pl.col("_best_idx") == 0).then(pl.col("KEGG_hmm_id"))
            .when(pl.col("_best_idx") == 1).then(pl.col("FOAM_hmm_id"))
            .when(pl.col("_best_idx") == 2).then(pl.col("Pfam_hmm_id_clean"))
            .when(pl.col("_best_idx") == 3).then(pl.col("dbCAN_hmm_id_clean"))
            .when(pl.col("_best_idx") == 4).then(pl.col("METABOLIC_hmm_id"))
            .when(pl.col("_best_idx") == 5).then(pl.col("CAMPER_hmm_id"))
            .when(pl.col("_best_idx") == 6).then(pl.col("PHROG_hmm_id"))
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("top_hit_hmm_id"),
            pl.when(pl.col("_best_idx") == 0).then(pl.col("KEGG_Description"))
            .when(pl.col("_best_idx") == 1).then(pl.col("FOAM_Description"))
            .when(pl.col("_best_idx") == 2).then(pl.col("Pfam_Description"))
            .when(pl.col("_best_idx") == 3).then(pl.col("dbCAN_Description"))
            .when(pl.col("_best_idx") == 4).then(pl.col("METABOLIC_Description"))
            .when(pl.col("_best_idx") == 5).then(pl.col("CAMPER_Description"))
            .when(pl.col("_best_idx") == 6).then(pl.col("PHROG_Description"))
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("top_hit_description"),
            pl.when(pl.col("_best_idx") == 0).then(pl.lit("KEGG"))
            .when(pl.col("_best_idx") == 1).then(pl.lit("FOAM"))
            .when(pl.col("_best_idx") == 2).then(pl.lit("Pfam"))
            .when(pl.col("_best_idx") == 3).then(pl.lit("dbCAN"))
            .when(pl.col("_best_idx") == 4).then(pl.lit("METABOLIC"))
            .when(pl.col("_best_idx") == 5).then(pl.lit("CAMPER"))
            .when(pl.col("_best_idx") == 6).then(pl.lit("PHROG"))
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("top_hit_db"),
        ]
    ).drop("_best_idx")

    tbl = (
        tbl.unique(subset=["Protein", "Genome", "Contig", "gene_number"], keep="first")
        .sort(["Genome", "Contig", "gene_number"])
    )
    return tbl


def _clean_text_expr(expr: pl.Expr) -> pl.Expr:
    expr = expr.cast(pl.Utf8, strict=False)
    return (
        pl.when(expr.is_null())
        .then(pl.lit(None, dtype=pl.Utf8))
        .otherwise(
            pl.when(expr.str.strip_chars() == "")
            .then(pl.lit(None, dtype=pl.Utf8))
            .otherwise(expr)
        )
    )


def _is_in_ids(id_expr: pl.Expr, ids: list[str]) -> pl.Expr:
    if not ids:
        return pl.lit(False)
    return id_expr.is_in(pl.lit(ids)).fill_null(False)


def _best_desc_for_ref_ids(df: pl.DataFrame, ref_ids: list[str]) -> pl.Expr:
    def _pfam_id_expr():
        if "Pfam_hmm_id_clean" in df.columns:
            return pl.col("Pfam_hmm_id_clean").cast(pl.Utf8, strict=False)
        if "Pfam_hmm_id" in df.columns:
            return pl.col("Pfam_hmm_id").cast(pl.Utf8, strict=False).str.replace(r"\.\d+$", "", literal=False)
        return pl.lit(None, dtype=pl.Utf8)

    def _dbcan_id_expr():
        if "dbCAN_hmm_id_clean" in df.columns:
            return pl.col("dbCAN_hmm_id_clean").cast(pl.Utf8, strict=False)
        if "dbCAN_hmm_id" in df.columns:
            return pl.col("dbCAN_hmm_id").cast(pl.Utf8, strict=False).str.replace(r"_(.*)", "", literal=False)
        return pl.lit(None, dtype=pl.Utf8)

    sources = [
        ("KEGG", "KEGG_score", "KEGG_Description", pl.col("KEGG_hmm_id").cast(pl.Utf8, strict=False) if "KEGG_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("FOAM", "FOAM_score", "FOAM_Description", pl.col("FOAM_hmm_id").cast(pl.Utf8, strict=False) if "FOAM_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("Pfam", "Pfam_score", "Pfam_Description", _pfam_id_expr()),
        ("dbCAN", "dbCAN_score", "dbCAN_Description", _dbcan_id_expr()),
        ("METABOLIC", "METABOLIC_score", "METABOLIC_Description", pl.col("METABOLIC_hmm_id").cast(pl.Utf8, strict=False) if "METABOLIC_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("CAMPER", "CAMPER_score", "CAMPER_Description", pl.col("CAMPER_hmm_id").cast(pl.Utf8, strict=False) if "CAMPER_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
        ("PHROG", "PHROG_score", "PHROG_Description", pl.col("PHROG_hmm_id").cast(pl.Utf8, strict=False) if "PHROG_hmm_id" in df.columns else pl.lit(None, dtype=pl.Utf8)),
    ]

    cand_exprs = []
    desc_cols = []
    for _, score_col, desc_col, id_expr in sources:
        if score_col not in df.columns or desc_col not in df.columns:
            continue
        cond = _is_in_ids(id_expr, ref_ids)
        score = pl.col(score_col).cast(pl.Float64, strict=False).fill_null(float("-inf"))
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
    top_desc_expr = _clean_text_expr(pl.col("top_hit_description")) if "top_hit_description" in df.columns else pl.lit(None, dtype=pl.Utf8)
    restricted_best_expr = _best_desc_for_ref_ids(df, ref_ids)
    df = df.with_columns(pl.coalesce([restricted_best_expr, top_desc_expr]).alias("Function"))

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
    key_cols = ["Protein", "Genome", "Contig", "gene_number"]
    for c in key_cols:
        if c not in annot_df.columns:
            raise ValueError(f"annot_table missing required key column: {c}")

    def _map(df: pl.DataFrame) -> pl.DataFrame:
        for c in key_cols + ["Function"]:
            if c not in df.columns:
                raise ValueError(f"category table missing required column: {c}")
        return df.select(key_cols + ["Function"])

    function_map = pl.concat([_map(metab_df), _map(phys_df), _map(reg_df)], how="vertical")
    function_map = function_map.unique(subset=key_cols, keep="first")

    annot_df = annot_df.join(function_map, on=key_cols, how="left")

    top_desc_expr = _clean_text_expr(pl.col("top_hit_description")) if "top_hit_description" in annot_df.columns else pl.lit(None, dtype=pl.Utf8)
    annot_df = annot_df.with_columns(pl.coalesce([_clean_text_expr(pl.col("Function")), top_desc_expr]).alias("Function"))

    after_candidates = ["Genome", "genome"]
    after_col = next((c for c in after_candidates if c in annot_df.columns), None)
    if after_col is not None:
        cols = [c for c in annot_df.columns if c != "Function"]
        idx = cols.index(after_col) + 1
        cols.insert(idx, "Function")
        annot_df = annot_df.select(cols)
    return annot_df


def _ensure_context_cols(df: pl.DataFrame, viral_true_col: str) -> pl.DataFrame:
    if viral_true_col not in df.columns:
        df = df.with_columns(pl.lit(False).alias(viral_true_col))
    df = df.with_columns(pl.col(viral_true_col).cast(pl.Boolean, strict=False).fill_null(False).alias(viral_true_col))

    viral_false_col = "protein_in_nonviral_region"
    if viral_false_col not in df.columns:
        df = df.with_columns((~pl.col(viral_true_col)).alias(viral_false_col))
    else:
        df = df.with_columns(pl.col(viral_false_col).cast(pl.Boolean, strict=False).fill_null(False).alias(viral_false_col))

    if "nonviral_region_pos" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("nonviral_region_pos"))
    else:
        df = df.with_columns(pl.col("nonviral_region_pos").cast(pl.Utf8, strict=False).alias("nonviral_region_pos"))

    if "protein_in_avg_array" not in df.columns:
        df = df.with_columns(pl.lit(False).alias("protein_in_avg_array"))
    else:
        df = df.with_columns(pl.col("protein_in_avg_array").cast(pl.Boolean, strict=False).fill_null(False).alias("protein_in_avg_array"))

    if "avg_array_length" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Int32).alias("avg_array_length"))
    else:
        df = df.with_columns(pl.col("avg_array_length").cast(pl.Int32, strict=False).alias("avg_array_length"))

    if "protein_in_or_adjacent_nonviral_region" not in df.columns:
        # Adjacent definition:
        # - any nonviral gene is True
        # - any viral gene immediately bordering a nonviral gene is True (prev/next), requiring gene_number consecutiveness
        # - first/last gene are True (conservative)
        need_cols = ["Genome", "Contig", "gene_number"]
        for c in need_cols:
            if c not in df.columns:
                raise ValueError(f"Cannot derive protein_in_or_adjacent_nonviral_region; missing required column: {c}")

        df = df.sort(["Genome", "Contig", "gene_number"])

        df = df.with_columns(
            [
                pl.col("gene_number").min().over(["Genome", "Contig"]).alias("_min_gene_number"),
                pl.col("gene_number").max().over(["Genome", "Contig"]).alias("_max_gene_number"),
            ]
        ).with_columns(
            [
                (pl.col("gene_number") == pl.col("_min_gene_number")).alias("_is_first_gene"),
                (pl.col("gene_number") == pl.col("_max_gene_number")).alias("_is_last_gene"),
            ]
        )

        df = df.with_columns(
            [
                pl.col(viral_false_col).shift(1).over(["Genome", "Contig"]).fill_null(False).alias("_prev_is_nonviral"),
                pl.col("gene_number").shift(1).over(["Genome", "Contig"]).alias("_prev_gene_number"),
                pl.col(viral_false_col).shift(-1).over(["Genome", "Contig"]).fill_null(False).alias("_next_is_nonviral"),
                pl.col("gene_number").shift(-1).over(["Genome", "Contig"]).alias("_next_gene_number"),
            ]
        )

        prev_consecutive = (pl.col("gene_number") - pl.col("_prev_gene_number") == 1).fill_null(False)
        next_consecutive = (pl.col("_next_gene_number") - pl.col("gene_number") == 1).fill_null(False)

        borders_prev_nonviral = (pl.col("_prev_is_nonviral") & prev_consecutive)
        borders_next_nonviral = (pl.col("_next_is_nonviral") & next_consecutive)

        df = df.with_columns(
            (
                pl.col(viral_false_col)
                | (pl.col(viral_true_col) & (borders_prev_nonviral | borders_next_nonviral))
                | pl.col("_is_first_gene")
                | pl.col("_is_last_gene")
            ).alias("protein_in_or_adjacent_nonviral_region")
        ).drop(
            [
                "_min_gene_number",
                "_max_gene_number",
                "_is_first_gene",
                "_is_last_gene",
                "_prev_is_nonviral",
                "_prev_gene_number",
                "_next_is_nonviral",
                "_next_gene_number",
            ]
        )
    else:
        df = df.with_columns(
            pl.col("protein_in_or_adjacent_nonviral_region").cast(pl.Boolean, strict=False).fill_null(False).alias("protein_in_or_adjacent_nonviral_region")
        )

    return df


def _attach_full_context_flags(
    df: pl.DataFrame,
    context_df: pl.DataFrame,
    viral_true_col: str,
) -> pl.DataFrame:
    # Attach context-derived flags computed on the full annotation table (all genes on contig),
    # so boundary/adjacency checks remain valid after subsetting to AMGs/APGs/AReGs.
    key_cols = ["Protein", "Genome", "Contig", "gene_number"]
    for c in key_cols:
        if c not in df.columns:
            raise ValueError(f"Cannot attach context flags; df missing required key column: {c}")
        if c not in context_df.columns:
            raise ValueError(f"Cannot attach context flags; context_df missing required key column: {c}")

    need_cols = [
        viral_true_col,
        "protein_in_nonviral_region",
        "protein_in_or_adjacent_nonviral_region",
        "nonviral_region_pos",
    ]
    missing = [c for c in need_cols if c not in context_df.columns]
    if missing:
        raise ValueError(f"Cannot attach context flags; context_df missing required columns: {missing}")

    ctx = context_df.select(key_cols + need_cols).unique(subset=key_cols, keep="first")

    # Avoid duplicate columns on join
    drop_existing = [c for c in need_cols if c in df.columns]
    if drop_existing:
        df = df.drop(drop_existing, strict=False)

    return df.join(ctx, on=key_cols, how="left")


def _context_remove_expr(filter_nonviral_regions: bool, filter_avg_arrays: bool, avg_array_len_limit: int | None) -> pl.Expr:
    nonviral = pl.lit(False)
    if filter_nonviral_regions:
        nonviral = pl.col("protein_in_or_adjacent_nonviral_region").cast(pl.Boolean, strict=False).fill_null(False)

    avg = pl.lit(False)
    if filter_avg_arrays:
        in_arr = pl.col("protein_in_avg_array").cast(pl.Boolean, strict=False).fill_null(False)
        if avg_array_len_limit is None:
            avg = in_arr
        else:
            avg = in_arr & (pl.col("avg_array_length").cast(pl.Int32, strict=False) >= int(avg_array_len_limit))

    return nonviral | avg


def _context_reason_expr(filter_nonviral_regions: bool, filter_avg_arrays: bool, avg_array_len_limit: int | None) -> pl.Expr:
    parts = []

    if filter_nonviral_regions:
        parts.append(
            pl.when(pl.col("protein_in_or_adjacent_nonviral_region").cast(pl.Boolean, strict=False).fill_null(False))
            .then(pl.lit("context|nonviral_region"))
            .otherwise(pl.lit(None, dtype=pl.Utf8))
        )

    if filter_avg_arrays:
        in_arr = pl.col("protein_in_avg_array").cast(pl.Boolean, strict=False).fill_null(False)
        if avg_array_len_limit is None:
            parts.append(
                pl.when(in_arr)
                .then(pl.lit("context|avg_array"))
                .otherwise(pl.lit(None, dtype=pl.Utf8))
            )
        else:
            parts.append(
                pl.when(in_arr & (pl.col("avg_array_length").cast(pl.Int32, strict=False) >= int(avg_array_len_limit)))
                .then(pl.lit(f"context|avg_array_ge_{int(avg_array_len_limit)}"))
                .otherwise(pl.lit(None, dtype=pl.Utf8))
            )

    if not parts:
        return pl.lit(None, dtype=pl.Utf8)

    if len(parts) == 1:
        return parts[0]

    a, b = parts[0], parts[1]
    return (
        pl.when(a.is_not_null() & b.is_not_null()).then(pl.concat_str([a, pl.lit(","), b]))
        .when(a.is_not_null()).then(a)
        .when(b.is_not_null()).then(b)
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )


def apply_context_filters_to_predictions(
    df: pl.DataFrame,
    filter_nonviral_regions: bool,
    filter_avg_arrays: bool,
    avg_array_len_limit: int | None,
    viral_true_col: str,
    context_df: pl.DataFrame,
) -> pl.DataFrame:
    df = _attach_full_context_flags(df, context_df=context_df, viral_true_col=viral_true_col)
    df = _ensure_context_cols(df, viral_true_col=viral_true_col)
    remove_expr = _context_remove_expr(filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit)
    return df.filter(~remove_expr)


def apply_context_filters_to_audit(
    audit_df: pl.DataFrame,
    filter_nonviral_regions: bool,
    filter_avg_arrays: bool,
    avg_array_len_limit: int | None,
    viral_true_col: str,
    context_df: pl.DataFrame,
) -> pl.DataFrame:
    audit_df = _attach_full_context_flags(audit_df, context_df=context_df, viral_true_col=viral_true_col)
    audit_df = _ensure_context_cols(audit_df, viral_true_col=viral_true_col)

    if "removed" not in audit_df.columns or "kept" not in audit_df.columns:
        raise ValueError("Audit table is missing required columns: removed, kept")

    if "remove_reason" not in audit_df.columns:
        audit_df = audit_df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("remove_reason"))
    if "keep_reason" not in audit_df.columns:
        audit_df = audit_df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("keep_reason"))

    prev_removed = pl.col("removed").cast(pl.Boolean, strict=False).fill_null(False)
    prev_kept = pl.col("kept").cast(pl.Boolean, strict=False).fill_null(False)
    prev_remove_reason = pl.col("remove_reason").cast(pl.Utf8, strict=False)
    prev_keep_reason = pl.col("keep_reason").cast(pl.Utf8, strict=False)

    missing_prev_reason = prev_remove_reason.is_null() | (prev_remove_reason.str.strip_chars() == "")
    ghost_removed_expr = (prev_removed & missing_prev_reason).alias("_GHOST_REMOVED")

    ctx_remove_expr = _context_remove_expr(filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit).alias("_CTX_REMOVE")
    ctx_reason_expr = _context_reason_expr(filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit).alias("_CTX_REASON")

    audit_df = audit_df.with_columns([ghost_removed_expr, ctx_remove_expr, ctx_reason_expr])

    ghost_n = int(audit_df.select(pl.col("_GHOST_REMOVED").sum()).item() or 0)
    if ghost_n > 0:
        logger.warning(f"Reverting {ghost_n:,} ghost-removed rows (removed=True but remove_reason NULL/blank) back to kept=True.")

    sanitized_kept_expr = (prev_kept | pl.col("_GHOST_REMOVED")).alias("_SANITIZED_KEPT")
    sanitized_removed_expr = (prev_removed & ~pl.col("_GHOST_REMOVED")).alias("_SANITIZED_REMOVED")
    audit_df = audit_df.with_columns([sanitized_kept_expr, sanitized_removed_expr])

    newly_removed_expr = (pl.col("_CTX_REMOVE") & pl.col("_SANITIZED_KEPT")).alias("_NEWLY_REMOVED")
    audit_df = audit_df.with_columns([newly_removed_expr])

    audit_df = audit_df.with_columns(
        [
            (pl.col("_SANITIZED_REMOVED") | pl.col("_NEWLY_REMOVED")).alias("removed"),
            (pl.col("_SANITIZED_KEPT") & ~pl.col("_NEWLY_REMOVED")).alias("kept"),
        ]
    )

    audit_df = audit_df.with_columns(
        [
            pl.when(pl.col("_NEWLY_REMOVED"))
            .then(pl.col("_CTX_REASON"))
            .when(pl.col("_SANITIZED_REMOVED"))
            .then(prev_remove_reason)
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("remove_reason"),
            pl.when(pl.col("kept"))
            .then(prev_keep_reason)
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("keep_reason"),
        ]
    )

    audit_df = audit_df.drop(
        [
            "_GHOST_REMOVED",
            "_CTX_REMOVE",
            "_CTX_REASON",
            "_SANITIZED_KEPT",
            "_SANITIZED_REMOVED",
            "_NEWLY_REMOVED",
        ]
    )

    validate_audit_table(audit_df, "context_filter_audit")
    return audit_df


def filter_flagged_annots(
    table: pl.DataFrame,
    flagged_ids: pl.DataFrame,
    valid_hmm_ids,
    filter_presets,
):
    if isinstance(filter_presets, str):
        presets = {p.strip().lower() for p in filter_presets.split(",") if p.strip()}
    else:
        presets = {str(p).strip().lower() for p in (filter_presets or [])}

    # sources = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "CAMPER", "PHROG", "VOG"]
    sources = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "CAMPER", "PHROG"]

    cols = flagged_ids.columns
    filter_cols = sorted([c for c in cols if c.startswith("filter_")])
    exceptions = sorted({c[len("filter_"):] for c in filter_cols})

    for exc in exceptions:
        fc = f"filter_{exc}"
        if fc not in flagged_ids.columns:
            flagged_ids = flagged_ids.with_columns(pl.lit(False).alias(fc))

    flagged_ids = flagged_ids.with_columns(_normalize_any_pfam_like_id_expr(pl.col("id")).alias("id"))

    filter_ids_active_by_exc = {}
    filter_ids_allowed_by_exc = {}

    for exc in exceptions:
        fc = f"filter_{exc}"
        ids = set(
            flagged_ids.filter(pl.col(fc).fill_null(False))
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

    def _id_expr_for_src(src: str) -> pl.Expr:
        if src == "Pfam":
            if "Pfam_hmm_id_clean" in table.columns:
                return pl.col("Pfam_hmm_id_clean").cast(pl.Utf8, strict=False)
            if "Pfam_hmm_id" in table.columns:
                return pl.col("Pfam_hmm_id").cast(pl.Utf8, strict=False).str.replace(r"\.\d+$", "", literal=False)
            return pl.lit(None, dtype=pl.Utf8)

        if src == "dbCAN":
            if "dbCAN_hmm_id_clean" in table.columns:
                return pl.col("dbCAN_hmm_id_clean").cast(pl.Utf8, strict=False)
            if "dbCAN_hmm_id" in table.columns:
                return pl.col("dbCAN_hmm_id").cast(pl.Utf8, strict=False).str.replace(r"_(.*)", "", literal=False)
            return pl.lit(None, dtype=pl.Utf8)

        id_col = f"{src}_hmm_id"
        if id_col in table.columns:
            return pl.col(id_col).cast(pl.Utf8, strict=False)
        return pl.lit(None, dtype=pl.Utf8)

    src_id_exprs = {src: _id_expr_for_src(src) for src in sources}

    def any_flag_match(ids_by_exc):
        exprs = []
        for src in sources:
            id_expr = src_id_exprs[src]
            for ids in ids_by_exc.values():
                if not ids:
                    continue
                exprs.append(_is_in_ids(id_expr, list(ids)))
        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    def first_flag_token(etype: str, ids_by_exc):
        pieces = []
        for exc in exceptions:
            ids = ids_by_exc.get(exc, set())
            if not ids:
                continue
            ids_list = list(ids)
            for src in sources:
                id_expr = src_id_exprs[src]
                m = _is_in_ids(id_expr, ids_list)
                tok = pl.concat_str(
                    [
                        pl.lit(f"{etype}|"),
                        pl.lit(f"{exc}|"),
                        pl.lit(f"{src}|"),
                        id_expr,
                    ]
                )
                pieces.append(pl.when(m).then(tok).otherwise(pl.lit(None, dtype=pl.Utf8)))
        return pl.coalesce(pieces) if pieces else pl.lit(None, dtype=pl.Utf8)

    filter_match_active = any_flag_match(filter_ids_active_by_exc).alias("FILTER_MATCH")

    filter_all_by_exc = {
        exc: (filter_ids_active_by_exc.get(exc, set()) | filter_ids_allowed_by_exc.get(exc, set()))
        for exc in exceptions
    }
    filter_match_all = any_flag_match(filter_all_by_exc).alias("FILTER_MATCH_ALL")
    filter_first_tok_all = first_flag_token("hard", filter_all_by_exc).alias("FILTER_REMOVE_TOKEN_ALL")
    filter_exc_match_all = any_flag_match(filter_ids_allowed_by_exc).alias("FILTER_EXCEPTION_ALL")
    filter_exc_tok_all = first_flag_token("hard", filter_ids_allowed_by_exc).alias("FILTER_EXCEPTION_TOKEN_ALL")

    table_with_flags = table.with_columns(
        [
            filter_match_active,
            filter_match_all,
            filter_first_tok_all,
            filter_exc_match_all,
            filter_exc_tok_all,
        ]
    )

    if "no_filter" in presets:
        table_filtered = table
        removed_expr = pl.lit(False)
    else:
        table_filtered = table_with_flags.filter(~pl.col("FILTER_MATCH")).drop("FILTER_MATCH")
        removed_expr = pl.col("FILTER_MATCH")

    kept_expr = ~removed_expr

    remove_reason = (
        pl.when(removed_expr)
        .then(pl.col("FILTER_REMOVE_TOKEN_ALL"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
        .alias("remove_reason")
    )

    def token_to_exception(expr: pl.Expr) -> pl.Expr:
        return expr.str.split("|").list.get(1)

    keep_no_filter = (
        pl.when(kept_expr & pl.lit("no_filter" in presets))
        .then(pl.lit("no_filter"))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )

    keep_exc_filter = (
        pl.when(kept_expr & pl.col("FILTER_EXCEPTION_ALL") & pl.col("FILTER_MATCH_ALL"))
        .then(pl.concat_str([pl.lit("exception:"), token_to_exception(pl.col("FILTER_EXCEPTION_TOKEN_ALL"))]))
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )

    keep_reason_expr = (
        pl.when(keep_no_filter.is_not_null()).then(keep_no_filter)
        .when(keep_exc_filter.is_not_null()).then(keep_exc_filter)
        .otherwise(pl.lit(None, dtype=pl.Utf8))
    )
    keep_reason_expr = pl.when(removed_expr).then(pl.lit(None, dtype=pl.Utf8)).otherwise(keep_reason_expr)

    audit = (
        table_with_flags.with_columns(
            [
                removed_expr.alias("removed"),
                kept_expr.alias("kept"),
                remove_reason,
                keep_reason_expr.alias("keep_reason"),
            ]
        )
        .drop(
            [
                "FILTER_MATCH",
                "FILTER_MATCH_ALL",
                "FILTER_REMOVE_TOKEN_ALL",
                "FILTER_EXCEPTION_ALL",
                "FILTER_EXCEPTION_TOKEN_ALL",
            ]
        )
    )

    validate_audit_table(audit, "filter_flagged_annots_audit")
    return table_filtered, audit


def _safe_bool_col(df: pl.DataFrame, name: str, default=False) -> pl.Expr:
    if name not in df.columns:
        return pl.lit(bool(default))
    return pl.col(name).cast(pl.Boolean, strict=False).fill_null(bool(default))


def _safe_utf8_col(df: pl.DataFrame, name: str) -> pl.Expr:
    if name not in df.columns:
        return pl.lit(None, dtype=pl.Utf8)
    return pl.col(name).cast(pl.Utf8, strict=False)


def _category_mask_expr(df: pl.DataFrame, ids: list[str]) -> pl.Expr:
    return (
        _is_in_ids(_safe_utf8_col(df, "KEGG_hmm_id"), ids)
        | _is_in_ids(_safe_utf8_col(df, "FOAM_hmm_id"), ids)
        | _is_in_ids(_safe_utf8_col(df, "Pfam_hmm_id_clean"), ids)
        | _is_in_ids(_safe_utf8_col(df, "dbCAN_hmm_id_clean"), ids)
        | _is_in_ids(_safe_utf8_col(df, "METABOLIC_hmm_id"), ids)
        | _is_in_ids(_safe_utf8_col(df, "CAMPER_hmm_id"), ids)
        | _is_in_ids(_safe_utf8_col(df, "PHROG_hmm_id"), ids)
    )


def filter_metabolism_annots(table: pl.DataFrame, metab_ids: list[str], flagged_metabolism_table: pl.DataFrame, filter_presets):
    table = table.filter(_category_mask_expr(table, metab_ids))
    table, audit = filter_flagged_annots(table, flagged_metabolism_table, metab_ids, filter_presets)
    table = table.unique()
    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])


def filter_physiology_annots(table: pl.DataFrame, phys_ids: list[str], flagged_phys_table: pl.DataFrame, filter_presets):
    table = table.filter(_category_mask_expr(table, phys_ids))
    table, audit = filter_flagged_annots(table, flagged_phys_table, phys_ids, filter_presets)
    table = table.unique()
    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])


def filter_regulation_annots(table: pl.DataFrame, reg_ids: list[str], flagged_reg_table: pl.DataFrame, filter_presets):
    table = table.filter(_category_mask_expr(table, reg_ids))
    table, audit = filter_flagged_annots(table, flagged_reg_table, reg_ids, filter_presets)
    table = table.unique()
    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])


def add_avg_array_features(annot_df: pl.DataFrame, member_keys: pl.DataFrame) -> pl.DataFrame:
    key_cols = ["Protein", "Genome", "Contig", "gene_number"]
    for c in key_cols:
        if c not in annot_df.columns:
            raise ValueError(f"annot_df missing required column for avg array features: {c}")
    for c in key_cols:
        if c not in member_keys.columns:
            raise ValueError(f"member_keys missing required column for avg array features: {c}")

    member_keys = member_keys.select(key_cols).unique(subset=key_cols, keep="first")

    df = (
        annot_df.join(
            member_keys.with_columns(pl.lit(True).alias("_is_member")),
            on=key_cols,
            how="left",
        )
        .with_columns(pl.col("_is_member").fill_null(False).alias("_is_member"))
        .sort(["Genome", "Contig", "gene_number"])
    )

    keys = ["Genome", "Contig"]

    df = df.with_columns(pl.col("_is_member").rle_id().over(keys).cast(pl.Int32).alias("_run_id"))
    df = df.with_columns(pl.len().over(keys + ["_run_id"]).cast(pl.Int32).alias("_run_len"))

    df = df.with_columns(
        [
            (pl.col("_is_member") & (pl.col("_run_len") >= 2)).alias("protein_in_avg_array"),
            pl.when(pl.col("_is_member") & (pl.col("_run_len") >= 2))
            .then(pl.col("_run_len"))
            .otherwise(pl.lit(None, dtype=pl.Int32))
            .alias("avg_array_length"),
        ]
    )

    return df.drop(["_is_member", "_run_id", "_run_len"])


def _parse_int_or_none(v):
    if v is None:
        return None
    if isinstance(v, bool):
        return None
    if isinstance(v, int):
        return v
    s = str(v).strip()
    if s == "" or s.lower() in {"none", "null", "na"}:
        return None
    return int(float(s))


def _parse_presets(v):
    if v is None:
        return []
    if isinstance(v, (list, tuple)):
        return [str(x).strip() for x in v if str(x).strip() != ""]
    return [p.strip() for p in str(v).split(",") if p.strip() != ""]


def main():
    save_to_parquet = bool(snakemake.params.save_to_parquet)
    filter_presets = _parse_presets(snakemake.params.filter_presets)

    filter_nonviral_regions = bool(snakemake.params.filter_nonviral_regions)
    filter_avg_arrays = bool(snakemake.params.filter_avg_arrays)
    avg_array_len_limit = _parse_int_or_none(snakemake.params.avg_array_len_limit)

    viral_true_col = str(snakemake.params.viral_true_col).strip()

    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Starting the curation of annotations for metabolism, physiology, and regulation...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    logger.debug(
        f"Context filters enabled: filter_nonviral_regions={filter_nonviral_regions}, filter_avg_arrays={filter_avg_arrays}, avg_array_len_limit={avg_array_len_limit}"
    )

    # Adjust input and output paths based on whether we're saving to Parquet or TSV
    input_table = _as_parquet_path_if_enabled(snakemake.params.context_table, save_to_parquet)
    out_metabolism_table = _as_parquet_path_if_enabled(snakemake.params.metabolism_table_out, save_to_parquet)
    out_metabolism_audit = _as_parquet_path_if_enabled(snakemake.params.metabolism_table_audit, save_to_parquet)
    out_physiology_table = _as_parquet_path_if_enabled(snakemake.params.physiology_table_out, save_to_parquet)
    out_physiology_audit = _as_parquet_path_if_enabled(snakemake.params.physiology_table_audit, save_to_parquet)
    out_regulation_table = _as_parquet_path_if_enabled(snakemake.params.regulation_table_out, save_to_parquet)
    out_regulation_audit = _as_parquet_path_if_enabled(snakemake.params.regulation_table_audit, save_to_parquet)
    all_annot_out_table = _as_parquet_path_if_enabled(snakemake.params.all_annot_out_table, save_to_parquet)

    if save_to_parquet:
        table = pl.read_parquet(input_table)
    else:
        table = pl.read_csv(input_table, separator="\t")

    pl.Config.set_tbl_cols(-1)
    pl.Config.set_tbl_rows(20)
    pl.Config.set_fmt_str_lengths(200)

    # HMM descriptions and reference tables are always TSV
    hmm_ref = snakemake.params.hmm_ref
    metabolism_ref = snakemake.params.metabolism_table
    physiology_ref = snakemake.params.physiology_table
    regulation_ref = snakemake.params.regulation_table
    flagged_metab_table_path = snakemake.params.flagged_amgs
    flagged_phys_table_path = snakemake.params.flagged_apgs
    flagged_reg_table_path = snakemake.params.flagged_aregs

    hmm_descriptions = pl.read_csv(hmm_ref, schema={"id": pl.Utf8, "db": pl.Utf8, "name": pl.Utf8}, separator="\t")
    hmm_descriptions = hmm_descriptions.select(["id", "db", "name"])

    metabolism_table = pl.read_csv(
        metabolism_ref,
        separator="\t",
        schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8},
    )
    physiology_table = pl.read_csv(
        physiology_ref,
        separator="\t",
        schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8},
    )
    regulation_table = pl.read_csv(
        regulation_ref,
        separator="\t",
        schema={"id": pl.Utf8, "V-score": pl.Float32, "VL-score": pl.Float32, "db": pl.Utf8, "name": pl.Utf8},
    )

    metabolism_table = normalize_reference_table_pfam_ids(metabolism_table, "metabolism_table")
    physiology_table = normalize_reference_table_pfam_ids(physiology_table, "physiology_table")
    regulation_table = normalize_reference_table_pfam_ids(regulation_table, "regulation_table")

    flagged_metab_table = pl.read_csv(flagged_metab_table_path, separator="\t")
    flagged_phys_table = pl.read_csv(flagged_phys_table_path, separator="\t")
    flagged_reg_table = pl.read_csv(flagged_reg_table_path, separator="\t")

    flagged_metab_table = normalize_reference_table_pfam_ids(flagged_metab_table, "flagged_metab_table")
    flagged_phys_table = normalize_reference_table_pfam_ids(flagged_phys_table, "flagged_phys_table")
    flagged_reg_table = normalize_reference_table_pfam_ids(flagged_reg_table, "flagged_reg_table")

    annot_table = summarize_annot_table(table, hmm_descriptions, viral_true_col)

    metab_ids = metabolism_table["id"].unique().to_list()
    phys_ids = physiology_table["id"].unique().to_list()
    reg_ids = regulation_table["id"].unique().to_list()

    key_cols = ["Protein", "Genome", "Contig", "gene_number"]

    logger.debug(f"Building membership keys for avg-array detection (pre-filtering), avg_array_len_limit={avg_array_len_limit}...")
    annot_table = annot_table.with_columns(
        [
            _category_mask_expr(annot_table, metab_ids).alias("_is_metab"),
            _category_mask_expr(annot_table, phys_ids).alias("_is_phys"),
            _category_mask_expr(annot_table, reg_ids).alias("_is_reg"),
        ]
    )
    member_keys = (
        annot_table
        .filter(pl.col("_is_metab") | pl.col("_is_phys") | pl.col("_is_reg"))
        .select(key_cols)
        .unique(subset=key_cols, keep="first")
    )

    if avg_array_len_limit is not None:
        logger.debug(f"Detecting arrays of AVGs in a row of {avg_array_len_limit:,} or more (pre-keyword filtering)...")
    else:
        logger.debug("Detecting arrays of AVGs in a row (pre-keyword filtering)...")
    annot_table = add_avg_array_features(annot_table.drop(["_is_metab", "_is_phys", "_is_reg"]), member_keys)

    # Compute nonviral boundary/adjacency flags on the full annotation table (all genes on contig), before subsetting.
    annot_table = _ensure_context_cols(annot_table, viral_true_col=viral_true_col)

    # Columns that are only useful internally for this step and should not be carried to the final outputs
    internal_drop_cols = [
        "gene_number",
        "window_avg_KEGG_VL-score_viral",
        "window_avg_Pfam_VL-score_viral",
        "window_avg_PHROG_VL-score_viral",
        "Pfam_hmm_id_clean",
        "dbCAN_hmm_id_clean",
    ]

    logger.info("Filtering metabolic annotations...")
    metabolism_table_out, metabolism_filter_audit = filter_metabolism_annots(
        annot_table, metab_ids, flagged_metab_table, filter_presets
    )
    logger.info("Filtering physiological annotations...")
    physiology_table_out, physiology_filter_audit = filter_physiology_annots(
        annot_table, phys_ids, flagged_phys_table, filter_presets
    )
    logger.info("Filtering regulation annotations...")
    regulation_table_out, regulation_filter_audit = filter_regulation_annots(
        annot_table, reg_ids, flagged_reg_table, filter_presets
    )

    metabolism_table_out = add_function_column_for_category(metabolism_table_out, metab_ids)
    physiology_table_out = add_function_column_for_category(physiology_table_out, phys_ids)
    regulation_table_out = add_function_column_for_category(regulation_table_out, reg_ids)

    metabolism_filter_audit = add_function_column_for_category(metabolism_filter_audit, metab_ids)
    physiology_filter_audit = add_function_column_for_category(physiology_filter_audit, phys_ids)
    regulation_filter_audit = add_function_column_for_category(regulation_filter_audit, reg_ids)

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

    if filter_nonviral_regions and filter_avg_arrays:
        if avg_array_len_limit is None:
            logger.info("Filtering AVGs encoded on non-viral contig regions or present in AVG arrays...")
        else:
            logger.info(f"Filtering AVGs encoded on non-viral contig regions or present in arrays of {avg_array_len_limit:,} or more AVGs in a row...")
    elif filter_avg_arrays:
        if avg_array_len_limit is None:
            logger.info("Filtering AVGs present in AVG arrays...")
        else:
            logger.info(f"Filtering AVGs present in arrays of {avg_array_len_limit:,} or more AVGs in a row...")
    elif filter_nonviral_regions:
        logger.info("Filtering AVGs encoded on non-viral contig regions...")
    else:
        logger.warning("No context filters enabled (AVG arrays, non-viral contig regions), skipping context-based filtering.")

    out_dfs["metabolism_table_out"] = apply_context_filters_to_predictions(
        out_dfs["metabolism_table_out"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )
    out_dfs["physiology_table_out"] = apply_context_filters_to_predictions(
        out_dfs["physiology_table_out"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )
    out_dfs["regulation_table_out"] = apply_context_filters_to_predictions(
        out_dfs["regulation_table_out"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )

    annot_df = out_dfs["annot_table"]
    annot_df = _ensure_context_cols(annot_df, viral_true_col=viral_true_col)

    annot_df = (
        annot_df.join(
            member_keys.with_columns(pl.lit(True).alias("_is_member")),
            on=key_cols,
            how="left",
        )
        .with_columns(pl.col("_is_member").fill_null(False).alias("_is_member"))
    )

    annot_remove = _context_remove_expr(filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit)
    annot_reason = _context_reason_expr(filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit)

    annot_df = annot_df.with_columns(
        [
            (pl.col("_is_member") & annot_remove).alias("filtered_from_final_predictions"),
            pl.when(pl.col("_is_member") & annot_remove)
            .then(annot_reason)
            .otherwise(pl.lit(None, dtype=pl.Utf8))
            .alias("filtered_from_final_reason"),
        ]
    ).drop("_is_member")

    out_dfs["annot_table"] = annot_df

    out_dfs["metabolism_filter_audit"] = apply_context_filters_to_audit(
        out_dfs["metabolism_filter_audit"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )
    out_dfs["physiology_filter_audit"] = apply_context_filters_to_audit(
        out_dfs["physiology_filter_audit"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )
    out_dfs["regulation_filter_audit"] = apply_context_filters_to_audit(
        out_dfs["regulation_filter_audit"], filter_nonviral_regions, filter_avg_arrays, avg_array_len_limit, viral_true_col, context_df=annot_table
    )

    validate_audit_table(out_dfs["metabolism_filter_audit"], "metabolism_filter_audit")
    validate_audit_table(out_dfs["physiology_filter_audit"], "physiology_filter_audit")
    validate_audit_table(out_dfs["regulation_filter_audit"], "regulation_filter_audit")

    for table_name, df in list(out_dfs.items()):
        df = df.drop([c for c in internal_drop_cols if c in df.columns], strict=False)

        score_cols = [c for c in df.columns if c.endswith("_score")]
        if score_cols:
            df = df.with_columns(
                [pl.when(pl.col(c) == -float("inf")).then(pl.lit(None, dtype=pl.Float64)).otherwise(pl.col(c)).alias(c) for c in score_cols]
            )

        if "audit" in table_name:
            keep_cols = {
                "Protein",
                "Contig",
                "Genome",
                "gene_number",
                "contig_left_end_gene_dist",
                "contig_right_end_gene_dist",
                "removed",
                "kept",
                "remove_reason",
                "keep_reason",
                "Function",
                viral_true_col,
                "nonviral_region_pos",
                "protein_in_avg_array",
                "avg_array_length",
            }
            keep_cols |= {c for c in df.columns if c.endswith("_Description")}
            keep_cols |= {c for c in df.columns if c.endswith("_hmm_id")}
            keep_cols |= {c for c in df.columns if c.endswith("_coverage")}
            keep_cols |= {c for c in df.columns if c.endswith("_score")}
            keep_cols |= {c for c in df.columns if c.endswith("_V-score")}
            keep_cols |= {c for c in df.columns if c.endswith("_VL-score")}
            keep_cols |= {c for c in df.columns if c in {"Viral_Flanking_Genes_Left_Dist", "Viral_Flanking_Genes_Right_Dist", "MGE_Flanking_Genes_Left_Dist", "MGE_Flanking_Genes_Right_Dist"}}

            df = df.select([c for c in df.columns if c in keep_cols])

        out_dfs[table_name] = df

    drop_cols_final = [
        "Circular_Contig",
        "contig_left_end_gene_dist",
        "contig_right_end_gene_dist",
        "protein_in_avg_array",
        "avg_array_length",
        "protein_in_or_adjacent_nonviral_region",
        "filtered_from_final_predictions",
        "filtered_from_final_reason",
    ]

    annot_table = out_dfs["annot_table"].drop(drop_cols_final, strict=False)
    metabolism_table_out = out_dfs["metabolism_table_out"].drop(drop_cols_final, strict=False)
    metabolism_filter_audit = out_dfs["metabolism_filter_audit"]
    physiology_table_out = out_dfs["physiology_table_out"].drop(drop_cols_final, strict=False)
    physiology_filter_audit = out_dfs["physiology_filter_audit"]
    regulation_table_out = out_dfs["regulation_table_out"].drop(drop_cols_final, strict=False)
    regulation_filter_audit = out_dfs["regulation_filter_audit"]

    if save_to_parquet:
        annot_table.write_parquet(all_annot_out_table)
        metabolism_table_out.write_parquet(out_metabolism_table)
        metabolism_filter_audit.write_parquet(out_metabolism_audit)
        physiology_table_out.write_parquet(out_physiology_table)
        physiology_filter_audit.write_parquet(out_physiology_audit)
        regulation_table_out.write_parquet(out_regulation_table)
        regulation_filter_audit.write_parquet(out_regulation_audit)
    else:
        annot_table.write_csv(all_annot_out_table, separator="\t")
        metabolism_table_out.write_csv(out_metabolism_table, separator="\t")
        metabolism_filter_audit.write_csv(out_metabolism_audit, separator="\t")
        physiology_table_out.write_csv(out_physiology_table, separator="\t")
        physiology_filter_audit.write_csv(out_physiology_audit, separator="\t")
        regulation_table_out.write_csv(out_regulation_table, separator="\t")
        regulation_filter_audit.write_csv(out_regulation_audit, separator="\t")

    logger.info(f"Total number of genes analyzed: {annot_table.shape[0]:,}")
    log_filter_summary(metabolism_filter_audit, "AMG", avg_array_len_limit)
    log_filter_summary(physiology_filter_audit, "APG", avg_array_len_limit)
    log_filter_summary(regulation_filter_audit, "AReG", avg_array_len_limit)
    logger.info(f"Final number of curated metabolic genes: {metabolism_table_out.shape[0]:,}")
    logger.info(f"Final number of curated physiology genes: {physiology_table_out.shape[0]:,}")
    logger.info(f"Final number of curated regulatory genes: {regulation_table_out.shape[0]:,}")
    logger.info("Curation of annotations completed.")


if __name__ == "__main__":
    main()
