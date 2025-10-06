#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
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

# Load KEGG and FOAM thresholds
KEGG_THRESHOLDS_PATH = snakemake.params.kegg_cutoff_file
if Path(KEGG_THRESHOLDS_PATH).exists():
    df = pl.read_csv(KEGG_THRESHOLDS_PATH)
    KEGG_THRESHOLDS = dict(zip(df["id"].to_list(), df["threshold"].to_list()))
FOAM_THRESHOLDS_PATH = snakemake.params.foam_cutoff_file
if Path(FOAM_THRESHOLDS_PATH).exists():
    df = pl.read_csv(FOAM_THRESHOLDS_PATH)
    FOAM_THRESHOLDS = dict(zip(df["id"].to_list(), df["cutoff_full"].to_list()))
    
def summarize_annot_table(table, hmm_descriptions):
    """
    Summarizes the table with gene annotations by selecting relevant columns,
    and merging with HMM descriptions.
    Returns: pl.DataFrame
    """
    table = table.with_columns([
        pl.min_horizontal([
            pl.col("KEGG_viral_left_dist"),
            pl.col("Pfam_viral_left_dist"),
            pl.col("PHROG_viral_left_dist")
        ]).alias("Viral_Flanking_Genes_Left_Dist"),
        pl.min_horizontal([
            pl.col("KEGG_viral_right_dist"),
            pl.col("Pfam_viral_right_dist"),
            pl.col("PHROG_viral_right_dist")
        ]).alias("Viral_Flanking_Genes_Right_Dist"),
        pl.min_horizontal([
            pl.col("KEGG_MGE_left_dist"),
            pl.col("Pfam_MGE_left_dist"),
            pl.col("PHROG_MGE_left_dist")
        ]).alias("MGE_Flanking_Genes_Left_Dist"),
        pl.min_horizontal([
            pl.col("KEGG_MGE_right_dist"),
            pl.col("Pfam_MGE_right_dist"),
            pl.col("PHROG_MGE_right_dist")
        ]).alias("MGE_Flanking_Genes_Right_Dist")
    ])
    # rest of your required_cols etc unchanged, except for new names
    required_cols = [
        "protein", "contig", "circular_contig", "genome", "gene_number",
        "KEGG_hmm_id", "FOAM_hmm_id", "Pfam_hmm_id", "dbCAN_hmm_id", "METABOLIC_hmm_id", "PHROG_hmm_id",
        "KEGG_score", "FOAM_score", "Pfam_score", "dbCAN_score", "METABOLIC_score", "PHROG_score",
        "KEGG_coverage", "FOAM_coverage", "Pfam_coverage", "dbCAN_coverage", "METABOLIC_coverage", "PHROG_coverage",
        "KEGG_V-score", "Pfam_V-score", "PHROG_V-score",
        "window_avg_KEGG_VL-score_viral", "window_avg_Pfam_VL-score_viral", "window_avg_PHROG_VL-score_viral",
        "Viral_Flanking_Genes_Left_Dist", "Viral_Flanking_Genes_Right_Dist",
        "MGE_Flanking_Genes_Left_Dist", "MGE_Flanking_Genes_Right_Dist",
        "Viral_Origin_Confidence"
    ]
    # fill missing
    for col in required_cols:
        if col not in table.columns:
            if col.endswith("_id"):
                dtype = pl.Utf8
            elif col.endswith("_score"):
                dtype = pl.Float64
            elif col.endswith("_coverage"):
                dtype = pl.Float64
            else:
                dtype = pl.Utf8
            table = table.with_columns(pl.lit(None, dtype=dtype).alias(col))
    table = table.select(required_cols)
    table = table.rename({"protein": "Protein"})

    # KEGG and FOAM joins (exact ID, all rows)
    table = table.join(hmm_descriptions, left_on="KEGG_hmm_id", right_on="id", how="left").rename({"name": "KEGG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    table = table.join(hmm_descriptions, left_on="FOAM_hmm_id", right_on="id", how="left").rename({"name": "FOAM_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # Pfam: join to normalized IDs
    table = table.with_columns([
        pl.col("Pfam_hmm_id").str.replace(r"\.\d+$", "", literal=False).alias("Pfam_hmm_id_norm")
    ])
    hmm_pfam = hmm_descriptions.filter(pl.col("db") == "Pfam").with_columns([
        pl.col("id").str.replace(r"\.\d+$", "", literal=False).alias("id_norm")
    ])
    table = table.join(hmm_pfam, left_on="Pfam_hmm_id_norm", right_on="id_norm", how="left").rename({"name": "Pfam_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    table = table.drop(["Pfam_hmm_id_norm", "id_norm"])

    # dbCAN: handle underscore normalization
    table = table.with_columns(pl.col("dbCAN_hmm_id").str.replace(r'_(.*)', '', literal=False).alias("dbCAN_hmm_id_no_underscore"))
    table = table.join(hmm_descriptions, left_on="dbCAN_hmm_id_no_underscore", right_on="id", how="left").rename({"name": "dbCAN_Description"})
    table = table.drop("dbCAN_hmm_id_no_underscore")
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # METABOLIC join
    table = table.drop(col for col in table.columns if col.endswith("_right"))
    table = table.join(hmm_descriptions, left_on="METABOLIC_hmm_id", right_on="id", how="left").rename({"name": "METABOLIC_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    # PHROG join
    table = table.drop(col for col in table.columns if col.endswith("_right"))
    table = table.join(hmm_descriptions, left_on="PHROG_hmm_id", right_on="id", how="left").rename({"name": "PHROG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")

    # Processing scores and coverage
    score_cols = ["KEGG_score", "FOAM_score", "Pfam_score", "dbCAN_score", "METABOLIC_score", "PHROG_score"]
    table = table.with_columns([
        pl.col(c).cast(pl.Float64).fill_null(float('-inf')).alias(c) for c in score_cols
    ])
    table = table.with_columns([
        pl.max_horizontal(score_cols).alias("max_score")
    ])
    table = table.with_columns(
        pl.when(pl.col("max_score").is_null())
        .then(None)
        .otherwise(
            pl.struct(score_cols).map_elements(
                lambda row: list(row.values()).index(max(row.values())),
                return_dtype=pl.Int64
            )
        ).alias("best_idx")
    )
    table = table.drop("max_score")
    table = table.with_columns([
        pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_hmm_id"))
        .when(pl.col("best_idx") == 1).then(pl.col("FOAM_hmm_id"))
        .when(pl.col("best_idx") == 2).then(pl.col("Pfam_hmm_id"))
        .when(pl.col("best_idx") == 3).then(pl.col("dbCAN_hmm_id"))
        .when(pl.col("best_idx") == 4).then(pl.col("METABOLIC_hmm_id"))
        .when(pl.col("best_idx") == 5).then(pl.col("PHROG_hmm_id"))
        .otherwise(pl.lit(None)).alias("top_hit_hmm_id"),

        pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_Description"))
        .when(pl.col("best_idx") == 1).then(pl.col("FOAM_Description"))
        .when(pl.col("best_idx") == 2).then(pl.col("Pfam_Description"))
        .when(pl.col("best_idx") == 3).then(pl.col("dbCAN_Description"))
        .when(pl.col("best_idx") == 4).then(pl.col("METABOLIC_Description"))
        .when(pl.col("best_idx") == 5).then(pl.col("PHROG_Description"))
        .otherwise(pl.lit(None)).alias("top_hit_description"),

        pl.when(pl.col("best_idx") == 0).then(pl.lit("KEGG"))
        .when(pl.col("best_idx") == 1).then(pl.lit("FOAM"))
        .when(pl.col("best_idx") == 2).then(pl.lit("Pfam"))
        .when(pl.col("best_idx") == 3).then(pl.lit("dbCAN"))
        .when(pl.col("best_idx") == 4).then(pl.lit("METABOLIC"))
        .when(pl.col("best_idx") == 5).then(pl.lit("PHROG"))
        .otherwise(pl.lit(None)).alias("top_hit_db"),
    ])
    table = table.drop(["best_idx"])

    # Final select/rename/output
    table = table.select([
        "Protein", "contig", "genome", "gene_number",
        "KEGG_V-score", "Pfam_V-score", "PHROG_V-score",
        "KEGG_hmm_id", "KEGG_Description", "KEGG_score", "KEGG_coverage",
        "FOAM_hmm_id", "FOAM_Description", "FOAM_score", "FOAM_coverage",
        "Pfam_hmm_id", "Pfam_Description", "Pfam_score", "Pfam_coverage",
        "dbCAN_hmm_id", "dbCAN_Description", "dbCAN_score", "dbCAN_coverage",
        "METABOLIC_hmm_id", "METABOLIC_Description", "METABOLIC_score", "METABOLIC_coverage",
        "PHROG_hmm_id", "PHROG_Description", "PHROG_score", "PHROG_coverage",
        "top_hit_hmm_id", "top_hit_description", "top_hit_db",
        "circular_contig", "Viral_Origin_Confidence",
        "Viral_Flanking_Genes_Left_Dist", "Viral_Flanking_Genes_Right_Dist",
        "MGE_Flanking_Genes_Left_Dist", "MGE_Flanking_Genes_Right_Dist"
    ])
    table = table.rename({
        "contig": "Contig",
        "genome": "Genome",
        "circular_contig": "Circular_Contig"
    })
    table = table.unique()
    return table.sort(["Genome", "Contig", "gene_number"])

def filter_false_substrings(table, false_substring_table, bypass_min_bitscore, bypass_min_cov, valid_hmm_ids, filter_presets):
    """
    Filter results to exclude false positives based on descriptions.
    - Special EC filters: distinguish between exact EC matches vs. class/subclass matches.
    - Special word-boundary filters for 'lysin' and 'ADP'
    Also returns an audit table with columns: removed, kept, remove_reason, keep_reason
    """

    # normalize presets
    if isinstance(filter_presets, str):
        presets = {p.strip().lower() for p in filter_presets.split(",") if p.strip()}
    else:
        presets = {str(p).strip().lower() for p in (filter_presets or [])}

    sources = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "PHROG"]
    desc_cols = [f"{src}_Description" for src in sources]
    specials = {"lysin", "adp"}

    def is_exact_ec(keyword):
        return bool(re.fullmatch(r"EC[:\s]\d+\.\d+\.\d+\.\d+", keyword))

    hard_meta_all = []
    soft_meta_all = []
    for kw, kw_type, kw_exc in zip(
        false_substring_table["substring"],
        false_substring_table["type"],
        false_substring_table["exception"],
    ):
        kw_s = str(kw)
        kw_lc = kw_s.strip().lower()
        exc = (kw_exc or "").strip().lower()
        if kw_lc in specials:
            pat = rf"(?i)\b{re.escape(kw_lc)}\b"
        elif is_exact_ec(kw_s):
            ec_number = kw_s.split()[-1] if " " in kw_s else kw_s.split(":")[-1]
            pat = rf"(?i)\bEC[:\s]{re.escape(ec_number)}\b"
        else:
            pat = rf"(?i){re.escape(kw_s)}"
        t = str(kw_type).strip().lower()
        if t == "hard":
            hard_meta_all.append((pat, kw_s, exc))
        elif t == "soft":
            soft_meta_all.append((pat, kw_s, exc))

    hard_meta_active = [(p, k, e) for (p, k, e) in hard_meta_all if not (e and e in presets)]
    if "no_soft_filter" in presets:
        soft_meta_active = []
    else:
        soft_meta_active = [(p, k, e) for (p, k, e) in soft_meta_all if not (e and e in presets)]

    def any_match_flag(meta_list):
        exprs = [pl.col(col).str.contains(pat, literal=False).fill_null(False)
                 for (pat, _, _) in meta_list for col in desc_cols]
        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    def first_matching_keyword(meta_list):
        exprs = []
        for (pat, kw, _) in meta_list:
            for col in desc_cols:
                exprs.append(
                    pl.when(pl.col(col).str.contains(pat, literal=False))
                    .then(pl.lit(kw))
                    .otherwise(pl.lit(None))
                )
        return pl.coalesce(exprs) if exprs else pl.lit(None, dtype=pl.Utf8)

    def exception_match_flag_and_token(meta_list):
        allowed = [(pat, kw, exc) for (pat, kw, exc) in meta_list if exc and exc in presets]
        if not allowed:
            return pl.lit(False), pl.lit(None, dtype=pl.Utf8), pl.lit(None, dtype=pl.Utf8)
        flag_exprs, kw_exprs, tok_exprs = [], [], []
        for (pat, kw, exc) in allowed:
            for col in desc_cols:
                m = pl.col(col).str.contains(pat, literal=False).fill_null(False)
                flag_exprs.append(m)
                kw_exprs.append(pl.when(m).then(pl.lit(kw)).otherwise(pl.lit(None)))
                tok_exprs.append(pl.when(m).then(pl.lit(exc)).otherwise(pl.lit(None)))
        return (pl.any_horizontal(flag_exprs) if flag_exprs else pl.lit(False),
                pl.coalesce(kw_exprs) if kw_exprs else pl.lit(None, dtype=pl.Utf8),
                pl.coalesce(tok_exprs) if tok_exprs else pl.lit(None, dtype=pl.Utf8))

    def build_soft_valid(meta_list, hmm_ids_in_group):
        exprs = []
        for (pat, _, _) in meta_list:
            for src in sources:
                desc_col = f"{src}_Description"
                score_col = f"{src}_score"
                cov_col = f"{src}_coverage"
                id_col = f"{src}_hmm_id"

                score_col_casted = pl.col(score_col).cast(pl.Float64).fill_null(float("-inf"))
                cov_col_casted = pl.col(cov_col).cast(pl.Float64).fill_null(float("-inf"))

                if src == "KEGG":
                    score_thresh_expr = (
                        pl.when(pl.col(id_col).is_in(list(KEGG_THRESHOLDS.keys())))
                        .then(pl.col(id_col).map_elements(lambda x: KEGG_THRESHOLDS.get(x, bypass_min_bitscore), return_dtype=pl.Float64))
                        .otherwise(pl.lit(bypass_min_bitscore))
                    )
                elif src == "FOAM":
                    score_thresh_expr = (
                        pl.when(pl.col(id_col).is_in(list(FOAM_THRESHOLDS.keys())))
                        .then(pl.col(id_col).map_elements(lambda x: FOAM_THRESHOLDS.get(x, bypass_min_bitscore), return_dtype=pl.Float64))
                        .otherwise(pl.lit(bypass_min_bitscore))
                    )
                else:
                    score_thresh_expr = pl.lit(bypass_min_bitscore)

                exprs.append(
                    pl.col(desc_col).str.contains(pat, literal=False).fill_null(False) &
                    pl.col(id_col).is_in(hmm_ids_in_group) &
                    score_col_casted.is_finite() & (score_col_casted >= score_thresh_expr) &
                    cov_col_casted.is_finite() & (cov_col_casted >= bypass_min_cov)
                )
        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    # active flags for actual filtering
    hard_match_active = any_match_flag(hard_meta_active).alias("HARD_MATCH")
    soft_match_active = any_match_flag(soft_meta_active).alias("SOFT_MATCH")
    soft_valid_active = build_soft_valid(soft_meta_active, valid_hmm_ids).alias("SOFT_VALID")

    # all-patterns audit helpers
    hard_match_all = any_match_flag(hard_meta_all).alias("HARD_MATCH_ALL")
    soft_match_all = any_match_flag(soft_meta_all).alias("SOFT_MATCH_ALL")
    soft_valid_all = build_soft_valid(soft_meta_all, valid_hmm_ids).alias("SOFT_VALID_ALL")
    hard_first_kw_all = first_matching_keyword(hard_meta_all).alias("HARD_REMOVE_KEYWORD_ALL")
    soft_first_kw_all = first_matching_keyword(soft_meta_all).alias("SOFT_REMOVE_KEYWORD_ALL")
    hard_exc_flag_all, _, hard_exc_tok_all = exception_match_flag_and_token(hard_meta_all)
    soft_exc_flag_all, _, soft_exc_tok_all = exception_match_flag_and_token(soft_meta_all)
    hard_exc_flag_all = hard_exc_flag_all.alias("HARD_EXCEPTION_ALL")
    soft_exc_flag_all = soft_exc_flag_all.alias("SOFT_EXCEPTION_ALL")
    hard_exc_tok_all = hard_exc_tok_all.alias("HARD_EXCEPTION_TOKEN_ALL")
    soft_exc_tok_all = soft_exc_tok_all.alias("SOFT_EXCEPTION_TOKEN_ALL")

    table_with_flags = table.with_columns([
        hard_match_active, soft_match_active, soft_valid_active,
        hard_match_all, soft_match_all, soft_valid_all,
        hard_first_kw_all, soft_first_kw_all,
        hard_exc_flag_all, soft_exc_flag_all,
        hard_exc_tok_all, soft_exc_tok_all,
    ])

    # actual filter unchanged
    if "no_filter" in presets:
        table_filtered = table
    else:
        table_filtered = table_with_flags.filter(
            (~pl.col("HARD_MATCH")) & ((~pl.col("SOFT_MATCH")) | pl.col("SOFT_VALID"))
        ).drop(["HARD_MATCH", "SOFT_MATCH", "SOFT_VALID"])

    # removal decision from active flags only
    if "no_filter" in presets:
        removed_expr = pl.lit(False)
    else:
        hard_remove_expr = pl.col("HARD_MATCH")
        soft_remove_expr = pl.col("SOFT_MATCH") & (~pl.col("SOFT_VALID"))
        removed_expr = hard_remove_expr | soft_remove_expr

    kept_expr = ~removed_expr

    # remove_reason tied to actual cause only (no mixing)
    remove_reason_expr = (
        pl.when(removed_expr & pl.col("HARD_MATCH"))
        .then(pl.col("HARD_REMOVE_KEYWORD_ALL"))
        .when(removed_expr & (pl.col("SOFT_MATCH") & (~pl.col("SOFT_VALID"))))
        .then(pl.col("SOFT_REMOVE_KEYWORD_ALL"))
        .otherwise(pl.lit(None))
    )

    # hypothetical remove_reason for kept rows (first real would-be cause)
    kept_side_remove_reason = (
        pl.when(pl.col("HARD_MATCH_ALL"))
        .then(pl.col("HARD_REMOVE_KEYWORD_ALL"))
        .when(pl.col("SOFT_MATCH_ALL") & (~pl.col("SOFT_VALID_ALL")))
        .then(pl.col("SOFT_REMOVE_KEYWORD_ALL"))
        .otherwise(pl.lit(None))
    )
    remove_reason = pl.when(kept_expr).then(kept_side_remove_reason).otherwise(remove_reason_expr)

    # keep_reason with strict precedence, applied only when kept
    # 1) no_filter
    # 2) no_soft_filter, but only if a soft removal would have happened
    # 3) exception:<token> for hard (must have hard match)
    # 4) exception:<token> for soft (must have soft match that would remove)
    # 5) soft filter minimum HMM bitscore and coverage met (soft match and thresholds met)
    keep_no_filter = pl.when(kept_expr & pl.lit("no_filter" in presets)).then(pl.lit("no_filter")).otherwise(pl.lit(None))
    keep_no_soft = pl.when(kept_expr & pl.lit("no_soft_filter" in presets) & (pl.col("SOFT_MATCH_ALL") & (~pl.col("SOFT_VALID_ALL")))).then(pl.lit("no_soft_filter")).otherwise(pl.lit(None))
    keep_exc_hard = pl.when(kept_expr & pl.col("HARD_EXCEPTION_ALL") & pl.col("HARD_MATCH_ALL")).then(pl.concat_str([pl.lit("exception:"), pl.col("HARD_EXCEPTION_TOKEN_ALL")])).otherwise(pl.lit(None))
    keep_exc_soft = pl.when(kept_expr & pl.col("SOFT_EXCEPTION_ALL") & (pl.col("SOFT_MATCH_ALL") & (~pl.col("SOFT_VALID_ALL")))).then(pl.concat_str([pl.lit("exception:"), pl.col("SOFT_EXCEPTION_TOKEN_ALL")])).otherwise(pl.lit(None))
    keep_soft_valid = pl.when(kept_expr & (pl.col("SOFT_MATCH_ALL") & pl.col("SOFT_VALID_ALL"))).then(pl.lit("soft filter minimum HMM bitscore and coverage met")).otherwise(pl.lit(None))

    # precedence: build a single reason string without mixing unrelated reasons
    keep_reason_expr = (
        pl.when(keep_no_filter.is_not_null()).then(keep_no_filter)
        .when(keep_no_soft.is_not_null()).then(keep_no_soft)
        .when(keep_exc_hard.is_not_null()).then(keep_exc_hard)
        .when(keep_exc_soft.is_not_null()).then(keep_exc_soft)
        .when(keep_soft_valid.is_not_null()).then(keep_soft_valid)
        .otherwise(pl.lit(None))
    )

    # if removed, keep_reason must be None (prevents mixing)
    keep_reason_expr = pl.when(removed_expr).then(pl.lit(None)).otherwise(keep_reason_expr)

    audit = (
        table_with_flags
        .with_columns([
            removed_expr.alias("removed"),
            kept_expr.alias("kept"),
            remove_reason.alias("remove_reason"),
            keep_reason_expr.alias("keep_reason"),
        ])
        .drop([
            "HARD_MATCH", "SOFT_MATCH", "SOFT_VALID",
            "HARD_MATCH_ALL", "SOFT_MATCH_ALL", "SOFT_VALID_ALL",
            "HARD_REMOVE_KEYWORD_ALL", "SOFT_REMOVE_KEYWORD_ALL",
            "HARD_EXCEPTION_ALL", "SOFT_EXCEPTION_ALL",
            "HARD_EXCEPTION_TOKEN_ALL", "SOFT_EXCEPTION_TOKEN_ALL",
        ])
    )

    return table_filtered, audit

def filter_metabolism_annots(table, metabolism_table, false_substring_table, bypass_min_bitscore, bypass_min_cov, filter_presets):
    """
    Identify metabolism-related genes based on input metabolism table.
    by checking any of the five HMM ID columns for membership in metabolism_table["id"].
    Also, apply false-substring filtering to remove non-metabolic genes.
    """
    metab_ids = metabolism_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(metab_ids) |
        pl.col("FOAM_hmm_id").is_in(metab_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(metab_ids) |
        pl.col("dbCAN_hmm_id").is_in(metab_ids) |
        pl.col("METABOLIC_hmm_id").is_in(metab_ids) |
        pl.col("PHROG_hmm_id").is_in(metab_ids)
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table, audit = filter_false_substrings(table, false_substring_table, bypass_min_bitscore, bypass_min_cov, metab_ids, filter_presets)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])

def filter_physiology_annots(table, physiology_table, false_phys_substrings, bypass_min_bitscore, bypass_min_cov, filter_presets):
    """
    Identify physiology-related genes based on input physiology table.
    by checking any of the five HMM ID columns for membership in physiology_table["id"].
    Also, apply false-substring filtering to remove non-physiological genes.
    """
    phys_ids = physiology_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(phys_ids) |
        pl.col("FOAM_hmm_id").is_in(phys_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(phys_ids) |
        pl.col("dbCAN_hmm_id").is_in(phys_ids) |
        pl.col("METABOLIC_hmm_id").is_in(phys_ids) |
        pl.col("PHROG_hmm_id").is_in(phys_ids)
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table, audit = filter_false_substrings(table, false_phys_substrings, bypass_min_bitscore, bypass_min_cov, phys_ids, filter_presets)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "gene_number"]), audit.sort(["Genome", "Contig", "gene_number"])

def filter_regulation_annots(table, regulation_table, false_reg_substrings, bypass_min_bitscore, bypass_min_cov, filter_presets):
    """
    Identify regulation-related genes based on input regulation table.
    by checking any of the five HMM ID columns for membership in regulation_table["id"].
    Also, apply false-substring filtering to remove non-regulatory genes.
    """
    reg_ids = regulation_table["id"].to_list()
    condition = (
        pl.col("KEGG_hmm_id").is_in(reg_ids) |
        pl.col("FOAM_hmm_id").is_in(reg_ids) |
        pl.col("Pfam_hmm_id_clean").is_in(reg_ids) |
        pl.col("dbCAN_hmm_id").is_in(reg_ids) |
        pl.col("METABOLIC_hmm_id").is_in(reg_ids) |
        pl.col("PHROG_hmm_id").is_in(reg_ids)
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table, audit = filter_false_substrings(table, false_reg_substrings, bypass_min_bitscore, bypass_min_cov, reg_ids, filter_presets)
    
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
    false_metab_substrings = snakemake.params.false_amgs
    false_phys_substrings = snakemake.params.false_apgs
    false_reg_substrings = snakemake.params.false_aregs
    filter_presets = list(snakemake.params.filter_presets.split(","))
    
    scaling_factor = max(snakemake.params.soft_keyword_bypass_scaling_factor, 1.0)
    bypass_min_bitscore = float(scaling_factor * snakemake.params.min_bitscore)
    bypass_min_cov = min(float(scaling_factor * snakemake.params.cov_fraction), 1.0)
    
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

    hmm_descriptions = pl.read_csv(hmm_ref, schema={"id": pl.Utf8, "db": pl.Utf8, "name": pl.Utf8})
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
    
    false_metab_substring_table = pl.read_csv(false_metab_substrings)
    false_phys_substring_table = pl.read_csv(false_phys_substrings)
    false_reg_substring_table = pl.read_csv(false_reg_substrings)
    
    annot_table = summarize_annot_table(table, hmm_descriptions)
    
    # Remove .X or .XX suffixes from top_hit_hmm_id for proper matching of Pfam hits
    annot_table = annot_table.with_columns(
        pl.col("top_hit_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("top_hit_hmm_id_clean"),
        pl.col("Pfam_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("Pfam_hmm_id_clean"),
    )
    
    drop_cols = ["gene_number", "window_avg_KEGG_VL-score_viral", "window_avg_Pfam_VL-score_viral", "window_avg_PHROG_VL-score_viral", "top_hit_hmm_id_clean", "Pfam_hmm_id_clean"]
    metabolism_table_out, metabolism_filter_audit = filter_metabolism_annots(annot_table, metabolism_table, false_metab_substring_table, bypass_min_bitscore, bypass_min_cov, filter_presets)
    physiology_table_out, physiology_filter_audit = filter_physiology_annots(annot_table, physiology_table, false_phys_substring_table, bypass_min_bitscore, bypass_min_cov, filter_presets)
    regulation_table_out, regulation_filter_audit = filter_regulation_annots(annot_table, regulation_table, false_reg_substring_table, bypass_min_bitscore, bypass_min_cov, filter_presets)
    out_dfs = {
        "annot_table": annot_table,
        "metabolism_table_out": metabolism_table_out,
        "metabolism_filter_audit": metabolism_filter_audit,
        "physiology_table_out": physiology_table_out,
        "physiology_filter_audit": physiology_filter_audit,
        "regulation_table_out": regulation_table_out,
        "regulation_filter_audit": regulation_filter_audit,
    }
    for table in out_dfs.keys():
        df = out_dfs[table]
        df = df.drop([col for col in df.columns if col in drop_cols])
        replacements = []
        extra_drop_cols = []
        for col in df.columns:
            if "audit" in table:
                if not col.endswith("_Description") and not col.endswith("_hmm_id") and not col.endswith("_coverage") and not col.endswith("_score") and not col in ["Protein", "Contig", "Genome", "removed", "kept", "remove_reason", "keep_reason"]:
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
        out_dfs[table] = df
        
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