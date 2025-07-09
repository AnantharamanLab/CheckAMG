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
    
    Args:
        input_table (pl.DataFrame): Polars DataFrame of genome context.
        hmm_descriptions (pl.DataFrame): Polars DataFrame containing a mapping of
        HMM names to their descriptions.
    
    Returns:
        pl.DataFrame: A summarized Polars DataFrame.
    """
    # Add logical columns
    table = table.with_columns(
        pl.when(
            pl.col("KEGG_verified_flank_up") | pl.col("Pfam_verified_flank_up") | pl.col("PHROG_verified_flank_up")
        ).then(pl.lit(True))
        .otherwise(pl.lit(False)).alias("Viral_Flanking_Genes_Upstream"),
        pl.when(
            pl.col("KEGG_verified_flank_down") | pl.col("Pfam_verified_flank_down") | pl.col("PHROG_verified_flank_down")
        ).then(pl.lit(True))
        .otherwise(pl.lit(False)).alias("Viral_Flanking_Genes_Downstream"),
        pl.when(
            pl.col("KEGG_MGE_flank") | pl.col("Pfam_MGE_flank") | pl.col("PHROG_MGE_flank")
        ).then(pl.lit(True))
        .otherwise(pl.lit(False)).alias("MGE_Flanking_Genes")
    )

    required_cols = [
        "protein", "contig", "circular_contig", "genome", "gene_number",
        "KEGG_hmm_id", "FOAM_hmm_id", "Pfam_hmm_id", "dbCAN_hmm_id", "METABOLIC_hmm_id", "PHROG_hmm_id",
        "KEGG_score", "FOAM_score", "Pfam_score", "dbCAN_score", "METABOLIC_score", "PHROG_score",
        "KEGG_coverage", "FOAM_coverage", "Pfam_coverage", "dbCAN_coverage", "METABOLIC_coverage", "PHROG_coverage",
        "KEGG_V-score", "Pfam_V-score", "PHROG_V-score",
        "window_avg_KEGG_VL-score_viral", "window_avg_Pfam_VL-score_viral", "window_avg_PHROG_VL-score_viral",
        "KEGG_verified_flank_up", "KEGG_verified_flank_down", "Pfam_verified_flank_up", "Pfam_verified_flank_down", "PHROG_verified_flank_up", "PHROG_verified_flank_down",
        "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream",
        "KEGG_MGE_flank", "Pfam_MGE_flank", "PHROG_MGE_flank", "MGE_Flanking_Genes",
        "Viral_Origin_Confidence"
    ]
    for col in required_cols:
        if col not in table.columns:
            if col.endswith("_id"):
                dtype = pl.Utf8
            elif col.endswith("_score"):
                dtype = pl.Float64
            elif col.endswith("_coverage"):
                dtype = pl.Float64
            elif "_VL-score_viral" in col:
                dtype = pl.Boolean
            elif col.endswith("_verified_flank") or col.endswith("_MGE_flank"):
                dtype = pl.Boolean
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
        "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes"
    ])
    table = table.rename({
        "contig": "Contig",
        "genome": "Genome",
        "circular_contig": "Circular_Contig"
    })
    table = table.unique()
    return table.sort(["Genome", "Contig", "gene_number"])

def filter_false_substrings(table, false_substring_table, bypass_min_bitscore, bypass_min_cov, valid_hmm_ids):
    """
    Filter results to exclude false positives based on descriptions.
    - Special EC filters: distinguish between exact EC matches vs. class/subclass matches.
    - Special word-boundary filters for 'lysin' and 'ADP'
    """
    sources = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "PHROG"]
    desc_cols = [f"{src}_Description" for src in sources]

    specials = {"lysin", "adp"}

    def is_exact_ec(keyword):
        return bool(re.fullmatch(r"EC[:\s]\d+\.\d+\.\d+\.\d+", keyword))

    hard_patterns = []
    soft_patterns = []
    for kw, kw_type in zip(false_substring_table["substring"], false_substring_table["type"]):
        kw_lc = kw.lower()
        if kw_lc in specials:
            pattern = rf"(?i)\b{re.escape(kw_lc)}\b"
        elif is_exact_ec(kw):
            ec_number = kw.split()[-1] if " " in kw else kw.split(":")[-1]
            pattern = rf"(?i)\bEC[:\s]{re.escape(ec_number)}\b"
        else:
            pattern = rf"(?i){re.escape(kw)}"

        if kw_type.lower() == "hard":
            hard_patterns.append(pattern)
        elif kw_type.lower() == "soft":
            soft_patterns.append(pattern)

    def build_pattern_flags(patterns, desc_cols):
        exprs = [
            pl.col(col).str.contains(pat, literal=False).fill_null(False)
            for pat in patterns for col in desc_cols
        ]
        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    def build_soft_valid_flags(patterns, table, hmm_ids_in_group):
        exprs = []
        for pat in patterns:
            for src in ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "PHROG"]:
                desc_col = f"{src}_Description"
                score_col = f"{src}_score"
                cov_col = f"{src}_coverage"
                id_col = f"{src}_hmm_id"

                score_col_casted = pl.col(score_col).cast(pl.Float64).fill_null(float("-inf"))
                cov_col_casted = pl.col(cov_col).cast(pl.Float64).fill_null(float("-inf"))

                # Source-specific threshold expression
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
                    score_col_casted.is_finite() &
                    (score_col_casted >= score_thresh_expr) &
                    cov_col_casted.is_finite() &
                    (cov_col_casted >= bypass_min_cov)
                )

        return pl.any_horizontal(exprs) if exprs else pl.lit(False)

    # Build flags
    table = table.with_columns([
        build_pattern_flags(hard_patterns, desc_cols).alias("HARD_MATCH"),
        build_pattern_flags(soft_patterns, desc_cols).alias("SOFT_MATCH"),
        build_soft_valid_flags(soft_patterns, table, valid_hmm_ids).alias("SOFT_VALID"),
    ])

    table_filtered = table.filter(
        (~pl.col("HARD_MATCH")) & ((~pl.col("SOFT_MATCH")) | pl.col("SOFT_VALID"))
    )

    return table_filtered.drop(["HARD_MATCH", "SOFT_MATCH", "SOFT_VALID"])

def filter_metabolism_annots(table, metabolism_table, false_substring_table, bypass_min_bitscore, bypass_min_cov):
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
    table = filter_false_substrings(table, false_substring_table, bypass_min_bitscore, bypass_min_cov, metab_ids)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "gene_number"])

def filter_physiology_annots(table, physiology_table, false_phys_substrings, bypass_min_bitscore, bypass_min_cov):
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
    table = filter_false_substrings(table, false_phys_substrings, bypass_min_bitscore, bypass_min_cov, phys_ids)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "gene_number"])

def filter_regulation_annots(table, regulation_table, false_reg_substrings, bypass_min_bitscore, bypass_min_cov):
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
    table = filter_false_substrings(table, false_reg_substrings, bypass_min_bitscore, bypass_min_cov, reg_ids)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "gene_number"])

def main():
    input_table  = snakemake.params.context_table
    hmm_ref = snakemake.params.hmm_ref
    metabolism_ref = snakemake.params.metabolism_table
    physiology_ref = snakemake.params.physiology_table
    regulation_ref = snakemake.params.regulation_table
    false_metab_substrings = snakemake.params.false_amgs
    false_phys_substrings = snakemake.params.false_apgs
    false_reg_substrings = snakemake.params.false_aregs
    
    scaling_factor = snakemake.params.soft_keyword_bypass_scaling_factor
    bypass_min_bitscore = float(scaling_factor * snakemake.params.min_bitscore)
    bypass_min_cov = min(float(scaling_factor * snakemake.params.cov_fraction), 1.0)
    
    out_metabolism_table = snakemake.params.metabolism_table_out
    out_physiology_table = snakemake.params.physiology_table_out
    out_regulation_table = snakemake.params.regulation_table_out
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
    out_dfs = {
        "annot_table": annot_table,
        "metabolism_table_out": filter_metabolism_annots(annot_table, metabolism_table, false_metab_substring_table, bypass_min_bitscore, bypass_min_cov),
        "physiology_table_out": filter_physiology_annots(annot_table, physiology_table, false_phys_substring_table, bypass_min_bitscore, bypass_min_cov),
        "regulation_table_out": filter_regulation_annots(annot_table, regulation_table, false_reg_substring_table, bypass_min_bitscore, bypass_min_cov)
    }
    for table in out_dfs.keys():
        df = out_dfs[table]
        df = df.drop([col for col in df.columns if col in drop_cols])
        replacements = []
        for col in df.columns:
            if col.endswith("_score"):
                replacements.append(
                    pl.when(pl.col(col) == -float("inf"))
                    .then(None)
                    .otherwise(pl.col(col))
                    .alias(col)
                )
        if replacements:
            df = df.with_columns(replacements)
        out_dfs[table] = df
        
    annot_table, metabolism_table_out, physiology_table_out, regulation_table_out = out_dfs["annot_table"], out_dfs["metabolism_table_out"], out_dfs["physiology_table_out"], out_dfs["regulation_table_out"]

    annot_table.write_csv(all_annot_out_table, separator="\t")
    metabolism_table_out.write_csv(out_metabolism_table, separator="\t")
    physiology_table_out.write_csv(out_physiology_table, separator="\t")
    regulation_table_out.write_csv(out_regulation_table, separator="\t")
    
    logger.info("Curation of annotations completed.")
    logger.info(f"Total number of genes analyzed: {annot_table.shape[0]:,}")
    logger.info(f"Number of curated metabolic genes: {metabolism_table_out.shape[0]:,}")
    logger.info(f"Number of curated physiology genes: {physiology_table_out.shape[0]:,}")
    logger.info(f"Number of curated regulatory genes: {regulation_table_out.shape[0]:,}")

if __name__ == "__main__":
    main()