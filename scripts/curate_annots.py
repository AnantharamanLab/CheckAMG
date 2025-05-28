#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

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

if snakemake.params.build_or_annotate =="build":
    print("========================================================================\n  Step 9/22: Curate the predicted functions based on genomic context   \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n  Step 9/22: Curate the predicted functions based on genomic context   \n========================================================================\n")
elif snakemake.params.build_or_annotate == "annotate":
    print("========================================================================\n   Step 9/11: Curate the predicted functions based on genomic context   \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n   Step 9/11: Curate the predicted functions based on genomic context   \n========================================================================\n")

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
    required_cols = [
        "protein",
        "contig",
        "circular_contig",
        "genome",
        "KEGG_hmm_id",
        "Pfam_hmm_id",
        "dbCAN_hmm_id",
        "METABOLIC_hmm_id",
        "PHROG_hmm_id",
        "KEGG_score",
        "Pfam_score",
        "dbCAN_score",
        "METABOLIC_score",
        "PHROG_score",
        "KEGG_V-score",
        "Pfam_V-score",
        "PHROG_V-score",
        "window_avg_KEGG_VL-score_viral",
        "window_avg_Pfam_VL-score_viral",
        "window_avg_PHROG_VL-score_viral",
        "KEGG_verified_flank_up",
        "KEGG_verified_flank_down",
        "Pfam_verified_flank_up",
        "Pfam_verified_flank_down",
        "PHROG_verified_flank_up",
        "PHROG_verified_flank_down",
        "KEGG_MGE_flank",
        "Pfam_MGE_flank",
        "PHROG_MGE_flank"
    ]
    for col in required_cols:
        if col not in table.columns:
            if col.endswith("_id"):
                dtype = pl.Utf8
            elif col.endswith("_score"):
                dtype = pl.Float64
            elif "_VL-score_viral" in col:
                dtype = pl.Boolean
            elif col.endswith("_verified_flank") or col.endswith("_MGE_flank"):
                dtype = pl.Boolean
            else:
                dtype = pl.Utf8 # default to string type
            table = table.with_columns(pl.lit(None, dtype=dtype).alias(col))
        
    table = table.select(required_cols)
    table = table.rename({"protein": "Protein"})
    
    # Join with KEGG descriptions
    table = table.join(hmm_descriptions, left_on="KEGG_hmm_id", right_on="id", how="left").rename({"name": "KEGG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    
    # Join with Pfam descriptions
    table = table.join(hmm_descriptions, left_on="Pfam_hmm_id", right_on="id", how="left").rename({"name": "Pfam_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    
    # Custom join with dbCAN descriptions to handle IDs with underscores
    table = table.with_columns(pl.col("dbCAN_hmm_id").str.replace(r'_(.*)', '', literal=False).alias("dbCAN_hmm_id_no_underscore"))
    table = table.join(hmm_descriptions, left_on="dbCAN_hmm_id_no_underscore", right_on="id", how="left").rename({"name": "dbCAN_Description"})
    table = table.drop("dbCAN_hmm_id_no_underscore")
    if "db_right" in table.columns:
        table = table.drop("db_right")
    
    # Join with METABOLIC descriptions
    table = table.join(hmm_descriptions, left_on="METABOLIC_hmm_id", right_on="id", how="left").rename({"name": "METABOLIC_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
        
    # Join with PHROG descriptions
    table = table.join(hmm_descriptions, left_on="PHROG_hmm_id", right_on="id", how="left").rename({"name": "PHROG_Description"})
    if "db_right" in table.columns:
        table = table.drop("db_right")
    
    # Mark genes within virus-like windows (KEGG, Pfam, or PHROG)
    table = table.with_columns(
        pl.when(pl.col("window_avg_KEGG_VL-score_viral") | pl.col("window_avg_Pfam_VL-score_viral") | pl.col("window_avg_PHROG_VL-score_viral"))
        .then(True)
        .otherwise(False)
        .alias("Virus_Like_Window")
    ])

    
    # Mark genes with viral flanking genes
    table = table.with_columns(
        pl.when(pl.col("KEGG_verified_flank_up") | pl.col("Pfam_verified_flank_up") | pl.col("PHROG_verified_flank_up"))
        .then(True)
        .otherwise(False)
        .alias("Viral_Flanking_Genes_Upstream")
    )
    table = table.with_columns(
        pl.when(pl.col("KEGG_verified_flank_down") | pl.col("Pfam_verified_flank_down") | pl.col("PHROG_verified_flank_down"))
        .then(True)
        .otherwise(False)
        .alias("Viral_Flanking_Genes_Downstream")
    )
    
    # Mark genes with MGE flanking genes
    table = table.with_columns(
        pl.when(pl.col("KEGG_MGE_flank") | pl.col("Pfam_MGE_flank") | pl.col("PHROG_MGE_flank"))
        .then(True)
        .otherwise(False)
        .alias("MGE_Flanking_Genes")
    )
    
    # Infer best HMM hits based on bit scores
    ## 1. Cast all relevant scores to Float64, replacing null with -inf
    table = (
        table
        .with_columns([
            pl.col("KEGG_score").cast(pl.Float64).fill_null(float('-inf')),
            pl.col("Pfam_score").cast(pl.Float64).fill_null(float('-inf')),
            pl.col("dbCAN_score").cast(pl.Float64).fill_null(float('-inf')),
            pl.col("METABOLIC_score").cast(pl.Float64).fill_null(float('-inf')),
            pl.col("PHROG_score").cast(pl.Float64).fill_null(float('-inf')),
        ])
    )

    # 3. Cast scores to Float64, but *also* compute the real max before filling nulls
    score_cols = [
        "KEGG_score", "Pfam_score",
        "dbCAN_score", "METABOLIC_score", "PHROG_score",
    ]
    table = (
        table
        .with_columns([
            # the real max, ignoring nulls entirely
            pl.max_horizontal(score_cols).alias("max_score"),
            # now fill nulls to index into them safely
            *(pl.col(c).cast(pl.Float64).fill_null(float("-inf")).alias(c) for c in score_cols)
        ])
    )

    # 4. Build best_idx but leave it null if no real scores existed
    table = table.with_columns(
        pl.when(pl.col("max_score").is_null())
        .then(None) # if all scores were null
        .otherwise(
            pl.struct(score_cols).map_elements(
                lambda row: list(row.values()).index(max(row.values())),
                return_dtype=pl.Int64
            )
        )
        .alias("best_idx")
    )

    table = table.drop("max_score")

    ## 4. Use `best_idx` to fill in the top_hit_* columns
    table = table.with_columns([
        # top_hit_hmm_id
        pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_hmm_id"))
        .when(pl.col("best_idx") == 1).then(pl.col("Pfam_hmm_id"))
        .when(pl.col("best_idx") == 2).then(pl.col("dbCAN_hmm_id"))
        .when(pl.col("best_idx") == 3).then(pl.col("METABOLIC_hmm_id"))
        .when(pl.col("best_idx") == 4).then(pl.col("PHROG_hmm_id"))
        .otherwise(pl.lit(None))
        .alias("top_hit_hmm_id"),

        # top_hit_description
        pl.when(pl.col("best_idx") == 0).then(pl.col("KEGG_Description"))
        .when(pl.col("best_idx") == 1).then(pl.col("Pfam_Description"))
        .when(pl.col("best_idx") == 2).then(pl.col("dbCAN_Description"))
        .when(pl.col("best_idx") == 3).then(pl.col("METABOLIC_Description"))
        .when(pl.col("best_idx") == 4).then(pl.col("PHROG_Description"))
        .otherwise(pl.lit(None))
        .alias("top_hit_description"),

        # top_hit_db
        pl.when(pl.col("best_idx") == 0).then(pl.lit("KEGG"))
        .when(pl.col("best_idx") == 1).then(pl.lit("Pfam"))
        .when(pl.col("best_idx") == 2).then(pl.lit("dbCAN"))
        .when(pl.col("best_idx") == 3).then(pl.lit("METABOLIC"))
        .when(pl.col("best_idx") == 4).then(pl.lit("PHROG"))
        .otherwise(pl.lit(None))
        .alias("top_hit_db"),
    ])

    ## 5. Finally, remove the helper columns
    table = table.drop(["best_idx"])
    
    # Select only relevant columns for output
    table = table.select([
        "Protein",
        "contig",
        "genome",
        "KEGG_V-score",
        "Pfam_V-score",
        "PHROG_V-score",
        "KEGG_hmm_id",
        "KEGG_Description",
        "Pfam_hmm_id",
        "Pfam_Description",
        "dbCAN_hmm_id",
        "dbCAN_Description",
        "METABOLIC_hmm_id",
        "METABOLIC_Description",
        "PHROG_hmm_id",
        "PHROG_Description",
        "top_hit_hmm_id",
        "top_hit_description",
        "top_hit_db",
        "circular_contig",
        "Virus_Like_Window",
        "Viral_Flanking_Genes_Upstream",
        "Viral_Flanking_Genes_Downstream",
        "MGE_Flanking_Genes"
    ])
    table = table.rename({
        "contig": "Contig",
        "genome": "Genome",
        "circular_contig": "Circular_Contig"
    })
    
    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])

def filter_false_substrings(table, false_substrings_desc):
    """
    Filter results to exclude false positives based on descriptions.
    """
    for substring in false_substrings_desc:
        table = table.filter(
            ~(
                pl.col("KEGG_Description").str.contains(substring).fill_null(False) |
                pl.col("Pfam_Description").str.contains(substring).fill_null(False) |
                pl.col("dbCAN_Description").str.contains(substring).fill_null(False) |
                pl.col("METABOLIC_Description").str.contains(substring).fill_null(False) |
                pl.col("PHROG_Description").str.contains(substring).fill_null(False)
            )
        )
    return table

def filter_metabolism_annots(table, metabolism_table, false_metab_substrings):
    """
    Identify metabolism-related genes based on input metabolism table.
    by checking any of the five HMM ID columns for membership in metabolism_table["id"].
    Also, apply false-substring filtering to remove non-metabolic genes.
    """
    condition = (
        pl.col("KEGG_hmm_id").is_in(metabolism_table["id"]) |
        pl.col("Pfam_hmm_id_clean").is_in(metabolism_table["id"]) |
        pl.col("dbCAN_hmm_id").is_in(metabolism_table["id"]) |
        pl.col("METABOLIC_hmm_id").is_in(metabolism_table["id"]) |
        pl.col("PHROG_hmm_id").is_in(metabolism_table["id"])
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table = filter_false_substrings(table, false_metab_substrings)
    
    # dcCAN/CAZyme-spescific false AMGs 
    false_substrings_dbcan_id = [
        "CBM", # Carbohydrate-binding modules usually match to tail peptides used for receptor binding
        "GT", # Glycosyltransferases
        "GH", # Glycoside hydrolases
    ]   
    for substring in false_substrings_dbcan_id:
        table = table.filter(~(pl.col("dbCAN_hmm_id").str.contains(substring).fill_null(False)))
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])

def filter_physiology_annots(table, physiology_table, false_phys_substrings):
    """
    Identify physiology-related genes based on input physiology table.
    by checking any of the five HMM ID columns for membership in physiology_table["id"].
    Also, apply false-substring filtering to remove non-physiological genes.
    """
    condition = (
        pl.col("KEGG_hmm_id").is_in(physiology_table["id"]) |
        pl.col("Pfam_hmm_id_clean").is_in(physiology_table["id"]) |
        pl.col("dbCAN_hmm_id").is_in(physiology_table["id"]) |
        pl.col("METABOLIC_hmm_id").is_in(physiology_table["id"]) |
        pl.col("PHROG_hmm_id").is_in(physiology_table["id"])
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table = filter_false_substrings(table, false_phys_substrings)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])


def filter_regulation_annots(table, regulation_table, false_reg_substrings):
    """
    Identify regulation-related genes based on input regulation table.
    by checking any of the five HMM ID columns for membership in regulation_table["id"].
    Also, apply false-substring filtering to remove non-regulatory genes.
    """
    condition = (
        pl.col("KEGG_hmm_id").is_in(regulation_table["id"]) |
        pl.col("Pfam_hmm_id_clean").is_in(regulation_table["id"]) |
        pl.col("dbCAN_hmm_id").is_in(regulation_table["id"]) |
        pl.col("METABOLIC_hmm_id").is_in(regulation_table["id"]) |
        pl.col("PHROG_hmm_id").is_in(regulation_table["id"])
    )
    table = table.filter(condition)
    
    # Apply false-substring filtering
    table = filter_false_substrings(table, false_reg_substrings)
    
    # Drop the temporary 'top_hit_hmm_id_clean' column
    table = table.drop("top_hit_hmm_id_clean")

    # Remove duplicates, if any (this happens sometimes if the input table also had duplciates)
    table = table.unique()
    
    return table.sort(["Genome", "Contig", "Protein"])

def main():
    input_table  = snakemake.params.context_table
    hmm_ref = snakemake.params.hmm_ref
    metabolism_ref = snakemake.params.metabolism_table
    physiology_ref = snakemake.params.physiology_table
    regulation_ref = snakemake.params.regulation_table
    false_metab_substrings = snakemake.params.false_amgs
    false_phys_substrings = snakemake.params.false_apgs
    false_reg_substrings = snakemake.params.false_aregs
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
    
    metabolism_table = pl.read_csv(metabolism_ref, separator="\t")
    physiology_table = pl.read_csv(physiology_ref, separator="\t")
    regulation_table = pl.read_csv(regulation_ref, separator="\t")
    
    false_metab_substrings_desc = []
    false_phys_substrings_desc = []
    false_reg_substrings_desc = []

    for textfile, substrings_list in zip(
        [false_metab_substrings, false_phys_substrings, false_reg_substrings],
        [false_metab_substrings_desc, false_phys_substrings_desc, false_reg_substrings_desc]
    ):
        with open(textfile, 'r') as file:
            for line in file.readlines():
                if line.startswith("#") or line.strip() == "":
                    continue
                substrings_list.append(line.strip().replace('"', ''))
    
    annot_table = summarize_annot_table(table, hmm_descriptions)
    
    # Remove .X or .XX suffixes from top_hit_hmm_id for proper matching of Pfam hits
    annot_table = annot_table.with_columns(
        pl.col("top_hit_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("top_hit_hmm_id_clean"),
        pl.col("Pfam_hmm_id").str.replace(r'\.\d+$', '', literal=False).alias("Pfam_hmm_id_clean"),
    )
    
    metabolism_table_out = filter_metabolism_annots(annot_table, metabolism_table, false_metab_substrings_desc)
    physiology_table_out = filter_physiology_annots(annot_table, physiology_table, false_phys_substrings_desc)
    regulation_table_out = filter_regulation_annots(annot_table, regulation_table, false_reg_substrings_desc)
    
    drop_cols = ["window_avg_KEGG_VL-score_viral", "window_avg_Pfam_VL-score_viral", "window_avg_PHROG_VL-score_viral", "top_hit_hmm_id_clean", "Pfam_hmm_id_clean"]
    for col in drop_cols:
        if col in annot_table.columns:
            annot_table = annot_table.drop(col)
        if col in metabolism_table_out.columns:
            metabolism_table_out = metabolism_table_out.drop(col)
        if col in physiology_table_out.columns:
            physiology_table_out = physiology_table_out.drop(col)
        if col in regulation_table_out.columns:
            regulation_table_out = regulation_table_out.drop(col)

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