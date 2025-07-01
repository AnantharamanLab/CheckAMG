#!/usr/bin/env python3

import os
import sys
import resource
import platform
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

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

print("========================================================================\n     Step 7/11: Merge functional annotations with protein metadata      \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n     Step 7/11: Merge functional annotations with protein metadata      \n========================================================================\n")

def prefix_columns(dataframe, prefix):
    # Select 'sequence' only once, and prefix other columns as needed
    cols_to_select = [pl.col("sequence")]
    cols_to_select.extend(
        pl.col(col).alias(f"{prefix}_{col}") for col in dataframe.columns if col not in ["sequence", "db", "score"]
    )
    return dataframe.select(cols_to_select)

def widen_hmm_results(hmm_results_df):
    # Ensure consistent ordering
    hmm_results_df = hmm_results_df.sort(["sequence", "db"])

    dbs = hmm_results_df.select("db").unique().to_series().to_list()

    hmm_results_df = hmm_results_df.with_columns([
        pl.when(pl.col("db") == db).then(pl.col("score")).otherwise(None).alias(f"{db}_score")
        for db in dbs
    ] + [
        pl.when(pl.col("db") == db).then(pl.col("hmm_id")).otherwise(None).alias(f"{db}_hmm_id")
        for db in dbs
    ] + [
        pl.when(pl.col("db") == db).then(pl.col("hmm_coverage")).otherwise(None).alias(f"{db}_coverage")
        for db in dbs
    ])

    # Drop original columns
    hmm_results_df = hmm_results_df.drop(["score", "hmm_id", "hmm_coverage", "db"])

    # Aggregate by sequence
    hmm_results_df_wide = hmm_results_df.group_by("sequence", maintain_order=True).agg([
        pl.max(col).alias(col) for col in hmm_results_df.columns if col != "sequence"
    ])

    return hmm_results_df_wide

def assign_db(db_path):
    if "KEGG" in str(db_path) or "kegg" in str(db_path) or "kofam" in str(db_path):
        return "KEGG"
    elif "FOAM" in str(db_path) or "foam" in str(db_path):
        return "FOAM"
    elif "Pfam" in str(db_path) or "pfam" in str(db_path):
        return "Pfam"
    elif "dbcan" in str(db_path) or "dbCAN" in str(db_path) or "dbCan" in str(db_path):
        return "dbCAN"
    elif "METABOLIC_custom" in str(db_path) or "metabolic_custom" in str(db_path):
        return "METABOLIC"
    elif "VOG" in str(db_path) or "vog" in str(db_path):
        return "VOG"
    elif "eggNOG" in str(db_path) or "eggnog" in str(db_path):
        return "eggNOG"
    elif "PHROG" in str(db_path) or "phrog" in str(db_path):
        return "PHROG"
    elif "user_custom" in str(db_path):
        return "user_custom"
    else:
        return None
    
def main():
    tsv = snakemake.params.gene_index
    vscores = snakemake.params.vscores
    all_hmm_results = snakemake.params.all_hmm_results
    db_dir = snakemake.params.db_dir
    output = snakemake.params.gene_index_annotated
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    logger.info("Processing of V-scores/VL-scores and HMM results starting...")
    
    # Load input dataframes
    logger.debug(f"Loading input dataframes {tsv}, {vscores}, and {all_hmm_results}")
    tsv_df = pl.read_csv(tsv, separator="\t").unique()
    vscores_df = pl.read_csv(vscores, separator="\t").unique()
    hmm_df = pl.read_csv(all_hmm_results, separator="\t").unique()
    
    hmm_df_wide = widen_hmm_results(hmm_df)

    # Merge HMM results with the main dataframe
    merged_df = tsv_df.join(hmm_df_wide, left_on='protein', right_on='sequence', how='left')
    
    # Split the DataFrame by 'db' value
    df_pfam = vscores_df.filter(pl.col("db") == "Pfam")
    df_kegg = vscores_df.filter(pl.col("db") == "KEGG")
    df_phrog = vscores_df.filter(pl.col("db") == "PHROG")

    # Add prefixes
    pfam_prefixed = prefix_columns(df_pfam, "Pfam")
    kegg_prefixed = prefix_columns(df_kegg, "KEGG")
    phrog_prefixed = prefix_columns(df_phrog, "PHROG")

    # Join on 'sequence'
    wide_df = pfam_prefixed.join(kegg_prefixed, on="sequence", how="full")
    cols_to_remove = [col for col in wide_df.columns if col.endswith('_right') and not col.endswith('_score_right')]
    wide_df = wide_df.drop([col for col in cols_to_remove if col in wide_df.columns])
    wide_df = wide_df.join(phrog_prefixed, on="sequence", how="full")
    
    # Merge V-scores and L-scores with input TSV
    logger.debug(f"Merging V-scores/VL-scores and annotations with input {tsv} and writing results to {output}")
    merged_df = merged_df.join(wide_df, left_on='protein', right_on='sequence', how='left')
    
    # Ensure all DBs are represented in columns, even if there were no hits
    for db_path in os.listdir(db_dir):
        db_name = assign_db(db_path)
        col_found = any(db_name in col for col in merged_df.columns)
        if not col_found:
            logger.debug(f"No {db_name} found in the merged DataFrame; adding columns")
            merged_df = merged_df.with_columns(
                pl.lit(None).alias(f"{db_name}_score"),
                pl.lit(None).alias(f"{db_name}_hmm_id")
            )
    
    # Remove unnecessary columns and those ending with '_right'
    cols_to_remove = [
        'ID', 'rank', 'protein_cluster_rep', 'genome_cluster'
    ]
    # Add all columns ending with '_right' (except those ending with '_score_right')
    cols_to_remove += [col for col in merged_df.columns if col.endswith('_right') and not col.endswith('_score_right')]
    
    # Drop specified columns
    merged_df = merged_df.drop([col for col in cols_to_remove if col in merged_df.columns])
    
    # Rename and handle '_score_right' columns
    for col in merged_df.columns:
        if col.endswith('_score_right'):
            original_col = col.replace('_score_right', '_score')
            # Check if the original column already exists
            if original_col in merged_df.columns:
                logger.debug(f"Dropping existing column '{original_col}' before renaming '{col}' to avoid duplicates.")
                # Drop the existing column to avoid duplicates
                merged_df = merged_df.drop(original_col)
            # Rename '_score_right' to '_score'
            merged_df = merged_df.rename({col: original_col})

    # Remove duplicates and sort the final DataFrame
    merged_df = merged_df.unique()
    merged_df = merged_df.sort(["genome", "contig", "gene_number", "protein"])
    
    # Write the final merged DataFrame to the output file
    merged_df.write_csv(output, separator="\t")
    logger.info("Processing of V-scores/VL-scores and HMM results completed.")

if __name__ == "__main__":
    main()
