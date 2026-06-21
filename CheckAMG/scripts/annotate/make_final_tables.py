#!/usr/bin/env python3

import os
import sys
from pathlib import Path
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import gc

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Write the final summarized auxiliary gene tables", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")


def _ensure_function_col(df: pl.DataFrame) -> pl.DataFrame:
    if "Function" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias("Function"))
    return df

def _build_function_lookup(metabolism_df: pl.DataFrame, physiology_df: pl.DataFrame, regulatory_df: pl.DataFrame) -> pl.DataFrame:
    # Build per-protein Function lookup from category tables
    metab_fn = _ensure_function_col(metabolism_df).select(["Protein", "Function"]).rename({"Function": "Function_metabolic"})
    phys_fn = _ensure_function_col(physiology_df).select(["Protein", "Function"]).rename({"Function": "Function_physiological"})
    reg_fn = _ensure_function_col(regulatory_df).select(["Protein", "Function"]).rename({"Function": "Function_regulatory"})

    fn_lookup = (
        metab_fn
        .join(phys_fn, on="Protein", how='full', coalesce=True)
        .join(reg_fn, on="Protein", how='full', coalesce=True)
    )
    return fn_lookup

# Function to classify proteins based on their presence in metabolic, physiology, and regulatory tables
def classify_proteins(final_df, metabolism_df, physiology_df, regulatory_df):
    logger.debug(f"Columns in final_df before classification: {final_df.columns}")
    logger.debug(f"Columns in metabolism_df: {metabolism_df.columns}")
    logger.debug(f"Columns in physiology_df: {physiology_df.columns}")
    logger.debug(f"Columns in regulatory_df: {regulatory_df.columns}")

    # Ensure column names are consistent
    required_col = "Protein"
    for df_name, df in zip(["metabolism_df", "physiology_df", "regulatory_df"], [metabolism_df, physiology_df, regulatory_df]):
        if required_col not in df.columns:
            logger.error(f"Column '{required_col}' not found in {df_name}")
            raise ValueError(f"Column '{required_col}' not found in {df_name}")

    # Extract unique protein lists for fast lookup
    metabolic_proteins = set(metabolism_df["Protein"].to_list())
    physiology_proteins = set(physiology_df["Protein"].to_list())
    regulatory_proteins = set(regulatory_df["Protein"].to_list())

    # Assign classifications
    def classify(protein):
        if protein in metabolic_proteins:
            return "metabolic"
        elif protein in physiology_proteins:
            return "physiological"
        elif protein in regulatory_proteins:
            return "regulatory"
        else:
            return "unclassified"

    final_df = _ensure_function_col(final_df)

    # Attach category-derived Function values
    fn_lookup = _build_function_lookup(metabolism_df, physiology_df, regulatory_df)

    final_df = (
        final_df
        .with_columns(
            pl.col("Protein").map_elements(classify, return_dtype=pl.Utf8).alias("classification")
        )
        .join(fn_lookup, on="Protein", how="left")
        # Populate Function based on inferred classification, fallback to existing Function, then top_hit_description
        .with_columns(
            pl.when(pl.col("classification") == "metabolic").then(pl.col("Function_metabolic"))
            .when(pl.col("classification") == "physiological").then(pl.col("Function_physiological"))
            .when(pl.col("classification") == "regulatory").then(pl.col("Function_regulatory"))
            .otherwise(pl.col("Function"))
            .alias("Function_tmp")
        )
        .with_columns(
            pl.coalesce([pl.col("Function_tmp"), pl.col("Function"), pl.col("top_hit_description")]).alias("Function")
        )
        .drop(["Function_tmp", "Function_metabolic", "Function_physiological", "Function_regulatory"])
        .select(
            [
                "Protein",
                "Contig",
                "Genome",
                "classification",
                "Viral_Origin_Confidence",
                "step5_in_merged_region",
                "Function",
                # "Circular_Contig", "Viral_Flanking_Genes_Left_Dist", "Viral_Flanking_Genes_Right_Dist", "MGE_Flanking_Genes_Left_Dist", "MGE_Flanking_Genes_Right_Dist" ## Debugging
                "KEGG_hmm_id",
                "KEGG_Description",
                "FOAM_hmm_id",
                "FOAM_Description",
                "Pfam_hmm_id",
                "Pfam_Description",
                "dbCAN_hmm_id",
                "dbCAN_Description",
                "METABOLIC_hmm_id",
                "METABOLIC_Description",
                "CAMPER_hmm_id",
                "CAMPER_Description",
                "PHROG_hmm_id",
                "PHROG_Description",
                # "VOG_hmm_id",
                # "VOG_Description",
                "top_hit_hmm_id",
                "top_hit_description",
                "top_hit_db",
            ]
        )
        .rename(
            {
                "classification": "Protein Classification",
                "Viral_Origin_Confidence": "Protein Viral Origin Confidence",
                "step5_in_merged_region": "Protein in Strict Viral Region",
                "KEGG_hmm_id": "KEGG KO",
                "KEGG_Description": "KEGG KO Name",
                "FOAM_hmm_id": "FOAM ID",
                "FOAM_Description": "FOAM Annotation",
                "Pfam_hmm_id": "Pfam Accession",
                "Pfam_Description": "Pfam Name",
                "dbCAN_hmm_id": "CAZy Family",
                "dbCAN_Description": "CAZy Activities",
                "METABOLIC_hmm_id": "METABOLIC db ID",
                "METABOLIC_Description": "METABOLIC Annotation",
                "CAMPER_hmm_id": "CAMPER ID",
                "CAMPER_Description": "CAMPER Annotation",
                "PHROG_hmm_id": "PHROG Number",
                "PHROG_Description": "PHROG Annotation",
                # "VOG_hmm_id": "VOG Name",
                # "VOG_Description": "VOG Annotation",
                "top_hit_hmm_id": "Best Scoring HMM",
                "top_hit_description": "Best Scoring HMM Annotation",
                "top_hit_db": "Best Scoring HMM Origin",
            }
        )
    )

    # Free memory
    del metabolism_df, physiology_df, regulatory_df, fn_lookup
    gc.collect()

    return final_df

# Function to join the auxiliary status, hmm dataframe, all genes dataframe, and classification for CheckAMG annotate mode
def merge_dataframes(all_genes_df, metabolism_df, physiology_df, regulatory_df):
    logger.debug(f"Columns in all_genes_df: {all_genes_df.columns}")
    logger.debug(f"all_genes_df shape: {all_genes_df.shape}")
    return classify_proteins(all_genes_df, metabolism_df, physiology_df, regulatory_df)

def main():
    save_to_parquet = snakemake.params.save_to_parquet
    all_genes_path = as_parquet_path_if_enabled(snakemake.params.all_genes_annotated, save_to_parquet)
    gene_index_path = as_parquet_path_if_enabled(snakemake.params.gene_index, save_to_parquet)
    metab_cat_ref_path = snakemake.params.metabolic_categories_ref
    phys_cat_ref_path = snakemake.params.physiology_categories_ref
    reg_cat_ref_path = snakemake.params.regulation_categories_ref
    metabolism_path = as_parquet_path_if_enabled(snakemake.params.metabolism_table, save_to_parquet)
    metab_cat_out_path = as_parquet_path_if_enabled(snakemake.params.metabolism_categories_out, save_to_parquet)
    physiology_path = as_parquet_path_if_enabled(snakemake.params.physiology_table, save_to_parquet)
    phys_cat_out_path = as_parquet_path_if_enabled(snakemake.params.physiology_categories_out, save_to_parquet)
    regulatory_path = as_parquet_path_if_enabled(snakemake.params.regulation_table, save_to_parquet)
    reg_cat_out_path = as_parquet_path_if_enabled(snakemake.params.regulation_categories_out, save_to_parquet)
    final_table_path = as_parquet_path_if_enabled(snakemake.params.final_table, save_to_parquet)
    mem_limit = snakemake.resources.mem
    threads = snakemake.threads
    set_memory_limit(mem_limit)

    logger.info(f"Generating the final table with proteins, annotations, and classifications...")

    if save_to_parquet:
        logger.debug("Reading 'all genes', 'gene index', and metabolism/physiology/regulatory gene tables as parquet...")
        all_genes_df = pl.read_parquet(all_genes_path)
        gene_index_df = pl.read_parquet(gene_index_path)
        metabolism_df = pl.read_parquet(metabolism_path)
        physiology_df = pl.read_parquet(physiology_path)
        regulatory_df = pl.read_parquet(regulatory_path)
    else:
        logger.debug("Reading 'all genes', 'gene index', and metabolism/physiology/regulatory gene tables as TSV...")
        all_genes_df = pl.read_csv(all_genes_path, separator='\t')
        gene_index_df = pl.read_csv(gene_index_path, separator='\t')
        metabolism_df = pl.read_csv(metabolism_path, separator='\t')
        physiology_df = pl.read_csv(physiology_path, separator='\t')
        regulatory_df = pl.read_csv(regulatory_path, separator='\t')

    logger.debug("Loading reference AVG category tables...")
    metab_cat_ref_df = pl.read_csv(metab_cat_ref_path, separator='\t')
    phys_cat_ref_df = pl.read_csv(phys_cat_ref_path, separator='\t')
    reg_cat_ref_df = pl.read_csv(reg_cat_ref_path, separator='\t')

    gene_index_df = gene_index_df.select(["protein"] + [
        col for col in gene_index_df.columns if col.endswith("_protein_cluster") or col.endswith("_genome_cluster")
    ]).rename({"protein": "Protein"})

    logger.debug("Merging 'all genes', metabolism, physiology, and regulatory genes dataframes...")
    final_df = merge_dataframes(all_genes_df, metabolism_df, physiology_df, regulatory_df)
    final_df = final_df.join(gene_index_df, on="Protein", how="left")

    # Sort and write the final table output
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    final_df = final_df.with_columns(
        pl.col("Protein Viral Origin Confidence").replace(confidence_order).cast(pl.Int32).alias("Protein Viral Origin Confidence_sort")
    )
    classification_order = {"metabolic": 0, "physiological": 1, "regulatory": 2, "unclassified": 3}
    final_df = final_df.with_columns(
        pl.col("Protein Classification").replace(classification_order).cast(pl.Int32).alias("Protein Classification_sort")
    ).sort(["Protein Classification_sort", "Protein Viral Origin Confidence_sort", "Protein"]).drop(["Protein Classification_sort", "Protein Viral Origin Confidence_sort"])

    if save_to_parquet:
        final_df.write_parquet(final_table_path)
    else:
        final_df.write_csv(final_table_path, separator='\t')
    logger.info(f"Final table written to {final_table_path}")

    logger.debug(f"Joining categories to metabolic proteins...")
    metab_cat_df = (
        final_df
        .filter(pl.col("Protein Classification") == "metabolic")
        .select([
            "Protein", "Contig", "Genome",
            "Protein Viral Origin Confidence", "Protein in Strict Viral Region",
            "Function"
        ])
        .join(
            metab_cat_ref_df
            .select(["name", "category_L1", "category_L2", "category_L3"])
            .unique(),
            left_on="Function", right_on="name", how="left"
            )
    )
    logger.info(f"{metab_cat_df.get_column('Protein').n_unique():,} metabolic proteins span:")
    logger.info(f"\t{metab_cat_df.get_column('category_L1').n_unique():2,} level 1 AMG categories")
    logger.info(f"\t{metab_cat_df.get_column('category_L2').n_unique():2,} level 2 AMG cateogires")
    logger.info(f"\t{metab_cat_df.get_column('category_L3').n_unique():2,} level 3 AMG cateogires")
    logger.debug(f"Writing metabolic category table to {metab_cat_out_path}...")
    if save_to_parquet:
        metab_cat_df.write_parquet(metab_cat_out_path)
    else:
        metab_cat_df.write_csv(metab_cat_out_path, separator='\t')

    logger.debug(f"Joining categories to physiological proteins...")
    phys_cat_df = (
        final_df
        .filter(pl.col("Protein Classification") == "physiological")
        .select([
            "Protein", "Contig", "Genome",
            "Protein Viral Origin Confidence", "Protein in Strict Viral Region",
            "Function"
        ])
        .join(
            phys_cat_ref_df
            .select(["name", "category_L1", "category_L2", "category_L3"])
            .unique(),
            left_on="Function", right_on="name", how="left"
            )
    )
    logger.info(f"{phys_cat_df.get_column('Protein').n_unique():,} physiological proteins span:")
    logger.info(f"\t{phys_cat_df.get_column('category_L1').n_unique():2,} level 1 APG categories")
    logger.info(f"\t{phys_cat_df.get_column('category_L2').n_unique():2,} level 2 APG cateogires")
    logger.info(f"\t{phys_cat_df.get_column('category_L3').n_unique():2,} level 3 APG cateogires")
    logger.debug(f"Writing physiological category table to {phys_cat_out_path}...")
    if save_to_parquet:
        phys_cat_df.write_parquet(phys_cat_out_path)
    else:
        phys_cat_df.write_csv(phys_cat_out_path, separator='\t')

    logger.debug(f"Joining categories to regulatory proteins...")
    reg_cat_df = (
        final_df
        .filter(pl.col("Protein Classification") == "regulatory")
        .select([
            "Protein", "Contig", "Genome",
            "Protein Viral Origin Confidence", "Protein in Strict Viral Region",
            "Function"
        ])
        .join(
            reg_cat_ref_df
            .select(["name", "category_L1", "category_L2", "category_L3"])
            .unique(),
            left_on="Function", right_on="name", how="left"
            )
    )
    logger.info(f"{reg_cat_df.get_column('Protein').n_unique():,} regulatory proteins span:")
    logger.info(f"\t{reg_cat_df.get_column('category_L1').n_unique():2,} level 1 AReG categories")
    logger.info(f"\t{reg_cat_df.get_column('category_L2').n_unique():2,} level 2 AReG cateogires")
    logger.info(f"\t{reg_cat_df.get_column('category_L3').n_unique():2,} level 3 AReG cateogires")
    logger.debug(f"Writing regulatory category table to {reg_cat_out_path}...")

    # # Intentional polars error to test logger
    # reg_cat_df = reg_cat_df.with_columns(pl.col("category_L1").cast(pl.Int32)) # DEBUGGING

    if save_to_parquet:
        reg_cat_df.write_parquet(reg_cat_out_path)
    else:
        reg_cat_df.write_csv(reg_cat_out_path, separator='\t')

    # Log classification summary
    logger.debug(f"Columns in final_df after classification: {final_df.columns}")
    logger.debug(f"Protein Classification value counts:\n{final_df['Protein Classification'].value_counts()}")
    logger.debug(f"Protein Viral Origin Confidence value counts:\n{final_df['Protein Viral Origin Confidence'].value_counts()}")

if __name__ == "__main__":
    main()
