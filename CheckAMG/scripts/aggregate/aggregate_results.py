#!/usr/bin/env python3

import os
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Aggregate results from the annotate and de-novo modules", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

def read_file(file_path):
    if file_path.endswith('.parquet'):
        return pl.read_parquet(file_path)
    elif file_path.endswith('.tsv'):
        return pl.read_csv(file_path, separator='\t')
    else:
        raise ValueError(f"Unsupported file format: {file_path}")

def main():
    save_to_parquet = snakemake.params.save_to_parquet
    annotate_results_path = snakemake.params.annotate_final_results
    annotate_context_path = snakemake.params.annotate_genomic_context
    annotate_metab_cats_path = snakemake.params.annotate_metab_cats
    annotate_phys_cats_path = snakemake.params.annotate_phys_cats
    annotate_reg_cats_path = snakemake.params.annotate_reg_cats
    denovo_predictions_path = snakemake.params.denovo_predictions
    output_main_results = as_parquet_path_if_enabled(snakemake.params.output_main_results, save_to_parquet)
    output_detailed_results = as_parquet_path_if_enabled(snakemake.params.output_detailed_results, save_to_parquet)
    output_categories = as_parquet_path_if_enabled(snakemake.params.output_categories, save_to_parquet)

    logger.info(f"Aggregation of CheckAMG annotate and de-novo results starting...")

    # Read the CheckAMG annotate genomic context resutls with LGBM predictions
    logger.debug(f"Reading CheckAMG annotate genomic context table with LGBM predictions {annotate_context_path}...")
    required_cols_from_context = [
        "genome", "contig", "protein",
        "LGBM_viral_prob", "Viral_Origin_Confidence"
    ]
    annotate_context = (
        read_file(annotate_context_path)
        .select(required_cols_from_context)
        .rename({
            "genome": "Genome",
            "contig": "Contig",
            "protein": "Protein",
        })
    )

    # Read the CheckAMG annotate results with functional annotations
    logger.debug(f"Reading CheckAMG annotate results from {annotate_results_path} and joining with the LGBM predictions...")
    required_cols_from_annotate = [
        "Protein", "Contig", "Genome",
        "Protein Classification", "Function"
    ]
    annotate_results = (
        read_file(annotate_results_path)
        .select(required_cols_from_annotate)
        .join(
            annotate_context,
            on=["Protein", "Contig", "Genome"],
            how="full", coalesce=True
        )
    )

    # Read the de-novo predictions
    logger.debug(f"Reading CheckAMG de-novo predictions from {denovo_predictions_path}...")
    required_cols_from_denovo = [
        "Protein", "Contig", "Genome",
        "AVG-like prob", "AVG-like confidence",
        "Viral prob", "Viral confidence",
        "Final AVG prob", "Final AVG confidence"
    ]
    denovo_results = read_file(denovo_predictions_path).select(required_cols_from_denovo)

    # Join the results
    # coalesce=True merges the join keys into single
    # Protein/Contig/Genome columns so de-novo-only proteins keep their keys
    # (otherwise their keys would live in the dropped *_right columns).
    logger.debug("Joining the annotate and de-novo results...")
    merged_pre = annotate_results.join(
        denovo_results, on=["Protein", "Contig", "Genome"], how="full", coalesce=True
    )

    # Build the detailed output
    logger.debug("Building the aggregated results table (detailed)...")
    rename_dict = {
        "Protein Classification": "Classification (annotate)",
        "Function": "Function (annotate)",
        "LGBM_viral_prob": "Viral Probability (annotate)",
        "Viral_Origin_Confidence": "Viral Confidence Level (annotate)",
        "Viral prob": "Viral Probability (de-novo)",
        "Viral confidence": "Viral Confidence Level (de-novo)",
        "AVG-like prob": "AVG-like Probability (de-novo)",
        "AVG-like confidence": "AVG-like Confidence Level (de-novo)",
        "Final AVG prob": "Final AVG Probability (de-novo)",
        "Final AVG confidence": "Final AVG Confidence Level (de-novo)",
    }
    merged_detailed = (
        merged_pre
        .rename(rename_dict)
        .select(
            [
                "Protein", "Contig", "Genome",
                "Viral Probability (annotate)","Viral Confidence Level (annotate)",
                "Viral Probability (de-novo)", "Viral Confidence Level (de-novo)",
                "Classification (annotate)", "Function (annotate)",
                "AVG-like Probability (de-novo)", "AVG-like Confidence Level (de-novo)",
                "Final AVG Probability (de-novo)", "Final AVG Confidence Level (de-novo)"
            ]
        )
    )

    # Save the detailed results
    logger.debug(f"Saving the (detailed) aggregated results to {output_detailed_results}")
    if save_to_parquet:
        merged_detailed.write_parquet(output_detailed_results)
    else:
        merged_detailed.write_csv(output_detailed_results, separator='\t')


    # Build the main output
    logger.debug("Building the aggregated results table (simplified)...")
    main_rename_dict = {
        "Protein Classification": "Classification",
        "Viral prob": "Viral Probability",
        "Viral confidence": "Viral Confidence Level",
        "AVG-like prob": "AVG-like Probability",
        "AVG-like confidence": "AVG-like Confidence Level",
        "Final AVG prob": "Final AVG Probability",
        "Final AVG confidence": "Final AVG Confidence Level",
    }
    merged_main = (
        merged_pre
        .rename(main_rename_dict)
        .select([
            "Protein", "Contig", "Genome",
            "Classification", "Function",
            "Final AVG Probability", "Final AVG Confidence Level"
        ])
    )

    # Save the main results
    logger.info(f"Saving the aggregated results to {output_main_results}")
    if save_to_parquet:
        merged_main.write_parquet(output_main_results)
    else:
        merged_main.write_csv(output_main_results, separator='\t')
    

    # Build the table with functional-categories
    required_cols_from_cats = [
        "Protein", "Contig", "Genome",
        "category_L1", "category_L2", "category_L3"
    ]
    logger.debug(f"Reading CheckAMG metabolic categories from {annotate_metab_cats_path}")
    metab_cats = read_file(annotate_metab_cats_path).select(required_cols_from_cats)

    logger.debug(f"Reading CheckAMG physiological categories from {annotate_phys_cats_path}")
    phys_cats = read_file(annotate_phys_cats_path).select(required_cols_from_cats)

    logger.debug(f"Reading CheckAMG regulatory categories from {annotate_reg_cats_path}")
    reg_cats = read_file(annotate_reg_cats_path).select(required_cols_from_cats)

    # Stack the three category sources, then LEFT JOIN them onto `merged`
    # so that every protein from annotate and/or de-novo is represented,
    # even if it has no mapped category.

    all_cats = pl.concat([metab_cats, phys_cats, reg_cats], how="vertical")
    merged_cats = merged_main.join(all_cats, on=["Protein", "Contig", "Genome"], how="left")

    # Assign placeholder category levels where none were mapped:
    #   - Function is null (e.g. de-novo-only proteins): "Unknown"
    #   - Function present but Classification is "unclassified" with no mapped
    #     category: "Other"
    # Proteins that already have a mapped category keep their values.

    category_levels = ["category_L1", "category_L2", "category_L3"]
    no_mapped_category = pl.col("category_L1").is_null()
    merged_cats = merged_cats.with_columns(
        [
            pl.when(pl.col("Function").is_null())
            .then(pl.lit("Unknown"))
            .when((pl.col("Classification") == "unclassified") & no_mapped_category)
            .then(pl.lit("Other"))
            .otherwise(pl.col(level))
            .alias(level)
            for level in category_levels
        ]
    )

    merged_cats = merged_cats.select(
        [
            "Protein", "Contig", "Genome",
            "Classification", "category_L1", "category_L2", "category_L3", "Function",
            "Final AVG Probability", "Final AVG Confidence Level"
        ]
    )

    # Save the categories
    logger.debug(f"Saving the aggregated results with functional categories to {output_categories}")
    if save_to_parquet:
        merged_cats.write_parquet(output_categories)
    else:
        merged_cats.write_csv(output_categories, separator='\t')


    logger.info(f"Aggregation of CheckAMG annotate and de-novo results completed.")

if __name__ == "__main__":
    main()

