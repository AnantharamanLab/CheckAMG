#!/usr/bin/env python3

import os
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled
from CheckAMG.scripts.annotate import summaries

AVG_TYPE = summaries.AVG_TYPE

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
    annotate_avg_path = snakemake.params.annotate_avg_predictions
    annotate_supported_conf = list(snakemake.params.annotate_supported_conf)
    denovo_supported_conf = list(snakemake.params.denovo_supported_conf)
    output_main_results = as_parquet_path_if_enabled(snakemake.params.output_main_results, save_to_parquet)
    output_detailed_results = as_parquet_path_if_enabled(snakemake.params.output_detailed_results, save_to_parquet)
    output_avg_predictions = as_parquet_path_if_enabled(snakemake.params.output_avg_predictions, save_to_parquet)
    output_genome_summary = as_parquet_path_if_enabled(snakemake.params.output_genome_summary, save_to_parquet)
    os.makedirs(os.path.dirname(output_main_results), exist_ok=True)

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
        "Protein Classification", "Function", "Protein in Strict Viral Region"
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
    # The confidence levels that count as support for each module are set in the Snakemake rule
    annotate_support = (
        pl.col("Classification").is_in(list(AVG_TYPE)).fill_null(False)
        & pl.col("Viral_Origin_Confidence").is_in(annotate_supported_conf).fill_null(False)
    )
    denovo_support = pl.col("Final AVG Confidence Level").is_in(denovo_supported_conf).fill_null(False)
    merged_main = (
        merged_pre
        .rename(main_rename_dict)
        .with_columns(
            pl.when(annotate_support & denovo_support).then(pl.lit("annotate and de-novo"))
            .when(annotate_support).then(pl.lit("annotate only"))
            .when(denovo_support).then(pl.lit("de-novo only"))
            .otherwise(pl.lit("not supported"))
            .alias("Evidence Class"),
            pl.col("Classification").replace_strict(AVG_TYPE, default=None, return_dtype=pl.Utf8).alias("AVG Type"),
            pl.col("Viral_Origin_Confidence").alias("Viral Origin Confidence (annotate)"),
        )
        .select([
            "Protein", "Contig", "Genome",
            "Evidence Class",
            "Classification", "AVG Type", "Function",
            "Viral Origin Confidence (annotate)",
            "Final AVG Probability", "Final AVG Confidence Level",
            "Protein in Strict Viral Region",
        ])
    )

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

    all_cats = summaries.collapse_categories(pl.concat([metab_cats, phys_cats, reg_cats], how="vertical"))
    merged_cats = merged_main.join(all_cats, on=["Protein", "Contig", "Genome"], how="left")

    # Assign placeholder category levels where none were mapped:
    #   - Function is null (e.g. de-novo-only proteins): "Unknown"
    #   - Function present but Classification is "unclassified" with no mapped
    #     category: "Other"
    # Proteins that already have a mapped category keep their values.

    category_levels = [summaries.CAT_TOP, summaries.CAT_PATHS]
    no_mapped_category = pl.col(summaries.CAT_TOP).is_null()
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

    all_proteins = merged_cats.select(
        [
            "Protein", "Contig", "Genome",
            "Evidence Class", "AVG Type", "Function",
            "Final AVG Confidence Level", "Final AVG Probability",
            "Viral Origin Confidence (annotate)",
            summaries.CAT_TOP, summaries.CAT_PATHS,
            "Protein in Strict Viral Region",
        ]
    )
    logger.debug(f"Saving the aggregated results for all proteins to {output_main_results}")
    summaries.write_table(all_proteins, output_main_results)

    annotate_avg = summaries.read_table(annotate_avg_path) if annotate_avg_path else None
    avg_pred = summaries.build_aggregate_avg_predictions(all_proteins, annotate_avg)
    summaries.write_table(avg_pred, output_avg_predictions)
    summaries.write_table(summaries.build_aggregate_genome_summary(all_proteins, avg_pred), output_genome_summary)
    logger.info(f"{avg_pred.height:,} AVGs supported by annotate, de-novo, or both written to {output_avg_predictions}")

    run_rows = [
        ("run", "CheckAMG version", summaries.checkamg_version()),
        ("run", "Module", "aggregate"),
        ("run", "Annotate all-protein table", annotate_results_path),
        ("run", "De-novo predictions", denovo_predictions_path),
        ("definitions", "Annotate support", "classified as AMG, APG, or AReG with viral-origin confidence " + " or ".join(annotate_supported_conf)),
        ("definitions", "De-novo support", "Final AVG confidence " + ", ".join(denovo_supported_conf)),
        ("inputs", "Proteins", all_proteins.height),
    ]
    ev = all_proteins.group_by(["Evidence Class", "AVG Type"]).len()
    for evidence in summaries.EVIDENCE_CLASSES:
        sub = ev.filter(pl.col("Evidence Class") == evidence)
        run_rows.append(("predictions", f"AVGs supported by {evidence}", int(sub["len"].sum())))
        for t in ["AMG", "APG", "AReG"]:
            n = int(sub.filter(pl.col("AVG Type") == t)["len"].sum())
            run_rows.append(("predictions", f"AVGs supported by {evidence}: {t}", n))
        if evidence == "de-novo only":
            n = int(sub.filter(pl.col("AVG Type").is_null())["len"].sum())
            run_rows.append(("predictions", f"AVGs supported by {evidence}: no annotate AVG type", n))
    for c in summaries.FINAL_CONFIDENCES:
        run_rows.append(("predictions", f"Supported AVGs with {c} final AVG confidence", avg_pred.filter(pl.col("Final AVG Confidence Level") == c).height))
    pl.DataFrame([(a, b, str(c)) for a, b, c in run_rows], schema=["Section", "Metric", "Value"], orient="row").write_csv(
        snakemake.params.output_run_summary, separator="\t"
    )
    summaries.write_aggregate_readme(snakemake.params.output_readme, all_proteins, avg_pred, annotate_supported_conf, denovo_supported_conf)
    logger.info(f"See {snakemake.params.output_readme} for a guide to the output files.")

    logger.info(f"Aggregation of CheckAMG annotate and de-novo results completed.")

if __name__ == "__main__":
    main()

