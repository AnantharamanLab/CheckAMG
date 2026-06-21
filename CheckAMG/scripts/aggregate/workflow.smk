import os

rule all:
    input:
        os.path.join(config["paths"]["output_dir"], "snakemake", "aggregate.done")

rule aggregate:
    output:
        touch(os.path.join(config["paths"]["output_dir"], "snakemake", "aggregate.done"))
    params:
        annotate_final_results = config["paths"]["annotate_final_results"],
        annotate_genomic_context = config["paths"]["annotate_genomic_context"],
        annotate_metab_cats = config["paths"]["annotate_metab_cats"],
        annotate_phys_cats = config["paths"]["annotate_phys_cats"],
        annotate_reg_cats = config["paths"]["annotate_reg_cats"],
        denovo_predictions = config["paths"]["denovo_predictions"],
        output_main_results = os.path.join(config["paths"]["output_dir"], "aggregated_results.tsv"),
        output_categories = os.path.join(config["paths"]["output_dir"], "aggregated_results_categories.tsv"),
        output_detailed_results = os.path.join(config["paths"]["output_dir"], "aggregated_results_detailed.tsv"),
        save_to_parquet = bool(config["save_to_parquet"]),
        debug = bool(config["debug"]),
        log = config["log"],
        step = (1,1)
    threads:
        config["threads"]
    resources:
        mem = config["mem_limit"]
    message:
        "Aggregating results from CheckAMG annotate and CheckAMG de-novo."
    script:
        os.path.join(config["paths"]["steps_dir"], "aggregate_results.py")