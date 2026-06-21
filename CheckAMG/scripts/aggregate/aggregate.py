#!/usr/bin/env python3

import yaml
import os
from pathlib import Path
from CheckAMG.scripts.common.runner_logging import setup_runner_logger, redirect_streams_to_log
from CheckAMG.scripts.common.set_up_snakemake import get_paths

def check_inputs(args, logger):
    """Validate that all required inputs exist and resolve their full paths.

    Returns a tuple ``(valid, resolved_paths)`` where ``valid`` is a bool
    indicating whether every required input was found, and ``resolved_paths``
    maps each config path key to the absolute path of the file that was found
    (whether the .tsv or .parquet variant). Keys for missing files are omitted.
    """
    valid = True
    resolved_paths = {}

    # Check annotate directory
    if not Path(args.annotate_dir).exists() or not Path(args.annotate_dir).is_dir():
        logger.error(f"CheckAMG annotate directory does not exist or is not a directory: {args.annotate_dir}")
        return False, resolved_paths

    results_dir = Path(args.annotate_dir).joinpath("results")

    if not results_dir.exists() or not results_dir.is_dir():
        logger.error(f"Missing results directory in annotate directory: {results_dir}")
        valid = False

    # Map each config path key to the candidate file names (tsv preferred, parquet fallback)
    required_annotate_files = {
        "annotate_final_results": ["final_results.tsv", "final_results.parquet"],
        "annotate_genomic_context": ["genes_genomic_context.tsv", "genes_genomic_context.parquet"],
        "annotate_metab_cats": ["metabolic_gene_categories.tsv", "metabolic_gene_categories.parquet"],
        "annotate_phys_cats": ["physiology_gene_categories.tsv", "physiology_gene_categories.parquet"],
        "annotate_reg_cats": ["regulation_gene_categories.tsv", "regulation_gene_categories.parquet"],
    }

    for config_key, file_options in required_annotate_files.items():
        found_path = None
        for filename in file_options:
            filepath = results_dir.joinpath(filename)
            if filepath.exists() and filepath.is_file():
                found_path = os.path.abspath(filepath)
                break
        if found_path is None:
            logger.error(
                f"Missing required annotate results file. Expected one of: "
                f"{', '.join(file_options)}"
            )
            valid = False
        else:
            resolved_paths[config_key] = found_path

    # Check de-novo directory
    if not Path(args.denovo_dir).exists() or not Path(args.denovo_dir).is_dir():
        logger.error(f"CheckMAG de-novo directory does not exist or is not a directory: {args.denovo_dir}")
        return False, resolved_paths

    denovo_predictions = None
    for filename in ["predictions.tsv", "predictions.parquet"]:
        filepath = Path(args.denovo_dir).joinpath(filename)
        if filepath.exists() and filepath.is_file():
            denovo_predictions = os.path.abspath(filepath)
            break
    if denovo_predictions is None:
        logger.error(
            "Missing required de-novo predictions file. "
            "Expected one of: predictions.tsv, predictions.parquet"
        )
        valid = False
    else:
        resolved_paths["denovo_predictions"] = denovo_predictions

    return valid, resolved_paths

def generate_config(args):
    """Generate a YAML config file based on provided arguments."""

    # Define the log file path and set up the logger
    log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_aggregate.log')
    # Structured logger for explicit log.info() / log/debug() / log.warning() / log.error() calls
    logger = setup_runner_logger(log_file_path, debug=args.debug)
    # Mirror ALL stdout + stderr (including third-party noise) to the same file
    redirect_streams_to_log(log_file_path)

    # Define the directories beneath output_dir
    files_dir, steps_dir, common_dir, workflow_path = get_paths("aggregate")

    valid, resolved_paths = check_inputs(args, logger)
    assert valid

    paths = {
        "files_dir": files_dir,
        "steps_dir": steps_dir,
        "common_dir": common_dir,
        "annotate_dir": os.path.abspath(args.annotate_dir),
        "denovo_dir": os.path.abspath(args.denovo_dir),
        "output_dir": os.path.abspath(args.output),
        # Full paths to the specific input files resolved by check_inputs
        # (consumed by workflow.smk as config["paths"][...])
        **resolved_paths,
    }
    
    # Create the output directory and snakemake subdirectory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    os.makedirs(os.path.join(args.output, "snakemake"), exist_ok=True)

    config = {
        "save_to_parquet" : args.save_to_parquet,
        "threads": args.threads,
        "mem_limit": args.mem,
        "debug": args.debug,
        "log": log_file_path,
        "paths": paths,

    }
    
    config_path = os.path.join(args.output, 'config_aggregate.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path, logger, log_file_path
