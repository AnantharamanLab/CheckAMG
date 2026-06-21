#!/usr/bin/env python3

import yaml
import os
from CheckAMG.scripts.common.runner_logging import setup_runner_logger, redirect_streams_to_log
from CheckAMG.scripts.common.set_up_snakemake import get_paths
from CheckAMG.scripts.common.db_paths import resolve_annotate_db

def generate_config(args):
    """Generate a YAML config file based on provided arguments."""

    # Define the log file path and set up the logger
    log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_annotate.log')
    logger = setup_runner_logger(log_file_path, args.debug)
    # Structured logger for explicit log.info() / log/debug() / log.warning() / log.error() calls
    logger = setup_runner_logger(log_file_path, debug=args.debug)
    # Mirror ALL stdout + stderr (including third-party noise) to the same file
    redirect_streams_to_log(log_file_path)
        
    # Define the directories beneath output_dir
    files_dir, steps_dir, common_dir, workflow_path = get_paths("annotate")

    paths = {
        "files_dir": files_dir,
        "steps_dir": steps_dir,
        "common_dir": common_dir,
        "db_dir": resolve_annotate_db(args.db_dir),
        "output_dir": os.path.abspath(args.output),
    }
    
    # Create the output directory and snakemake subdirectory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    os.makedirs(os.path.join(args.output, "snakemake"), exist_ok=True)
    
    # Update with full paths
    for key, path in list(paths.items()):
        if not os.path.isabs(path):
            paths[key] = os.path.join(args.output, path)

    # Ensure bins is an absolute path
    bins_abs = os.path.abspath(args.input_bins) if args.input_bins else None
    bins_prot_abs = os.path.abspath(args.input_bin_proteins) if args.input_bin_proteins else None

    # List all genomic fasta files in bins and construct their absolute paths
    if args.input_type == "nucl":
        bin_fna_files = ' '.join(
            os.path.join(bins_abs, fasta) for fasta in os.listdir(bins_abs) if (
                fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".fna") or\
                fasta.endswith(".fasta.gz") or fasta.endswith(".fa.gz") or fasta.endswith(".fna.gz")
                )
            ) if bins_abs else ''
        bin_faa_files = []
    if args.input_type == "prot":
        bin_fna_files = []
        bin_faa_files = ' '.join(
            os.path.join(bins_prot_abs, fasta) for fasta in os.listdir(bins_prot_abs) if (
                fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".faa") or\
                fasta.endswith(".fasta.gz") or fasta.endswith(".fa.gz") or fasta.endswith(".faa.gz")
                )
            ) if bins_prot_abs else ''

    config = {
        "input_type": args.input_type,
        "input_single_contigs": os.path.abspath(args.input_contigs) if (args.input_type == "nucl" and args.input_contigs) else "",
        "input_bins" : bin_fna_files,
        "input_single_contig_prots": os.path.abspath(args.input_proteins) if (args.input_type == "prot" and args.input_proteins) else "",
        "input_bin_prots" : bin_faa_files,
        "min_cds" : args.min_orf,
        "min_len": args.min_len,
        "threads": args.threads,
        "mem_limit": args.mem,
        "debug": args.debug,
        "log": log_file_path,
        "paths": paths,
        "annotation_percent_threshold" : args.min_annot,
        "window_size" : args.window_size,
        "minimum_flank_vscore" : args.min_flank_Vscore,
        "minimum_window_avg_vlscore" : args.min_window_avg_VL_score,
        "use_hallmark" : args.use_hallmark,
        "cov_fraction" : args.cov_fraction,
        "min_bitscore" : args.bit_score,
        "min_bitscore_fraction_heuristic" : args.bitscore_fraction_heuristic,
        "max_evalue" : args.evalue,
        "keep_full_hmm_results" : args.keep_full_hmm_results,
        "save_to_parquet" : args.save_to_parquet,
        "min_amg_weight" : args.min_amg_weight,
        "filter_presets" : args.filter_presets,
        "filter_nonviral_regions" : args.filter_ambig_regions,
        "filter_avg_arrays" : args.filter_avg_arrays,
        "avg_array_len_limit" : args.avg_array_limit
    }
    
    config_path = os.path.join(args.output, 'config_annotate.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path, logger, log_file_path

