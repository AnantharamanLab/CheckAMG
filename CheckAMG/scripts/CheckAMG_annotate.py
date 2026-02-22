#!/usr/bin/env python3

import subprocess
import yaml
import os
import platform
from importlib.resources import files as resource_files
import logging
from CheckAMG.scripts.checkAMG_ASCII import ASCII

# Access the packaged files and scripts directories
scripts_dir = os.path.abspath(os.path.dirname(__file__))
try:
    files_dir = str(resource_files("CheckAMG").joinpath("files"))
except ModuleNotFoundError as e:
    raise RuntimeError("Package data not found. Is 'CheckAMG/files' included in your package?") from e

def log_command_args(args):
    parts = ["checkamg", "annotate"]

    for key, value in vars(args).items():
        if key in {"command", "func", "_cli_argv"}:
            continue
        if callable(value):
            continue

        dashed = key.replace("_", "-")
        flag = f"--{dashed}"
        no_flag = f"--no-{dashed}"

        # Always print booleans as --flag / --no-flag
        if isinstance(value, bool):
            parts.append(flag if value else no_flag)
            continue

        # Skip empty values
        if value is None or value == "" or value == [] or value == "None":
            continue

        # Lists/tuples
        if isinstance(value, (list, tuple)):
            for v in value:
                if v is None or v == "" or v == "None":
                    continue
                parts.extend([flag, str(v)])
        else:
            parts.extend([flag, str(value)])

    return " ".join(parts)
        
def setup_logger(log_file_path, debug):
    """Sets up the logger to write to both console and a file."""
    # Create a custom logger
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG if debug else logging.INFO)

    # Remove any existing handlers to avoid duplicate logs
    if logger.hasHandlers():
        logger.handlers.clear()

    # Create handlers for both console and file
    console_handler = logging.StreamHandler()
    file_handler = logging.FileHandler(log_file_path)

    # Set log format without milliseconds
    formatter = logging.Formatter(
        '%(asctime)s | %(levelname)s | %(message)s', 
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(formatter)
    file_handler.setFormatter(formatter)

    # Add handlers to the logger
    logger.addHandler(console_handler)
    logger.addHandler(file_handler)

    return logger

def create_output_dir(output_dir):
    """Create the output directory if it doesn't already exist."""
    os.makedirs(output_dir, exist_ok=True)

def generate_config(args):
    """Generate a YAML config file based on provided arguments."""

    # Define the log file path and set up the logger
    log_file_path = os.path.abspath(os.path.join(args.output, 'CheckAMG_annotate.log'))
    
    # Initialize the logger
    logger = setup_logger(log_file_path, args.debug)
    
    # Define the directories beneath output_dir
    paths = {
        "scripts_dir" : scripts_dir,
        "files_dir" : files_dir,
        "db_dir" : os.path.abspath(args.db_dir),
        "output_dir" : os.path.abspath(args.output)
    }
    
    # Create the output directory and snakemake subdirectory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    os.makedirs(os.path.join(args.output, "snakemake"), exist_ok=True)
    
    # Update with full paths
    for key, path in paths.items():
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
        "filter_presets" : args.filter_presets,
        "filter_nonviral_regions" : args.filter_ambig_regions,
        "filter_avg_arrays" : args.filter_avg_arrays,
        "avg_array_len_limit" : args.avg_array_limit
    }
    
    config_path = os.path.join(args.output, 'config_annotate.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path

def run_snakemake(config_path, args, checkamg_version):
    """Run the Snakemake pipeline using the generated config file."""

    logger = logging.getLogger()
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write(f"{ASCII}\n") # Write the ASCII art to the log file directly
            handler.flush() # Ensure the ASCII is written immediately
    logger.info(f"CheckAMG version {checkamg_version}")
    logger.info("Starting CheckAMG annotate...")
    
    current_os = platform.system()
    if current_os == "Darwin":
        logger.warning(
            f"The detected OS is {current_os}, which means no hard memory limit can be set. "
            "This should be fine, but there may be problems/crashes if you are working with very large inputs that exceed your available memory."
        )
    elif current_os == "Windows":
        logger.error(
            f"The detected OS is {current_os}, which is not supported. Exiting..."
        )
        raise OSError("Windows is not supported for CheckAMG.")
    
    logger.info(f"Command issued: {log_command_args(args)}")
    logger.debug(f"The input type is {args.input_type}")

    # Execute the snakemake workflow
    try:
        if args.debug:
            snakemake_command = [
                "snakemake", "--snakefile", os.path.join(scripts_dir, "CheckAMG_annotate.smk"),
                "--nolock", "--configfile", config_path, "--directory", args.output, "--cores",
                str(args.threads) , "--rerun-triggers", "input",
                "--keep-incomplete",
                "--ignore-incomplete", # Debugging, for when the order of rules have been modified but old outputs were saved
                "--verbose", "all"
            ]
        else:
            snakemake_command = [
                "snakemake", "--snakefile", os.path.join(scripts_dir, "CheckAMG_annotate.smk"),
                "--nolock", "--configfile", config_path, "--directory", args.output, "--cores",
                str(args.threads) , "--rerun-triggers", "input",
                "--keep-incomplete",
                "--ignore-incomplete", # Debugging, for when the order of rules have been modified but old outputs were saved
                "--quiet", "all"
            ]
        subprocess.run(snakemake_command, check=True)
        log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_annotate.log')
        print("========================================================================\n              The CheckAMG annotate pipeline is complete               \n========================================================================")
        with open(log_file_path, "a") as log:
            log.write("========================================================================\n              The CheckAMG annotate pipeline is complete               \n========================================================================\n")

    except subprocess.CalledProcessError as e:
        logger.error("CheckAMG annotate ended prematurely with an error!")