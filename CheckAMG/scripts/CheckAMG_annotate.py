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
    params_string = "checkamg annotate "
    for arg, value in vars(args).items():
        if arg == "command":
            continue
        if isinstance(value, bool):
            if value:
                params_string += f"--{arg} "
        else:
            if value is not None and value != "" and value != [] and value != "None":
                params_string += f"--{arg} {value} "
    return params_string
        
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
    
    # Create the output directory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    # Update with full paths
    for key, path in paths.items():
        paths[key] = os.path.join(args.output, path)

    # Ensure vmags is an absolute path
    vmags_abs = os.path.abspath(args.vmags) if args.vmags else None

    # List all genomic fasta files in vmags and construct their absolute paths
    if args.input_type == "nucl":
        vmag_fna_files = ' '.join(os.path.join(vmags_abs, fasta) for fasta in os.listdir(vmags_abs) if (fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".fna"))) if vmags_abs else ''
        vmag_faa_files = []
    if args.input_type == "prot":
        vmag_fna_files = []
        vmag_faa_files = ' '.join(os.path.join(vmags_abs, fasta) for fasta in os.listdir(vmags_abs) if (fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".faa"))) if vmags_abs else ''

    config = {
        "input_type": args.input_type,
        "input_single_contig_genomes": os.path.abspath(args.genomes) if (args.input_type == "nucl" and args.genomes) else "",
        "input_vmag_fastas" : vmag_fna_files,
        "input_single_contig_prots": os.path.abspath(args.proteins) if (args.input_type == "prot" and args.proteins) else "",
        "input_vmag_prots" : vmag_faa_files,
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
        "max_flank_length" : args.max_flank,
        "use_hallmark" : args.use_hallmark,
        "cov_fraction" : args.cov_fraction,
        "min_bitscore" : args.bit_score,
        "min_bitscore_fraction_heuristic" : args.bitscore_fraction_heuristic,
        "max_evalue" : args.evalue,
        "soft_keyword_bypass_scaling_factor" : args.scaling_factor,
    }
    
    config_path = os.path.join(args.output, 'config_annotate.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path

def run_snakemake(config_path, args):
    """Run the Snakemake pipeline using the generated config file."""

    logger = logging.getLogger()
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write(f"{ASCII}\n") # Write the ASCII art to the log file directly
            handler.flush() # Ensure the ASCII is written immediately
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
        log_file_path = os.path.join(os.path.abspath(args.output), 'PhAuxDB_build.log')
        print("========================================================================\n              The CheckAMG annotate pipeline is complete               \n========================================================================")
        with open(log_file_path, "a") as log:
            log.write("========================================================================\n              The CheckAMG annotate pipeline is complete               \n========================================================================\n")

    except subprocess.CalledProcessError as e:
        logger.error("CheckAMG annotate ended prematurely with an error!")