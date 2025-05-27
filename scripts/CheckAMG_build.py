import subprocess
import yaml
import os
import platform
import logging
from scripts import checkAMG_ASCII

# Automatically determine the base directory of the current script
scripts_dir = os.path.abspath(os.path.dirname(__file__))
base_dir = os.path.abspath(os.path.join(scripts_dir, '..'))
files_dir = os.path.join(base_dir, 'files')

def log_command_args(args):
    params_string = "CheckAMG build "
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
    
    # Define the log file path
    log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_build.log')

    # Initialize the logger
    logger = setup_logger(log_file_path, args.debug)
    
    # Define the directories beneath output_dir
    paths = {
        "scripts_dir": scripts_dir,
        "files_dir": files_dir,
        "db_dir": os.path.abspath(args.db_dir),
        "output_dir": os.path.abspath(args.output)
    }
    # Create the output directory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)

    # Update with full paths
    for key, path in paths.items():
        paths[key] = os.path.join(args.output, path)
    
    # Ensure vmags is an absolute path
    vmags_abs = os.path.abspath(args.vmags) if args.vmags else None

    # List all genomic fasta files in vmags and construct their absolute paths
    vmag_fasta_files = ' '.join(os.path.join(vmags_abs, fasta) for fasta in os.listdir(vmags_abs) if (fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".fna"))) if vmags_abs else ''
    
    # Split the comma-separated list of cluster ranks and validate
    valid_ranks = {"class", "family", "genus", "species"}
    cluster_ranks = args.cluster_ranks.split(",")
    for rank in cluster_ranks:
        if rank not in valid_ranks:
            raise ValueError(f"Invalid taxonomic rank specified: {rank}. Valid choices are {valid_ranks}")
    confidence_levels = args.confidence_levels.split(",")
    for level in confidence_levels:
        if level not in ["high", "medium", "low"]:
            raise ValueError(f"Invalid confidence level specified: {level}. Valid choices are ['high', 'medium', 'low']")
        
    try:
        min_lscores = [float(lscore) for lscore in (args.min_window_avg_Lscores).split(",")]
        min_lscores = {"KEGG": min_lscores[0], "Pfam": min_lscores[1], "PHROG": min_lscores[2]}
    except ValueError:
        raise ValueError(f"Invalid format for --min_window_avg_Lscores: {args.min_window_avg_Lscores} It should be a comma-separated list of floating points.")
    if len(min_lscores) != 3:
        raise ValueError(f"Please provide exactly three values for --min_window_avg_Lscores (KEGG, Pfam, PHROG).")
    
        
    config = {
        "input_single_contig_genomes": os.path.abspath(args.genomes),
        "input_vmag_fastas": vmag_fasta_files,
        "min_cds" : args.min_orf,
        "min_len": args.min_len,
        "threads": args.threads,
        "mem_limit": args.mem,
        "debug": args.debug,
        "log": log_file_path,
        "mmseqs_params": {
            "sensitivity": args.sensitivity,
            "min_seq_id": args.min_seq_id,
            "cov_fraction": args.cov_fraction
        },
        "paths": paths,
        "cluster_ranks": cluster_ranks,
        "confidence_levels": confidence_levels,
        "exclude_singletons": args.exclude_singletons,
        "annotation_percent_threshold": args.min_annot,
        "min_window_avg_lscores": min_lscores,
        "window_size": args.window_size,
        "minimum_flank_vscore": args.min_flank_Vscore,
        "max_flank_length": args.max_flank,
        "use_hallmark": args.use_hallmark,
        "cov_fraction": args.cov_fraction,
        "hmm_acc_prefix": args.hmm_prefix
    }

    logger.debug(f"vMAG fasta files: {vmag_fasta_files}")
    logger.debug(f"Single contig genome files: {args.genomes}")
    
    config_path = os.path.join(args.output, 'config_build.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    
    return config_path

def run_snakemake(config_path, args):
    """Run the Snakemake pipeline using the generated config file."""

    logger = logging.getLogger()
    for handler in logger.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.stream.write(f"{checkAMG_ASCII.ASCII}\n") # Write the ASCII art to the log file directly
            handler.flush() # Ensure the ASCII is written immediately
    logger.info("Starting CheckAMG build...")
    
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

    # Execute the snakemake workflow
    try:
        if args.debug:
            snakemake_command = [
                "snakemake", "--snakefile", os.path.join(scripts_dir, "CheckAMG_build.smk"),
                "--nolock", "--configfile", config_path, "--directory", args.output, "--cores",
                str(args.threads), "--rerun-triggers", "input",
                "--keep-incomplete",
                "--ignore-incomplete", # Debugging, for when the order of rules have been modified but old outputs were saved
                "--verbose", "all",
            ]
        else:
            snakemake_command = [
                "snakemake", "--snakefile", os.path.join(scripts_dir, "CheckAMG_build.smk"),
                "--nolock", "--configfile", config_path, "--directory", args.output, "--cores",
                str(args.threads), "--rerun-triggers", "input",
                "--keep-incomplete",
                "--ignore-incomplete", # Debugging, for when the order of rules have been modified but old outputs were saved
                "--quiet", "all"
            ]
        
        subprocess.run(snakemake_command, check=True) # check=True raises an error if the command fails
        logger.info("CheckAMG build complete.")
        
    except subprocess.CalledProcessError as e:
        logger.error("CheckAMG build ended prematurely with an error!")