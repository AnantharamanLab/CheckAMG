#!/usr/bin/env python3

import yaml
import os
from pathlib import Path
from CheckAMG.scripts.common.runner_logging import setup_runner_logger, redirect_streams_to_log
from CheckAMG.scripts.common.set_up_snakemake import get_paths
from CheckAMG.scripts.common.resources import clamp_gpu_devices

def generate_config(args):
    """Generate a YAML config file based on provided arguments."""

    # Define the log file path and set up the logger
    log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_train.log')
    # Structured logger for explicit log.info() / log.debug() / log.warning() / log.error() calls
    logger = setup_runner_logger(log_file_path, debug=args.verbose)
    # Mirror ALL stdout + stderr (including third-party noise) to the same file
    redirect_streams_to_log(log_file_path)

    # Define the directories beneath output_dir
    files_dir, steps_dir, common_dir, workflow_path = get_paths("train_pst")

    effective_devices, gpu_clamp_warning = clamp_gpu_devices(args.accelerator, args.devices)

    train_embed_output = os.path.join(args.output, Path(args.train_file).with_suffix(".PST-EMBED.h5").name) if args.save_train_embed else None
    logdir = args.output

    paths = {
        # Snakemake
        "steps_dir": steps_dir,
        "common_dir": common_dir,
        # Inputs
        "train_file": os.path.abspath(args.train_file) if args.train_file else None,
        "test_file":  os.path.abspath(args.test_file) if args.test_file else None,
        "model_ckpt":  os.path.abspath(args.model_ckpt) if args.model_ckpt else None,
        # Outputs
        "output_dir": os.path.abspath(args.output),
        "train_embed_output": os.path.abspath(train_embed_output) if train_embed_output else None,
        "logdir": os.path.abspath(logdir),
    }

    # Create the output directory and snakemake subdirectory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    os.makedirs(os.path.join(args.output, "snakemake"), exist_ok=True)

    config = {
        # Parameters
        "batch_size": args.batch_size,
        "learning_rate": args.lr,
        "max_epochs": args.max_epochs,
        "dropout": args.dropout,
        "margin": args.margin,
        "knn": args.knn,
        "max_ap_pairs": args.max_ap_pairs,
        "positive_mining_strategy": args.positive_mining_strategy,
        "negative_mining_strategy": args.negative_mining_strategy,
        "opposite_class_negatives": args.opposite_class_negatives,
        "class_weighting": args.class_weighting,
        "context_size": args.context_size,

        # Resources
        "accelerator": args.accelerator,
        "devices": effective_devices,
        "threads": effective_devices,
        "mem_limit": args.mem,

        # Other
        "debug": args.verbose,
        "log": log_file_path,
        "seed": args.seed,

        # Paths
        "paths": paths,

    }

    config_path = os.path.join(args.output, 'config_train.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path, logger, log_file_path, gpu_clamp_warning
