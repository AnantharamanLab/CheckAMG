#!/usr/bin/env python3

import yaml
import os
import glob
import json
from CheckAMG.scripts.common.runner_logging import setup_runner_logger, redirect_streams_to_log
from CheckAMG.scripts.common.set_up_snakemake import get_paths
from CheckAMG.scripts.common.resources import clamp_gpu_devices
from CheckAMG.scripts.common.db_paths import resolve_denovo_db


def _discover_reference_files(db_dir):
    """Locate de-novo reference files within a CheckAMG de-novo database directory.

    The de-novo DB is its own directory. --db-dir may point directly at it, or at a
    parent that contains it as an immediate subdirectory (e.g. when the annotate and
    de-novo DBs are downloaded side by side). The legacy 'checkamg_pst' subdirectory
    is still honored as a fallback for older combined-DB layouts.
    """
    def _subdirs(d):
        try:
            return [os.path.join(d, s) for s in sorted(os.listdir(d)) if os.path.isdir(os.path.join(d, s))]
        except OSError:
            return []

    ordered = [db_dir, os.path.join(db_dir, "checkamg_pst"), *_subdirs(db_dir)]
    search_dirs = []
    for d in ordered:
        if os.path.isdir(d) and d not in search_dirs:
            search_dirs.append(d)

    def first(patterns, exclude=()):
        for d in search_dirs:
            for pat in patterns:
                for hit in sorted(glob.glob(os.path.join(d, pat))):
                    base = os.path.basename(hit)
                    if any(ex in base for ex in exclude):
                        continue
                    return os.path.abspath(hit)
        return None

    return {
        "model_ckpt": first(["*.ckpt"]),
        "train_index_file": first(["*.index.faiss", "*.faiss"]),
        "train_labels_file": first(["*.labels.h5"]),
        "train_embed_file": first(["*.PST-EMBED.h5"], exclude=[".labels.", ".index."]),
        "train_data_file": first(["*.graphfmt.h5"]),
    }


def resolve_reference_data(args, logger):
    """Fill in the model checkpoint and training reference files from --db-dir
    when they are not explicitly provided.

    Reference priority: (train index + train labels) > train embeddings > train
    data. The model checkpoint is always resolved this way. Explicit CLI args
    always take precedence over discovery.
    """
    explicit_train = any(
        [args.train_data_file, args.train_embed_file, args.train_index_file, args.train_labels_file]
    )

    if not args.db_dir:
        if not args.model_ckpt:
            logger.error("No --db-dir and no --model-ckpt provided; cannot locate the de-novo model checkpoint.")
            raise FileNotFoundError("Missing --model-ckpt (and no --db-dir to discover it from).")
        if not explicit_train:
            logger.error(
                "No --db-dir and no training reference provided. Provide one of "
                "(--train-index-file and --train-labels-file), --train-embed-file, or --train-data-file."
            )
            raise FileNotFoundError("Missing de-novo training reference (and no --db-dir to discover it from).")
        return

    args.db_dir = resolve_denovo_db(args.db_dir)
    found = _discover_reference_files(args.db_dir)

    if not args.model_ckpt:
        if found["model_ckpt"]:
            args.model_ckpt = found["model_ckpt"]
            logger.debug(f"Using model checkpoint from --db-dir: {args.model_ckpt}")
        else:
            logger.error(f"Could not find a model checkpoint (*.ckpt) under --db-dir '{args.db_dir}'.")
            raise FileNotFoundError(f"No model checkpoint (*.ckpt) found under '{args.db_dir}'.")

    if not explicit_train:
        if found["train_index_file"] and found["train_labels_file"]:
            args.train_index_file = found["train_index_file"]
            args.train_labels_file = found["train_labels_file"]
            logger.debug(
                f"Using training index/labels from --db-dir: {args.train_index_file}, {args.train_labels_file}"
            )
        elif found["train_embed_file"]:
            args.train_embed_file = found["train_embed_file"]
            logger.debug(f"Using training embeddings from --db-dir: {args.train_embed_file}")
        elif found["train_data_file"]:
            args.train_data_file = found["train_data_file"]
            logger.debug(f"Using training data from --db-dir: {args.train_data_file}")
        else:
            logger.error(
                f"Could not find any de-novo training reference under --db-dir '{args.db_dir}' "
                "(expected an index+labels pair, a PST embeddings file, or a graph-formatted training file)."
            )
            raise FileNotFoundError(f"No de-novo training reference found under '{args.db_dir}'.")

def check_inputs(args, logger):
    if args.query_proteins:
        if os.path.isfile(args.query_proteins):
            query_input = "faa"
        else:
            logger.error(f"Provided --query-proteins path '{args.query_proteins}' does not exist or is not a valid file.")
            raise ValueError(f"Provided --query-proteins path '{args.query_proteins}' does not exist or is not a valid file.")
    if args.query_bin_proteins:
        if os.path.isdir(args.query_bin_proteins):
            query_input = "faa"
        else:
            logger.error(f"Provided --query-bin-proteins path '{args.query_bin_proteins}' does not exist or is not a valid directory.")
            raise ValueError(f"Provided --query-bin-proteins path '{args.query_bin_proteins}' does not exist or is not a valid directory.")
    if args.query_contigs:
        if os.path.isfile(args.query_contigs):
            query_input = "fna"
        else:
            logger.error(f"Provided --query-contigs path '{args.query_contigs}' does not exist or is not a valid file.")
            raise ValueError(f"Provided --query-contigs path '{args.query_contigs}' does not exist or is not a valid file.")
    if args.query_bins:
        if os.path.isdir(args.query_bins):
            query_input = "fna"
        else:
            logger.error(f"Provided --query-bins path '{args.query_bins}' does not exist or is not a valid directory.")
            raise ValueError(f"Provided --query-bins path '{args.query_bins}' does not exist or is not a valid directory.")
    
    if (args.query_proteins or args.query_bin_proteins) and (args.query_contigs or args.query_bins):
        logger.error("Conflicting query input: both protein and nucleotide query inputs provided. Please provide only one type of query input (protein or nucleotide).")
        raise ValueError("Conflicting query input: both protein and nucleotide query inputs provided. Please provide only one type of query input (protein or nucleotide).")

    if args.query_file_esm2:
        if os.path.isfile(args.query_file_esm2):
            query_input = "esm2"
        else:
            logger.error(f"Provided --query-file-esm2 path '{args.query_file_esm2}' does not exist or is not a valid file.")
            raise ValueError(f"Provided --query-file-esm2 path '{args.query_file_esm2}' does not exist or is not a valid file.")
    if args.query_file:
        if os.path.isfile(args.query_file):
            query_input = "pst_embed"
        else:
            logger.error(f"Provided --query-file path '{args.query_file}' does not exist or is not a valid file.")
            raise ValueError(f"Provided --query-file path '{args.query_file}' does not exist or is not a valid file.")
    
    if args.query_file_esm2 and args.query_file:
        logger.error("Conflicting query input: both --query-file-esm2 and --query-file provided. Please provide only one of these inputs.")
        raise ValueError("Conflicting query input: both --query-file-esm2 and --query-file provided. Please provide only one of these inputs.")

    if ((args.query_proteins or args.query_bin_proteins) or (args.query_contigs or args.query_bins)) and (args.query_file_esm2 or args.query_file):
        logger.error("Conflicting query input: both raw sequence inputs and precomputed embedding inputs provided. Please provide only one type of query input (raw sequences or precomputed embeddings).")
        raise ValueError("Conflicting query input: both raw sequence inputs and precomputed embedding inputs provided. Please provide only one type of query input (raw sequences or precomputed embeddings).")

    if not (args.query_proteins or args.query_bin_proteins or args.query_contigs or args.query_bins or args.query_file_esm2 or args.query_file):
        logger.error("No valid query input file provided. Please provide one of: (--query-proteins and/or --query-bin-proteins), --query-file-esm2, or --query-file.")
        raise ValueError("No valid query input file provided. Please provide one of (--query-proteins and/or --query-bin-proteins), --query-file-esm2, or --query-file.")

    if getattr(args, "esm2_ckpt_dir", None):
        if not os.path.isdir(args.esm2_ckpt_dir):
            logger.error(f"Provided --esm2-ckpt-dir path '{args.esm2_ckpt_dir}' does not exist or is not a valid directory.")
            raise ValueError(f"Provided --esm2-ckpt-dir path '{args.esm2_ckpt_dir}' does not exist or is not a valid directory.")
        required = ["esm2_t30_150M_UR50D.pt", "esm2_t30_150M_UR50D-contact-regression.pt"]
        missing = [f for f in required if not os.path.isfile(os.path.join(args.esm2_ckpt_dir, f))]
        if missing:
            logger.error(f"Provided --esm2-ckpt-dir '{args.esm2_ckpt_dir}' is missing required ESM2 checkpoint files: {missing}")
            raise FileNotFoundError(f"Provided --esm2-ckpt-dir '{args.esm2_ckpt_dir}' is missing required ESM2 checkpoint files: {missing}")

    return query_input


def generate_config(args):
    """Generate a YAML config file based on provided arguments."""

    # Define the log file path and set up the logger
    log_file_path = os.path.join(os.path.abspath(args.output), 'CheckAMG_denovo.log')
    # Structured logger for explicit log.info() / log/debug() / log.warning() / log.error() calls
    logger = setup_runner_logger(log_file_path, debug=args.debug)
    # Mirror ALL stdout + stderr (including third-party noise) to the same file
    redirect_streams_to_log(log_file_path)

    # Define the directories beneath output_dir
    files_dir, steps_dir, common_dir, workflow_path = get_paths("denovo")

    # Load the PST probability -> confidence cutoffs that are bundled with the package.
    pst_thresholds_path = os.path.join(files_dir, "pst_thresholds.json")
    logger.debug(f"Reading PST confidence thresholds from {pst_thresholds_path}")
    with open(pst_thresholds_path) as f:
        pst_thresholds = json.load(f)

    # Ensure bins is an absolute path
    bins_abs = os.path.abspath(args.query_bins) if args.query_bins else None
    bins_prot_abs = os.path.abspath(args.query_bin_proteins) if args.query_bin_proteins else None

    query_input = check_inputs(args, logger)

    resolve_reference_data(args, logger)

    effective_devices, gpu_clamp_warning = clamp_gpu_devices(args.accelerator, args.devices)

    # List all genomic fasta files in bins and construct their absolute paths
    if query_input == "fna":
        bin_fna_files = ' '.join(
            os.path.join(bins_abs, fasta) for fasta in os.listdir(bins_abs) if (
                fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".fna") or\
                fasta.endswith(".fasta.gz") or fasta.endswith(".fa.gz") or fasta.endswith(".fna.gz")
                )
            ) if bins_abs else ''
        bin_faa_files = []
    elif query_input == "faa":
        bin_fna_files = []
        bin_faa_files = ' '.join(
            os.path.join(bins_prot_abs, fasta) for fasta in os.listdir(bins_prot_abs) if (
                fasta.endswith(".fasta") or fasta.endswith(".fa") or fasta.endswith(".faa") or\
                fasta.endswith(".fasta.gz") or fasta.endswith(".fa.gz") or fasta.endswith(".faa.gz")
                )
            ) if bins_prot_abs else ''
    else:
        bin_fna_files = []
        bin_faa_files = []
        
    paths = {
        "steps_dir": steps_dir,
        "common_dir": common_dir,

        "db_dir": os.path.abspath(args.db_dir) if args.db_dir else None,
        "model_ckpt": os.path.abspath(args.model_ckpt) if args.model_ckpt else None,
        "esm2_ckpt_dir": os.path.abspath(args.esm2_ckpt_dir) if getattr(args, "esm2_ckpt_dir", None) else None,

        "input_single_contig_fna": os.path.abspath(args.query_contigs) if (query_input == "fna" and args.query_contigs) else "",
        "input_bin_fna" : bin_fna_files,
        "input_bin_fna_dir" : bins_abs if query_input == "fna" else None,
        "input_single_contig_faa": os.path.abspath(args.query_proteins) if (query_input == "faa" and args.query_proteins) else "",
        "input_bin_faa" : bin_faa_files,
        "input_bin_faa_dir" : bins_prot_abs if query_input == "faa" else None,

        "query_fasta":  os.path.abspath(args.fasta_file) if args.fasta_file else None,
        "query_esm2_file":  os.path.abspath(args.query_file_esm2) if args.query_file_esm2 else None,
        "query_file_pst":  os.path.abspath(args.query_file) if args.query_file else None,

        "train_data_file":  os.path.abspath(args.train_data_file) if args.train_data_file else None,
        "train_embed_file":  os.path.abspath(args.train_embed_file) if args.train_embed_file else None,
        "train_index_file":  os.path.abspath(args.train_index_file) if args.train_index_file else None,
        "train_labels_file":  os.path.abspath(args.train_labels_file) if args.train_labels_file else None,

        "output_dir": os.path.abspath(args.output),
    }
    
    # Create the output directory and snakemake subdirectory
    os.makedirs(os.path.abspath(args.output), exist_ok=True)
    os.makedirs(os.path.join(args.output, "snakemake"), exist_ok=True)

    config = {
        "query_input": query_input,

        "pst_thresholds": pst_thresholds,

        "knn": args.knn,
        "batch_size": args.batch_size,
        "nn_batch_size": args.nn_batch_size,

        "num_index_cells": args.num_index_cells,
        "num_probe_cells": args.num_probe_cells,
        "seed": args.seed,

        "accelerator": args.accelerator,
        "devices": effective_devices,
        "batch_size": args.batch_size,

        "threads": effective_devices,
        "mem_limit": args.mem,
        "debug": args.debug,
        "log": log_file_path,
        "paths": paths,

    }
    
    config_path = os.path.join(args.output, 'config_denovo.yaml')
    with open(config_path, 'w') as file:
        yaml.dump(config, file, default_flow_style=False)
    return config_path, logger, log_file_path, gpu_clamp_warning
