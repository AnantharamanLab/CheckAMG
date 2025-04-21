#!/usr/bin/env python3

import os
import shutil
from pathlib import Path
from pyfastatools import Parser, write_fasta
import sys
import logging
import subprocess
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

def set_memory_limit(limit_in_gb):
    current_os = platform.system()
    if current_os == "Linux":
        limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))

log_level = logging.DEBUG if snakemake.params.debug else logging.INFO
log_file = snakemake.params.log
logging.basicConfig(
    level=log_level,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(log_file, mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()

print("========================================================================\n             Step 19/22: Cluster filtered protein sequences             \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n             Step 19/22: Cluster filtered protein sequences             \n========================================================================\n")

def create_mmseqs_database(faa_path, db_path, log):
    """
    Converts a FASTA file(s) to an MMseqs2 database format, logging output.

    Args:
        faa_path (list): List of path to the input amino-acid FASTA file.
        db_path (str): Path to the output MMseqs2 database.
        log (file object): File object for logging the process.
    """
    cmd = ["mmseqs", "createdb", faa_path, db_path]
    result = subprocess.run(
        cmd, 
        stdout=log, 
        stderr=subprocess.PIPE, 
        text=True
    )
    if result.returncode != 0:
        log.write(result.stderr)
        log.flush()
        raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)

def run_mmseqs_search(db_path, result_path, tmp_dir, cov_fraction, search_sensitivity, num_cpus, log):
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = [
        "mmseqs", "search",
        db_path, db_path, result_path, tmp_dir,
        "-c", str(cov_fraction),
        "--start-sens", str(round(search_sensitivity * 0.5)), # Slightly slower, slightly more sensitive search, but still faster than just using the maximum sensitivity
        "--sens-steps", "3",
        "-s", str(search_sensitivity),
        "--cov-mode", "0",
        "--alignment-mode", "3",
        "--threads", str(num_cpus)
    ]
    result = subprocess.run(cmd, stdout=log, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        log.write(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, result.args)

def run_mmseqs_cluster(db_path, cluster_out, tmp_dir, min_seq_id, cov_fraction, sensitivity, num_cpus, log):
    os.makedirs(tmp_dir, exist_ok=True)
    cmd = [
        "mmseqs", "cluster",
        db_path, cluster_out, tmp_dir,
        "--cluster-mode", "0",
        "--min-seq-id", str(min_seq_id),
        "-c", str(cov_fraction),
        "-s", str(sensitivity),
        "--threads", str(num_cpus)
    ]
    result = subprocess.run(cmd, stdout=log, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        log.write(result.stderr)
        raise subprocess.CalledProcessError(result.returncode, result.args)

def run_mmseqs_convertalis(db_path, search_db, tsv_output, threads, log):
    cmd_convertalis = [
        "mmseqs", "convertalis",
        db_path, db_path,
        search_db, tsv_output,
        "--format-output", "query,target,fident,alnlen,qlen,tlen,qcov,tcov,bits",
        "--threads", str(threads)
    ]
    result = subprocess.run(cmd_convertalis, stdout=log, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.write(result.stderr)
        log.flush()
        raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)


def run_mmseqs_createtsv(db_path, cluster_output, tsv_output, threads, log):
    """
    Converts an MMseqs2 cluster output to a TSV format, logging output.

    Args:
        db_path (str): Path to the input MMseqs2 database.
        cluster_output (str): Path to the MMseqs2 cluster output.
        log (file object): File object for logging the process.
    """
    createtsv_command = [
        "mmseqs", "createtsv",
        db_path, cluster_output, tsv_output,
        "--threads", str(threads)
    ]
    result = subprocess.run(createtsv_command, stdout=log, stderr=subprocess.PIPE, text=True)

    if result.returncode != 0:
        log.write(result.stderr)
        log.flush()
        raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)

def filter_search_ses_by_bitscore(tsv, min_bits):
    logger.info(f"Filtering search results with bitscore >= {min_bits}...")
    cols = ["query","target","fident","alnlen","qlen","tlen","qcov","tcov","bits"]

    lf = (
        pl.scan_csv(
            tsv,
            separator="\t",
            has_header=False,
            new_columns=cols,
            schema_overrides={"query": pl.Utf8, "target": pl.Utf8, "bits": pl.Float64},
            infer_schema_length=0,
            ignore_errors=True
        )
        .select(["query", "target", "bits"])
        .filter(pl.col("bits") >= min_bits)
    )

    query_ids_lf  = lf.select(pl.col("query")).unique()
    target_ids_lf = lf.select(pl.col("target")).unique()

    queries = query_ids_lf.collect().get_column("query").to_list()
    targets = target_ids_lf.collect().get_column("target").to_list()

    seqs = set(queries)
    seqs.update(targets)

    logger.info(f"Collected {len(seqs):,} unique sequences passing the filter.")
    return seqs
    
def main():
    input_prot_dir = snakemake.params.input_prot_dir
    wdir = snakemake.params.wdir
    confidence_levels = snakemake.params.confidence_levels
    min_seq_id = snakemake.params.min_seq_id
    cov_fraction = snakemake.params.cov_fraction
    min_score = snakemake.params.min_score
    sensitivity = snakemake.params.sensitivity
    num_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    logger.info("Clustering of filtered proteins starting...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    
    # Build paths for each confidence level file
    input_faa_paths = [os.path.join(input_prot_dir, f"{level}_confidence_viral.faa") for level in confidence_levels]
    logger.debug(f"Input protein files: {input_faa_paths}")
    assert len(input_faa_paths) > 0, "No input protein files found. Please check the input directory."

    if len(input_faa_paths) == 1:
        filtered_prots_path = input_faa_paths[0]
        logger.debug(f"Using {filtered_prots_path} as the input protein file.")
    elif len(input_faa_paths) == 2:
        filtered_prots_path = os.path.join(input_prot_dir, f"{"_".join(confidence_levels)}_confidence_viral.faa")
        unique_records = {}
        for faa_path in input_faa_paths:
            # Parse each fasta file
            for record in Parser(faa_path):
                # Use the record's id as a key for uniqueness;
                # if a record with the same ID is already seen, it is skipped.
                unique_records[record.header.name] = record.seq
        # Write the union (non-redundant records) to the output fasta file.
        with open(filtered_prots_path, "w") as outfile:
            for record in unique_records.values():
                write_fasta(record, outfile)
        logger.debug(f"Combined protein file created: {filtered_prots_path}")
        logger.debug(f"Using {filtered_prots_path} as the input protein file.")
    elif len(input_faa_paths) == 3:
        filtered_prots_path = os.path.join(input_prot_dir, "filtered_proteins.faa")
        logger.debug(f"Using {filtered_prots_path} as the input protein file.")
    else:
        filtered_prots_path = None
        logger.error("More than 3 or less than 1 confidence levels provided. Please check the input files.")
        raise ValueError("More than 3 or less than 1 confidence levels provided. Please check the input files.")
    
    os.makedirs(wdir, exist_ok=True)
    
    log = os.path.join(wdir, "cluster_filtered_createdb.log")
    
    pre_db_path    = os.path.join(wdir, "filtered_prots.db")
    search_tsv     = os.path.join(wdir, "pre_search_res.tsv")
    search_res     = os.path.join(wdir, "search_res")
    search_tmp     = os.path.join(wdir, "search_tmp")
    
    db_path        = os.path.join(wdir, "filtered_by_bitscore.db")
    subset_fa      = os.path.join(wdir, "filtered_by_bitscore.faa")
    
    cluster_out    = os.path.join(wdir, "filtered_clusters")
    cluster_tmp    = os.path.join(wdir, "cluster_tmp")
    tsv_out        = os.path.join(wdir, "filtered_clusters.tsv")
    
    with open (log, "a") as log:
        # MMseqs2 database
        logger.info("Creating the MMseqs2 database from the input sequences...")
        create_mmseqs_database(filtered_prots_path, pre_db_path, log)
        
        # All-vs-all search with bit-score cutoff
        logger.info(f"Running MMseqs2 search on the input sequences...")
        run_mmseqs_search(pre_db_path, search_res, search_tmp, cov_fraction=cov_fraction, search_sensitivity=sensitivity, num_cpus=num_cpus, log=log)\
        
        # Convert search results to TSV for filtering
        logger.debug("Converting the MMseqs2 search output to TSV...")
        run_mmseqs_convertalis(pre_db_path, search_res, search_tsv, num_cpus, log)
        
        # Read the search result, calculate & filter bitscore
        logger.debug("Filtering the search results by bitscore...")
        filtered_seqs = filter_search_ses_by_bitscore(search_tsv, min_bits=min_score)

        # Write subset FASTA
        logger.debug(f"Writing the bitscore-filtered sequences to {filtered_prots_path}...")
        with open(subset_fa, "w") as out:
            for record in Parser(filtered_prots_path):
                if record.header.name in filtered_seqs:
                    write_fasta(record, out)

        # Create DB of subset
        logger.info(f"Creating an MMseqs2 database for the filtered sequences...")
        create_mmseqs_database(subset_fa, db_path, log)
        
        # Greedy Set Cover clustering on filtered graph
        logger.info(f"Clustering the bitscore-filtered sequences with MMseqs2 with a minimum coverage of {cov_fraction}...")
        run_mmseqs_cluster(db_path, cluster_out, cluster_tmp, min_seq_id, cov_fraction, sensitivity, num_cpus, log)

        # Convert clusters to TSV
        logger.info("Converting the MMseqs2 cluster output to TSV...")
        run_mmseqs_createtsv(db_path, cluster_out, tsv_out, num_cpus, log)
    
    # Clean up the working directory
    lookup_path = Path(str(pre_db_path) + ".lookup")
    for item in os.listdir(wdir):
        item_path = os.path.join(wdir, item)
        if os.path.abspath(item_path) != os.path.abspath(tsv_out):
            if os.path.abspath(item_path) != os.path.abspath(lookup_path):
                if os.path.isfile(item_path) or os.path.islink(item_path):
                    os.unlink(item_path)
                elif os.path.isdir(item_path):
                    shutil.rmtree(item_path)

    logger.info("Clustering of filtered proteins completed.")
    
if __name__ == "__main__":
    main()
