#!/usr/bin/env python3

import os
import sys
import resource
import platform
import logging
from pyfastatools import Parser, write_fasta
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
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

print("========================================================================\n      Step 17/22: Extract proteins from filtered protein clusters       \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n      Step 17/22: Extract proteins from filtered protein clusters       \n========================================================================\n")

def extract_sequences_from_single_contig_fasta(single_contig_fasta, protein_ids_set):
    """
    Uses pyfastatools Parser for the single-contig FASTA file.
    Returns a dictionary {protein_id: record} for protein_ids in protein_ids_set.
    This streaming approach avoids building a full index and is faster
    when a large fraction of records need to be extracted.
    """
    extracted = {}
    for record in Parser(single_contig_fasta):
        if record.header.name in protein_ids_set:
            extracted[record.header.name] = record
    return extracted

def extract_sequences_from_vmag_fasta(faa_path, protein_ids_set):
    """
    Reads a vMAG FASTA file and returns a dictionary of {protein_id: record}
    for proteins in protein_ids_set.
    """
    extracted = {}
    for record in Parser(faa_path):
        if record.header.name in protein_ids_set:
            extracted[record.header.name] = record
    return extracted

def parallel_extract_vmag_fastas(vmag_faa_paths, protein_ids_set, max_workers):
    """
    Processes the list of vMAG FASTA files in parallel (using threads) and merges
    the results into a nonredundant dictionary.
    """
    extracted = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(extract_sequences_from_vmag_fasta, path, protein_ids_set): path
                   for path in vmag_faa_paths}
        for future in futures:
            result = future.result()
            extracted.update(result)
    return extracted

def write_extracted_proteins_to_fasta(extracted, out_fasta_path):
    """
    Writes the extracted protein sequences to a single FASTA file.
    """
    with open(out_fasta_path, "w") as out_fasta:
        for protein_id, record in extracted.items():
            write_fasta(record, out_fasta)
    logger.info(f"Extracted proteins written to {out_fasta_path}")
    
def main():
    outdir_path = snakemake.params.outdir
    aux_scores_path = snakemake.params.aux_scores
    protein_dir = snakemake.params.protein_dir
    vmag_proteins_subdir = snakemake.params.vmag_proteins_subdir
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    logger.info("Extracting proteins that were assigned auxiliary scores...")

    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)
    
    # Determine FASTA file paths:
    # Assume one large FASTA file is in protein_dir as 'single_contig_proteins.faa'
    single_contig_fasta = None
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        single_contig_fasta = os.path.join(protein_dir, 'single_contig_proteins.faa')
    # All vMAG FASTA files from the vMAG proteins subdirectory:
    vmag_fastas = []
    if os.path.exists(vmag_proteins_subdir) and os.path.isdir(vmag_proteins_subdir):
        vmag_fastas = [os.path.join(vmag_proteins_subdir, f)
                       for f in os.listdir(vmag_proteins_subdir)
                       if f.endswith('.faa')]
    if single_contig_fasta is None and not vmag_fastas:
        logger.error("No protein files found. Exiting.")
        sys.exit(1)
    logger.debug(f"Found single-contig FASTA: {single_contig_fasta}" if single_contig_fasta else "No single-contig FASTA found.")
    logger.debug(f"Found {len(vmag_fastas):,} vMAG FASTA files.")

    # Read auxiliary scores table to get the set of protein IDs with valid aux scores
    aux_scores = pl.read_csv(aux_scores_path, separator="\t", has_header=True)
    # aux_scores = aux_scores.filter(
    #             (pl.col("quantile_weighted").is_not_null()) &
    #             (pl.col("aux_score_weighted").is_not_null())
    #         )
    protein_ids = set(aux_scores.select("protein").to_series().to_list())
    logger.info(f"Extracting {len(protein_ids):,} proteins with valid auxiliary scores from FASTA files.")

    # Extract proteins from the single-contig FASTA
    extracted = {}
    if single_contig_fasta:
        logger.info(f"Extracting sequences from the single-contig FASTA file {single_contig_fasta}")
        extracted.update(extract_sequences_from_single_contig_fasta(single_contig_fasta, protein_ids))
    # Extract proteins from the vMAG FASTA files in parallel
    if vmag_fastas:
        logger.info(f"Extracting sequences from individual vMAG FASTA files in {vmag_proteins_subdir}")
        vmag_extracted = parallel_extract_vmag_fastas(vmag_fastas, protein_ids, max_workers=threads)
        extracted.update(vmag_extracted)

    # Write the extracted proteins to a single FASTA file
    out_fasta_path = os.path.join(outdir_path, "filtered_proteins.faa")
    write_extracted_proteins_to_fasta(extracted, out_fasta_path)
    # Check if the output file was created successfully
    if os.path.exists(out_fasta_path) and os.path.getsize(out_fasta_path) > 0:
        logger.debug(f"Filtered proteins successfully extracted to {out_fasta_path}")
    else:
        logger.error(f"Failed to create output file: {out_fasta_path}")
        raise RuntimeError(f"Failed to create output file: {out_fasta_path}")
        
    logger.info("Filtered protein extraction complete.")

if __name__ == "__main__":
    main()