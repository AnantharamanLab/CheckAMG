#!/usr/bin/env python3

import os
import sys
import shutil
from pathlib import Path
import logging
import resource
import platform
from metapyrodigal.load_balancer import SingleFileLoadBalancer, load_balancer
from metapyrodigal.orf_finder import OrfFinder
from typing import Optional

def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    except (ValueError, OSError, AttributeError) as e:
        logger.warning(f"Unable to set memory limit. Error: {e}")

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
print("========================================================================\n   Step 3/11: Predict and translate ORFs in genomes with pyrodigal-GV  \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n   Step 3/11: Predict and translate ORFs in genomes with pyrodigal-GV  \n========================================================================\n")

def run_metapyrodigal(input_fasta: Optional[str],
                      vmag_fasta_dir: Optional[str],
                      output_dir: str,
                      vmag_proteins_subdir: str,
                      single_contig_prots: str,
                      threads: int) -> None:
    """
    Run metapyrodigal-GV on provided single contig and/or vMAG fasta files.
    Depending on the provided inputs, this function will execute metapyrodigal on
    vMAG fasta files located in a directory or on a single contig genome fasta file.
    
    Arguments:
        input_fasta: Path to the single contig genome fasta file.
        vmag_fasta_dir: Path to the directory containing vMAG fasta files.
        output_dir: Directory where output files will be saved.
        vmag_proteins_subdir: Subdirectory name for vMAG protein outputs.
        single_contig_prots: Output file path for single contig protein results.
        threads: Number of CPU threads to leverage.
    """
    orf_finder = OrfFinder(virus_mode=True)

    # Process vMAG fasta files if a directory is provided
    if vmag_fasta_dir:
        output_vmags = Path(output_dir) / vmag_proteins_subdir
        output_vmags.mkdir(parents=True, exist_ok=True)
        vmag_path = Path(vmag_fasta_dir)
        files = list(vmag_path.glob("*.fna")) + list(vmag_path.glob("*.fasta"))
        
        with load_balancer(files, orf_finder=orf_finder, allow_unordered=True, n_threads=threads) as balancer:
            balancer.submit_to_pool(files, output_vmags, False)

    # Process single contig genome if provided
    if input_fasta:
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        input_path = Path(input_fasta)
        with load_balancer([input_path], orf_finder=orf_finder, allow_unordered=True, n_threads=threads) as balancer:
            balancer.submit_to_pool([input_path], Path(output_dir), False)

        output_faa = Path(output_dir) / (input_path.stem + ".faa")
        if not output_faa.exists():
            logger.error(f"The expected output .faa file {output_faa} was not created.")
            raise RuntimeError(f"Error: The expected output .faa file {output_faa} does not exist.")
        elif output_faa.stat().st_size == 0:
            logger.error(f"The output .faa file {output_faa} is empty. metapyrodigal may have failed.")
            raise RuntimeError(f"Error: The output .faa file {output_faa} is empty.")

        # Only attempt renaming if a valid destination path is provided
        if single_contig_prots is not None:
            logger.debug(f"Renaming {output_faa} to {single_contig_prots}")
            shutil.move(str(output_faa), single_contig_prots)

            if not os.path.exists(single_contig_prots):
                logger.error(f"The output file {single_contig_prots} was not created after renaming.")
                raise RuntimeError(f"Error: The output file {single_contig_prots} does not exist.")
            elif os.path.getsize(single_contig_prots) == 0:
                logger.error(f"The renamed file {single_contig_prots} is empty.")
                raise RuntimeError(f"Error: The renamed file {single_contig_prots} is empty.")
        
def main():
    input_single_contig_genomes = snakemake.params.input_single_contig_genomes
    input_vmag_fastas = snakemake.params.input_vmag_fastas
    wdir = snakemake.params.wdir
    output_dir = snakemake.params.output_dir
    output_single_contig_prots = snakemake.params.single_contig_prots
    vmag_proteins_subdir = snakemake.params.vmag_proteins_subdir
    n_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Metapyrodigal-GV run starting...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    
    if not os.path.exists(wdir):
        os.makedirs(wdir, exist_ok=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    if os.path.exists(input_single_contig_genomes) and os.path.exists(input_vmag_fastas):
        logger.info(f"Running metapyrodigal-GV on single contig genomes and vMAGs with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=input_single_contig_genomes,
            vmag_fasta_dir=input_vmag_fastas,
            output_dir=output_dir,
            vmag_proteins_subdir=vmag_proteins_subdir,
            single_contig_prots=output_single_contig_prots,
            threads=n_cpus
            )
        if not os.path.exists(output_single_contig_prots) or os.path.getsize(output_single_contig_prots) == 0:
            logger.error(f"The output file {output_single_contig_prots} was not created or is empty.")
            raise RuntimeError(f"Error: The output file {output_single_contig_prots} does not exist or is empty.")
        if not Path(vmag_proteins_subdir).exists():
            logger.error(f"The output directory {vmag_proteins_subdir} was not created.")
            raise RuntimeError(f"Error: The output directory {vmag_proteins_subdir} does not exist.")
        folder_files = list(Path(vmag_proteins_subdir).glob("*"))
        if not folder_files:
            logger.error(f"No files were found in the directory {vmag_proteins_subdir}.")
            raise RuntimeError(f"Error: The output directory {vmag_proteins_subdir} is empty.")
        if folder_files[0].stat().st_size == 0:
            logger.error(f"The first file {folder_files[0]} in the directory {vmag_proteins_subdir} is empty.")
            raise RuntimeError(f"Error: The first file in {vmag_proteins_subdir} is empty.")
    elif os.path.exists(input_single_contig_genomes) and not os.path.exists(input_vmag_fastas):
        logger.info(f"Running metapyrodigal-GV on single contig genomes with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=input_single_contig_genomes,
            vmag_fasta_dir=None,
            output_dir=output_dir,
            vmag_proteins_subdir=None,
            single_contig_prots=output_single_contig_prots,
            threads=n_cpus
            )
        if not os.path.exists(output_single_contig_prots) or os.path.getsize(output_single_contig_prots) == 0:
            logger.error(f"The output file {output_single_contig_prots} was not created or is empty.")
            raise RuntimeError(f"Error: The output file {output_single_contig_prots} does not exist or is empty.")
    elif not os.path.exists(input_single_contig_genomes) and os.path.exists(input_vmag_fastas):
        logger.info(f"Running metapyrodigal-GV on vMAGs with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=None,
            vmag_fasta_dir=input_vmag_fastas,
            output_dir=output_dir,
            vmag_proteins_subdir=vmag_proteins_subdir,
            single_contig_prots=None,
            threads=n_cpus
            )
        if not Path(vmag_proteins_subdir).exists():
            logger.error(f"The output directory {vmag_proteins_subdir} was not created.")
            raise RuntimeError(f"Error: The output directory {vmag_proteins_subdir} does not exist.")
        folder_files = list(Path(vmag_proteins_subdir).glob("*"))
        if not folder_files:
            logger.error(f"No files were found in the directory {vmag_proteins_subdir}.")
            raise RuntimeError(f"Error: The output directory {vmag_proteins_subdir} is empty.")
        if folder_files[0].stat().st_size == 0:
            logger.error(f"The first file {folder_files[0]} in the directory {vmag_proteins_subdir} is empty.")
            raise RuntimeError(f"Error: The first file in {vmag_proteins_subdir} is empty.")
    else:
        logger.error("No input genomes or vMAGs provided to run Pyrodigal-GV on.")
        raise FileNotFoundError("Error: No input genomes or vMAGs provided to run Pyrodigal-GV on.")

    logger.info("Metapyrodigal-GV run completed.")

if __name__ == "__main__":
    main()
