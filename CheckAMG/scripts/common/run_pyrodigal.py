#!/usr/bin/env python3

import os
import shutil
from pathlib import Path
from metapyrodigal.load_balancer import SingleFileLoadBalancer, load_balancer
from metapyrodigal.orf_finder import OrfFinder
from typing import Optional

from CheckAMG.scripts.common.fasta_io import collect_nucleotide_inputs
from CheckAMG.scripts.common.validate_fasta import validate_fasta_inputs
from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Predict ORFs and translate into proteins with pyrodigal-GV", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

def run_metapyrodigal(input_fasta: Optional[str],
                      bin_fasta_dir: Optional[str],
                      output_dir: str,
                      bin_proteins_subdir: str,
                      single_contig_prots: str,
                      threads: int) -> None:
    """
    Run metapyrodigal-GV on provided single contig and/or binned fasta files.
    Depending on the provided inputs, this function will execute metapyrodigal on
    binned fasta files located in a directory or on a single contig genome fasta file.
    
    Arguments:
        input_fasta: Path to the single contig genome fasta file.
        bin_fasta_dir: Path to the directory containing binned fasta files.
        output_dir: Directory where output files will be saved.
        bin_proteins_subdir: Subdirectory name for binned protein outputs.
        single_contig_prots: Output file path for single contig protein results.
        threads: Number of CPU threads to leverage.
    """
    orf_finder = OrfFinder(virus_mode=True)

    # Process binned fasta files if a directory is provided
    if bin_fasta_dir:
        output_bins = Path(output_dir) / bin_proteins_subdir
        output_bins.mkdir(parents=True, exist_ok=True)
        bins_path = Path(bin_fasta_dir)
        files = list(bins_path.glob("*.fna")) + list(bins_path.glob("*.fasta"))
        
        with load_balancer(files, orf_finder=orf_finder, allow_unordered=True, n_threads=threads) as balancer:
            balancer.submit_to_pool(files, output_bins, False)

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
    input_single_contigs = snakemake.params.input_single_contigs
    input_bins = snakemake.params.input_bins
    wdir = snakemake.params.wdir
    output_dir = snakemake.params.output_dir
    output_single_contig_prots = snakemake.params.single_contig_prots
    bin_proteins_subdir = snakemake.params.bin_proteins_subdir
    starting_point = snakemake.params.starting_point
    n_cpus = snakemake.threads

    logger.info("Metapyrodigal-GV run starting...")
    
    if not os.path.exists(wdir):
        os.makedirs(wdir, exist_ok=True)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    if starting_point:
        input_files, bin_files, has_single_contig, has_bin_files = collect_nucleotide_inputs(
            input_single_contigs, input_bins, logger,
        )
        assert input_files, "No valid input files found."
        validate_fasta_inputs(input_files, input_type="nucl", logger=logger, label="nucleotide")
    else:
        has_single_contig = input_single_contigs and os.path.exists(input_single_contigs)
        has_bin_files = input_bins and os.path.exists(input_bins)

    if has_single_contig and has_bin_files:
        logger.info(f"Running metapyrodigal-GV on contigs and bins with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=input_single_contigs,
            bin_fasta_dir=input_bins,
            output_dir=output_dir,
            bin_proteins_subdir=bin_proteins_subdir,
            single_contig_prots=output_single_contig_prots,
            threads=n_cpus
            )
        if not os.path.exists(output_single_contig_prots) or os.path.getsize(output_single_contig_prots) == 0:
            logger.error(f"The output file {output_single_contig_prots} was not created or is empty.")
            raise RuntimeError(f"Error: The output file {output_single_contig_prots} does not exist or is empty.")
        if not Path(bin_proteins_subdir).exists():
            logger.error(f"The output directory {bin_proteins_subdir} was not created.")
            raise RuntimeError(f"Error: The output directory {bin_proteins_subdir} does not exist.")
        folder_files = list(Path(bin_proteins_subdir).glob("*"))
        if not folder_files:
            logger.error(f"No files were found in the directory {bin_proteins_subdir}.")
            raise RuntimeError(f"Error: The output directory {bin_proteins_subdir} is empty.")
        if folder_files[0].stat().st_size == 0:
            logger.error(f"The first file {folder_files[0]} in the directory {bin_proteins_subdir} is empty.")
            raise RuntimeError(f"Error: The first file in {bin_proteins_subdir} is empty.")
    elif os.path.exists(input_single_contigs) and not os.path.exists(input_bins):
        logger.info(f"Running metapyrodigal-GV on contigs with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=input_single_contigs,
            bin_fasta_dir=None,
            output_dir=output_dir,
            bin_proteins_subdir=None,
            single_contig_prots=output_single_contig_prots,
            threads=n_cpus
            )
        if not os.path.exists(output_single_contig_prots) or os.path.getsize(output_single_contig_prots) == 0:
            logger.error(f"The output file {output_single_contig_prots} was not created or is empty.")
            raise RuntimeError(f"Error: The output file {output_single_contig_prots} does not exist or is empty.")
    elif not os.path.exists(input_single_contigs) and os.path.exists(input_bins):
        logger.info(f"Running metapyrodigal-GV on bins with {n_cpus} CPUs")
        run_metapyrodigal(
            input_fasta=None,
            bin_fasta_dir=input_bins,
            output_dir=output_dir,
            bin_proteins_subdir=bin_proteins_subdir,
            single_contig_prots=None,
            threads=n_cpus
            )
        if not Path(bin_proteins_subdir).exists():
            logger.error(f"The output directory {bin_proteins_subdir} was not created.")
            raise RuntimeError(f"Error: The output directory {bin_proteins_subdir} does not exist.")
        folder_files = list(Path(bin_proteins_subdir).glob("*"))
        if not folder_files:
            logger.error(f"No files were found in the directory {bin_proteins_subdir}.")
            raise RuntimeError(f"Error: The output directory {bin_proteins_subdir} is empty.")
        if folder_files[0].stat().st_size == 0:
            logger.error(f"The first file {folder_files[0]} in the directory {bin_proteins_subdir} is empty.")
            raise RuntimeError(f"Error: The first file in {bin_proteins_subdir} is empty.")
    else:
        logger.error("No input contigs or bins provided to run Pyrodigal-GV on.")
        raise FileNotFoundError("Error: No input contigs or bins provided to run Pyrodigal-GV on.")

    logger.info("Metapyrodigal-GV run completed.")

if __name__ == "__main__":
    main()
