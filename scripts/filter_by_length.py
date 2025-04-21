#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
import multiprocessing as mp
from pyfastatools import Parser
from textwrap import wrap

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

if snakemake.params.build_or_annotate =="build":
    print("========================================================================\n            Step 1/22: Filter the input sequences by length             \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n            Step 1/22: Filter the input sequences by length             \n========================================================================\n")
elif snakemake.params.build_or_annotate == "annotate":
    print("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================\n")

# Function to filter genomes by length
def filter_genomes_by_length(genome_records, min_length):
    """Filter genomes by length."""
    filtered_records = [(header, seq) for (header, seq) in genome_records if len(seq) >= min_length]
    return filtered_records

# Parallel processing function
def parallel_processing(input_files, min_length, num_workers):
    # Read genome records from all input files
    genome_records = []
    for input_file in input_files:
        records = [(record.header.name, record.seq) for record in Parser(input_file)]
        genome_records.append((input_file, records))
    logger.info(f"Number of input sequences: {sum(len(records) for _, records in genome_records):,}")
    
    # Process each file in parallel
    with mp.Pool(processes=num_workers) as pool:
        results = pool.starmap(
            filter_genomes_by_length,
            [(records, min_length) for _, records in genome_records]
        )

    return list(zip([file for file, _ in genome_records], results))

def main():
    input_fasta = snakemake.params.input_single_contig_genomes
    input_vmag_fastas = snakemake.params.input_vmag_fastas
    output_folder = snakemake.params.output
    min_length = snakemake.params.min_len
    mem_limit = snakemake.resources.mem
    num_workers = snakemake.threads
    set_memory_limit(mem_limit)

    logger.info("Genome length filtering starting...")

    if input_fasta:
        input_files = [input_fasta] + input_vmag_fastas.split()
    else:
        input_files = input_vmag_fastas.split()
    
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    if len(input_files) > 1:
        vmag_output_folder = os.path.join(output_folder, "vMAG_fna")
        if not os.path.exists(vmag_output_folder):
            os.makedirs(vmag_output_folder)

    filtered_genomes = parallel_processing(input_files, min_length, num_workers)

    single_contig_output_file = os.path.join(output_folder, "single_contig_genomes.fna")
    for input_file, records in filtered_genomes:
        if input_file == input_fasta:
            with open(single_contig_output_file, "w") as output_handle:
                for record in records:
                    output_handle.write(f">{record[0]}\n")
                    for line in wrap(record[1], width=75):
                        output_handle.write(f"{line}\n")
        else:
            if len(input_files) > 1:
                output_file = os.path.join(vmag_output_folder, os.path.basename(input_file))
                if output_file.endswith(".fasta"):
                    output_file = output_file.replace(".fasta", ".fna")
                elif output_file.endswith(".fa"):
                    output_file = output_file.replace(".fa", ".fna")                    
                with open(output_file, "w") as output_handle:
                    for record in records:
                        output_handle.write(f">{record[0]}\n")
                        for line in wrap(record[1], width=75):
                            output_handle.write(f"{line}\n")

    logger.info(f"Number of sequences filtered by length: {sum(len(records) for _, records in filtered_genomes):,}")
    logger.info("Genome length filtering completed.")

if __name__ == "__main__":
    main()
