#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
import multiprocessing as mp
from tqdm import tqdm
from pyfastatools import Parser, Record, Header, write_fasta

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

print("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================\n")

# Function to filter genomes by length
def filter_single_record_by_length(args):
    name, seq, min_length = args
    return (name, seq) if len(seq) >= min_length else None

# Parallel processing function
def parallel_processing(input_files, min_length, num_workers):
    all_records = []
    seq_to_source = {}
    genome_names = set()

    # Preprocess all records
    for input_file in input_files:
        for record in Parser(input_file):
            name = record.header.name
            seq = record.seq
            all_records.append((name, seq))
            seq_to_source[name] = (input_file, record)
            if input_file == snakemake.params.input_single_contig_genomes:
                genome_names.add(name)
        if input_file != snakemake.params.input_single_contig_genomes:
            genome_names.add(input_file)

    logger.info(f"Number of input sequences: {len(all_records):,} ({len(genome_names):,} genomes)")

    # Prepare args for multiprocessing
    args = [(name, seq, min_length) for name, seq in all_records]

    # Filter
    with mp.Pool(processes=num_workers) as pool:
        filtered = list(
            tqdm(
                pool.imap_unordered(filter_single_record_by_length, args),
                total=len(args),
                desc="Filtering sequences",
                unit="sequence"
            )
        )

    # Remove filtered-out entries
    filtered = [item for item in filtered if item is not None]

    # Rebuild full Record objects and group by file
    grouped = {}
    for name, _ in filtered:
        source_file, record = seq_to_source[name]
        grouped.setdefault(source_file, []).append(record)

    return list(grouped.items())

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
    genome_names_filtered = set()
    for input_file, records in filtered_genomes:
        if input_file == input_fasta:
            with open(single_contig_output_file, "w") as output_handle:
                for record in records:
                    genome_names_filtered.add(record.header.name)
                    write_fasta(record, output_handle)
        else:
            genome_names_filtered.add(input_file)
            output_file = os.path.join(vmag_output_folder, os.path.basename(input_file))
            if output_file.endswith(".fasta"):
                output_file = output_file.replace(".fasta", ".fna")
            elif output_file.endswith(".fa"):
                output_file = output_file.replace(".fa", ".fna")                    
            with open(output_file, "w") as output_handle:
                for record in records:
                    write_fasta(record, output_handle)

    logger.info(f"Number of sequences filtered by length: {sum(len(records) for _, records in filtered_genomes):,} ({len(genome_names_filtered):,} genomes)")
    logger.info("Genome length filtering completed.")

if __name__ == "__main__":
    main()
