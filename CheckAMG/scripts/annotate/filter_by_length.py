#!/usr/bin/env python3

import os
import multiprocessing as mp
from tqdm import tqdm

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.common.fasta_io import (
    _SimpleRecord,
    collect_nucleotide_inputs,
    iter_fasta_records,
    output_basename_for_bin,
    write_fasta_custom,
)
from CheckAMG.scripts.common.validate_fasta import validate_fasta_inputs

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(
    title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Filter the input sequences by length",
    log_path=log_file,
    append_raw=False,
)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")


def _filter_record_by_length(args):
    name, seq, min_length = args
    return (name, seq) if len(seq) >= min_length else None


def parallel_processing(input_files, single_contig_path, min_length, num_workers):
    all_records = []
    name_to_source = {}
    genome_names = set()

    for input_file in input_files:
        for record in iter_fasta_records(input_file):
            name = record.header.name
            seq = record.seq
            all_records.append((name, seq))
            name_to_source[name] = input_file
            if input_file == single_contig_path:
                genome_names.add(name)
        if input_file != single_contig_path:
            genome_names.add(input_file)

    logger.info(f"Number of input sequences: {len(all_records):,} ({len(genome_names):,} genomes)")

    args = [(name, seq, min_length) for name, seq in all_records]

    with mp.Pool(processes=num_workers) as pool:
        filtered = list(
            tqdm(
                pool.imap_unordered(_filter_record_by_length, args),
                total=len(args),
                desc="Filtering sequences",
                unit="sequence",
            )
        )
    filtered = [item for item in filtered if item is not None]

    grouped = {}
    for name, seq in filtered:
        grouped.setdefault(name_to_source[name], []).append((name, seq))
    return list(grouped.items())


def main():
    input_fasta = snakemake.params.input_single_contigs
    input_bins = snakemake.params.input_bins
    output_folder = snakemake.params.output_parent
    output_single_contigs = snakemake.params.output_single_contigs
    output_bins_folder = snakemake.params.output_bins
    min_length = snakemake.params.min_len
    num_workers = snakemake.threads

    logger.info("Contig length filtering starting...")

    input_files, _bin_files, has_single_contig, has_bin_files = collect_nucleotide_inputs(
        input_fasta, input_bins, logger,
    )
    validate_fasta_inputs(input_files, input_type="nucl", logger=logger, label="nucleotide")

    os.makedirs(output_folder, exist_ok=True)
    if has_bin_files:
        os.makedirs(output_bins_folder, exist_ok=True)

    single_contig_path = input_files[0] if has_single_contig else None
    filtered_genomes = parallel_processing(input_files, single_contig_path, min_length, num_workers)

    genome_names_filtered = set()
    for input_file, records in filtered_genomes:
        if has_single_contig and input_file == single_contig_path:
            with open(output_single_contigs, "w", buffering=1024 * 1024) as output_handle:
                for name, seq in records:
                    genome_names_filtered.add(name)
                    write_fasta_custom(_SimpleRecord(name, seq), output_handle)
        else:
            genome_names_filtered.add(input_file)
            output_file = os.path.join(
                output_bins_folder, output_basename_for_bin(input_file, target_ext=".fna")
            )
            with open(output_file, "w", buffering=1024 * 1024) as output_handle:
                for name, seq in records:
                    write_fasta_custom(_SimpleRecord(name, seq), output_handle)

    logger.info(
        f"Number of sequences filtered by length: "
        f"{sum(len(records) for _, records in filtered_genomes):,} "
        f"({len(genome_names_filtered):,} genomes)"
    )
    logger.info("Contig length filtering completed.")


if __name__ == "__main__":
    main()
