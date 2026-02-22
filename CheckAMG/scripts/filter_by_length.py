#!/usr/bin/env python3

import os
import sys
import gzip
import logging
import resource
import multiprocessing as mp
from tqdm import tqdm
from pyfastatools import Parser

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
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode="a"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger()

print("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n            Step 1/11: Filter the input sequences by length             \n========================================================================\n")

class _SimpleHeader:
    __slots__ = ("name",)

    def __init__(self, name):
        self.name = name

class _SimpleRecord:
    __slots__ = ("header", "seq")

    def __init__(self, name, seq):
        self.header = _SimpleHeader(name)
        self.seq = seq

def _iter_fasta_records_gz(path):
    name = None
    seq_chunks = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield _SimpleRecord(name, "".join(seq_chunks))
                header = line[1:].rstrip("\n\r")
                name = header.split(None, 1)[0] if header else ""
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield _SimpleRecord(name, "".join(seq_chunks))

def iter_fasta_records(path):
    if path.endswith(".gz"):
        return _iter_fasta_records_gz(path)
    return Parser(path)

def _find_first_existing(paths):
    for p in paths:
        if p and os.path.exists(p):
            return p
    return None

# Custom fasta writing function that handles very large individual sequences well
def write_fasta_custom(record, handle, small_line=75, big_chunk=100000):
    seq = record.seq
    length = len(seq)
    header = str(record.header.name)
    handle.write(f">{header}\n")
    if length <= 1000000:
        for i in range(0, length, small_line):
            handle.write(seq[i:i + small_line] + "\n")
    else:
        for i in range(0, length, big_chunk):
            handle.write(seq[i:i + big_chunk] + "\n")

def filter_single_record_by_length(args):
    name, seq, min_length = args
    return (name, seq) if len(seq) >= min_length else None

def parallel_processing(input_files, min_length, num_workers, single_contig_path):
    all_records = []
    name_to_source = {}
    genome_names = set()

    for input_file in input_files:
        for record in iter_fasta_records(input_file):
            name = record.header.name
            seq = record.seq
            all_records.append((name, seq))
            name_to_source[name] = (input_file, name, seq)

            if input_file == single_contig_path:
                genome_names.add(name)

        if input_file != single_contig_path:
            genome_names.add(input_file)

    logger.info(f"Number of input sequences: {len(all_records):,} ({len(genome_names):,} genomes)")

    args = [(name, seq, min_length) for name, seq in all_records]

    with mp.Pool(processes=num_workers) as pool:
        filtered = list(
            tqdm(
                pool.imap_unordered(filter_single_record_by_length, args),
                total=len(args),
                desc="Filtering sequences",
                unit="sequence",
            )
        )

    filtered = [item for item in filtered if item is not None]

    grouped = {}
    for name, seq in filtered:
        source_file, _, _ = name_to_source[name]
        grouped.setdefault(source_file, []).append((name, seq))

    return list(grouped.items())

def main():
    input_fasta = snakemake.params.input_single_contigs
    input_bins = snakemake.params.input_bins
    output_folder = snakemake.params.output_parent
    output_single_contigs = snakemake.params.output_single_contigs
    output_bins_folder = snakemake.params.output_bins
    min_length = snakemake.params.min_len
    mem_limit = snakemake.resources.mem
    num_workers = snakemake.threads
    set_memory_limit(mem_limit)

    logger.info("Contig length filtering starting...")

    input_files = []
    if input_fasta:
        input_fasta = _find_first_existing([input_fasta])
        input_files.append(input_fasta)

    if input_bins:
        for p in input_bins.split():
            found = _find_first_existing([p])
            if found:
                input_files.append(found)

    if not input_files:
        raise FileNotFoundError("No input FASTA files found (including .gz).")

    os.makedirs(output_folder, exist_ok=True)

    if len(input_files) > 1:
        os.makedirs(output_bins_folder, exist_ok=True)

    filtered_genomes = parallel_processing(input_files, min_length, num_workers, input_fasta)

    genome_names_filtered = set()

    for input_file, records in filtered_genomes:
        if input_file == input_fasta:
            with open(output_single_contigs, "w", buffering=1024 * 1024) as output_handle:
                for name, seq in records:
                    genome_names_filtered.add(name)
                    write_fasta_custom(_SimpleRecord(name, seq), output_handle)
        else:
            genome_names_filtered.add(input_file)
            out_base = os.path.basename(input_file)
            if out_base.endswith(".gz"):
                out_base = out_base[:-3]
            output_file = os.path.join(output_bins_folder, out_base)
            if output_file.endswith(".fasta"):
                output_file = output_file.replace(".fasta", ".fna")
            elif output_file.endswith(".fa"):
                output_file = output_file.replace(".fa", ".fna")

            with open(output_file, "w", buffering=1024 * 1024) as output_handle:
                for name, seq in records:
                    write_fasta_custom(_SimpleRecord(name, seq), output_handle)

    logger.info(
        f"Number of sequences filtered by length: {sum(len(records) for _, records in filtered_genomes):,} ({len(genome_names_filtered):,} genomes)"
    )
    logger.info("Contig length filtering completed.")

if __name__ == "__main__":
    main()
