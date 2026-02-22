#!/usr/bin/env python3

import os
import sys
import gzip
import logging
import resource
import shutil
import multiprocessing as mp
from collections import defaultdict
from functools import partial
from pathlib import Path
import re
import shlex
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

print("========================================================================\n        Step 4/11: Filter the input sequences by number of ORFs         \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n        Step 4/11: Filter the input sequences by number of ORFs         \n========================================================================\n")

_PRODIGAL_HEADER_SPLIT_RE = re.compile(r"\s+#\s+")

def _is_prodigal_header_line(line: str) -> bool:
    # Prodigal protein header format is:
    # >seq_id # start # end # strand # ID=...
    if not line.startswith(">"):
        return False
    header = line[1:].strip()
    parts = _PRODIGAL_HEADER_SPLIT_RE.split(header)
    if len(parts) < 5:
        return False

    # parts[1:4] should be ints; strand usually -1/1 (or +/- in some tools, but Prodigal is -1/1)
    try:
        int(parts[1])
        int(parts[2])
        int(parts[3])
    except ValueError:
        return False

    # 5th field usually contains "ID=" somewhere
    if "ID=" not in parts[4]:
        return False

    return True

def validate_prodigal_formatted_faa(faa_file: str, label: str, max_headers: int = 3) -> None:
    """
    Fail fast if a protein FASTA does not have Prodigal-formatted headers.
    """
    if faa_file.endswith(".gz"):
        # Don't rely on Parser supporting gz; validate by inspecting first header line.
        with gzip.open(faa_file, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    if not _is_prodigal_header_line(line):
                        raise RuntimeError(
                            f"{label} FASTA headers are not Prodigal-formatted (first header invalid): {faa_file}"
                        )
                    return
        raise RuntimeError(f"{label} FASTA appears empty or has no headers: {faa_file}")

    # .faa (plain): use pyfastatools prodigal header parser
    try:
        n = 0
        for _ in Parser(faa_file).all_prodigal_headers():
            n += 1
            if n >= max_headers:
                break
        if n == 0:
            raise RuntimeError(f"{label} FASTA appears empty or has no records: {faa_file}")
    except RuntimeError as e:
        logger.error(f"Error parsing {label} .faa file {faa_file}: {e}. Headers are not in Prodigal format.")
        raise RuntimeError(f"Error parsing {label} .faa file {faa_file}: {e}. Headers are not in Prodigal format.")
    
class _SimpleHeader:
    __slots__ = ("name", "desc")

    def __init__(self, name, desc):
        self.name = name
        self.desc = desc

class _SimpleRecord:
    __slots__ = ("header", "seq")

    def __init__(self, name, desc, seq):
        self.header = _SimpleHeader(name, desc)
        self.seq = seq

def _iter_fasta_records_gz(path):
    name = None
    desc = ""
    seq_chunks = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield _SimpleRecord(name, desc, "".join(seq_chunks))
                header = line[1:].rstrip("\n\r")
                if not header:
                    name = ""
                    desc = ""
                else:
                    parts = header.split(None, 1)
                    name = parts[0]
                    desc = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield _SimpleRecord(name, desc, "".join(seq_chunks))

def iter_fasta_records(path):
    if path.endswith(".gz"):
        return _iter_fasta_records_gz(path)
    return Parser(path)

def _is_faa(path):
    return path.endswith(".faa") or path.endswith(".faa.gz")

def count_cds_and_filter_by_contig(file_path, min_cds):
    """
    Groups amino-acid sequences by contig and filters contigs based on a minimum CDS threshold.
    Returns a dictionary of contigs with CDS counts meeting the threshold.
    """
    contig_cds_counts = {}
    for record in iter_fasta_records(file_path):
        contig_id = record.header.name.rsplit("_", 1)[0]
        contig_cds_counts.setdefault(contig_id, []).append(record)

    return {contig_id: cds for contig_id, cds in contig_cds_counts.items() if len(cds) >= min_cds}

def process_nucl_input(input_single_contig_prots, input_bin_prots_dir):
    """
    Processes input files when the input type is 'nucl' (nucleotide-based inputs translated to proteins).
    Checks if single-contig and/or bin protein files exist and organizes input paths.
    """
    has_single_contig = input_single_contig_prots is not None

    has_bin_files = os.path.exists(input_bin_prots_dir) and os.path.isdir(input_bin_prots_dir) and any(
        _is_faa(os.path.join(input_bin_prots_dir, f)) and os.path.isfile(os.path.join(input_bin_prots_dir, f))
        for f in os.listdir(input_bin_prots_dir)
    )

    if not has_single_contig and not has_bin_files:
        logger.error("No valid nucleotide-based input files found.")
        raise FileNotFoundError("No valid nucleotide-based input files found for amino-acid sequence filtering.")

    input_files = []
    bin_files = []

    if has_single_contig:
        if not os.path.isfile(input_single_contig_prots):
            raise FileNotFoundError(f"Single-contig protein FASTA not found: {input_single_contig_prots}")
        input_files.append(input_single_contig_prots)

    if has_bin_files:
        bin_files = sorted(
            os.path.join(input_bin_prots_dir, faa)
            for faa in os.listdir(input_bin_prots_dir)
            if _is_faa(faa) and os.path.isfile(os.path.join(input_bin_prots_dir, faa))
        )
        input_files.extend(bin_files)

    return input_files, bin_files, has_single_contig, (len(bin_files) > 0)

def process_prot_input(input_single_contig_prots, input_bin_prots):
    """
    prot mode:
      - input_single_contig_prots: one FASTA where each contig is a genome
      - input_bin_prots: either a directory OR a list of FASTA paths (list/tuple) OR a whitespace-separated string
        where each FASTA file is a genome (may contain multiple contigs)
    """

    def _clean(p):
        if p is None:
            return None
        if isinstance(p, Path):
            p = str(p)
        else:
            p = str(p)
        return os.path.expandvars(os.path.expanduser(p)).strip()

    def _is_faa_path(p):
        return p.endswith(".faa") or p.endswith(".faa.gz")

    def _list_faa(dir_path):
        return sorted(
            os.path.join(dir_path, f)
            for f in os.listdir(dir_path)
            if _is_faa_path(f) and os.path.isfile(os.path.join(dir_path, f))
        )

    # single contig
    input_single_contig_prots = _clean(input_single_contig_prots)
    has_single_contig = bool(input_single_contig_prots)
    if has_single_contig and not os.path.isfile(input_single_contig_prots):
        raise FileNotFoundError(f"Single-contig protein FASTA not found: {repr(input_single_contig_prots)}")

    # bins
    bin_files = []
    if input_bin_prots:
        # Case A: list/tuple already
        if isinstance(input_bin_prots, (list, tuple)):
            candidates = [_clean(x) for x in input_bin_prots if x]
        else:
            s = _clean(input_bin_prots)

            logger.debug(f"input_bin_prots raw repr: {repr(s)}")

            # Case B: it's a real directory path
            if os.path.isdir(s):
                bin_files = _list_faa(s)
                candidates = []
            # Case C: it's a single file path
            elif os.path.isfile(s) and _is_faa_path(s):
                candidates = [s]
            else:
                # Case D: it's a whitespace-separated string of paths
                # Use shlex so quoted paths still work if they ever appear
                candidates = [os.path.expandvars(os.path.expanduser(x)).strip() for x in shlex.split(s)]

        # Expand any directories inside candidates; keep any .faa/.faa.gz files
        for c in candidates:
            if not c:
                continue
            if os.path.isdir(c):
                bin_files.extend(_list_faa(c))
            elif os.path.isfile(c) and _is_faa_path(c):
                bin_files.append(c)
            else:
                raise FileNotFoundError(f"Bin input is not a directory or .faa/.faa.gz file: {repr(c)}")

    # de-dup while preserving order
    seen = set()
    bin_files = [p for p in bin_files if not (p in seen or seen.add(p))]

    has_bin_files = len(bin_files) > 0

    if not has_single_contig and not has_bin_files:
        raise FileNotFoundError("No valid input protein files found (single-contig and bins both missing).")

    input_files = []
    if has_single_contig:
        input_files.append(input_single_contig_prots)
    if has_bin_files:
        input_files.extend(bin_files)

    logger.debug(f"Detected single-contig FASTA: {has_single_contig}")
    logger.debug(f"Detected bin FASTAs: {len(bin_files):,}")
    if bin_files:
        logger.debug("First 5 bin FASTAs: " + ", ".join(bin_files[:5]))

    return input_files, bin_files, has_single_contig, has_bin_files

def get_totals_before_filtering(input_files, single_contig_file, bin_files):
    """
    Report totals split by source:
      - single-contig FASTA: each contig is a genome
      - bins: each bin file is a genome
    """
    # Single-contig totals
    sc_genomes = 0
    sc_contigs = 0
    sc_orfs = 0

    # Bin totals
    bin_genomes = 0
    bin_contigs = 0
    bin_orfs = 0

    # Single-contig
    if single_contig_file:
        contig_counts = {}
        for record in iter_fasta_records(single_contig_file):
            contig_id = record.header.name.rsplit("_", 1)[0]
            contig_counts.setdefault(contig_id, 0)
            contig_counts[contig_id] += 1

        sc_contigs = len(contig_counts)
        sc_genomes = sc_contigs
        sc_orfs = sum(contig_counts.values())

    # Bins
    for bf in bin_files:
        contig_counts = {}
        for record in iter_fasta_records(bf):
            contig_id = record.header.name.rsplit("_", 1)[0]
            contig_counts.setdefault(contig_id, 0)
            contig_counts[contig_id] += 1

        if contig_counts:
            bin_genomes += 1
            bin_contigs += len(contig_counts)
            bin_orfs += sum(contig_counts.values())

    total_genomes = sc_genomes + bin_genomes
    total_contigs = sc_contigs + bin_contigs
    total_orfs = sc_orfs + bin_orfs

    logger.info(f"Single-contig input: {sc_contigs:,} contigs ({sc_genomes:,} genomes), {sc_orfs:,} ORFs")
    logger.info(f"Bin inputs: {bin_contigs:,} contigs ({bin_genomes:,} genomes), {bin_orfs:,} ORFs")
    logger.info(f"Total before filtering: {total_contigs:,} contigs ({total_genomes:,} genomes), {total_orfs:,} ORFs")

    return total_genomes, total_contigs, total_orfs

def process_bin_file(bin_file, output_folder, min_num_sequences):
    filtered_contigs = count_cds_and_filter_by_contig(bin_file, min_num_sequences)
    out_base = os.path.basename(bin_file)
    if out_base.endswith(".gz"):
        out_base = out_base[:-3]
    output_file = os.path.join(output_folder, out_base)

    if not filtered_contigs:
        if os.path.exists(output_file):
            os.remove(output_file)
        return 0, 0

    buffer = []
    chunk_size = 10000
    total_orfs = 0

    with open(output_file, "w", buffering=1024 * 1024) as output_handle:
        for records in filtered_contigs.values():
            for record in records:
                total_orfs += 1
                desc = record.header.desc
                header = f">{record.header.name} {desc}\n" if desc else f">{record.header.name}\n"
                buffer.append(f"{header}{record.seq}\n")
                if len(buffer) >= chunk_size:
                    output_handle.write("".join(buffer))
                    buffer = []
        if buffer:
            output_handle.write("".join(buffer))

    return len(filtered_contigs), total_orfs

def filter_and_save_bin_proteins(bin_files, bins_output_folder, min_num_sequences):
    if not bin_files:
        return

    os.makedirs(bins_output_folder, exist_ok=True)

    logger.info(f"Processing {len(bin_files)} bins...")

    with mp.Pool(processes=mp.cpu_count()) as pool:
        func = partial(process_bin_file, output_folder=bins_output_folder, min_num_sequences=min_num_sequences)
        results = list(tqdm(pool.imap_unordered(func, bin_files), total=len(bin_files), desc="Filtering bins", unit="file"))

    total_contigs = sum(r[0] for r in results)
    total_orfs = sum(r[1] for r in results)
    logger.info(f"Filtered bins: {total_contigs:,} contigs and {total_orfs:,} ORFs retained")

def extract_and_filter_contig_group(args):
    contig_id, records, min_num_sequences = args
    if len(records) >= min_num_sequences:
        return contig_id, records
    return None

def filter_and_save_single_contig_proteins(input_file, output_folder, min_num_sequences):
    contig_cds = defaultdict(list)
    for record in iter_fasta_records(input_file):
        contig_id = record.header.name.rsplit("_", 1)[0]
        contig_cds[contig_id].append(record)

    args = [
        (contig_id, [(rec.header.name, rec.header.desc, rec.seq) for rec in records], min_num_sequences)
        for contig_id, records in contig_cds.items()
    ]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        filtered_results = list(
            tqdm(
                pool.imap_unordered(extract_and_filter_contig_group, args),
                total=len(args),
                desc="Filtering contigs",
                unit="contig",
            )
        )

    kept = [res for res in filtered_results if res is not None]

    single_contig_output_file = os.path.join(output_folder, "single_contig_proteins.faa")

    if not kept:
        if os.path.exists(single_contig_output_file):
            os.remove(single_contig_output_file)
        logger.info("No contigs met the CDS threshold; no output FASTA written.")
        return

    total_orfs = 0
    buffer = []
    chunk_size = 10000

    with open(single_contig_output_file, "w", buffering=1024 * 1024) as output_handle:
        for _, records_data in kept:
            for name, desc, seq in records_data:
                total_orfs += 1
                header = f">{name} {desc}\n" if desc else f">{name}\n"
                buffer.append(f"{header}{seq}\n")
                if len(buffer) >= chunk_size:
                    output_handle.write("".join(buffer))
                    buffer = []
        if buffer:
            output_handle.write("".join(buffer))

    logger.info(f"Filtered contigs: {len(kept):,} contigs and {total_orfs:,} ORFs retained")

def get_totals_after_filtering(input_files, min_num_sequences, single_contig_file):
    """
    Calculate total contigs and amino-acid sequences after filtering.
    """
    total_filtered_genomes = 0
    total_filtered_contigs = 0
    total_filtered_amino_acid_sequences = 0

    for file_path in input_files:
        filtered_contigs = count_cds_and_filter_by_contig(file_path, min_num_sequences)

        n_contigs = len(filtered_contigs)
        n_orfs = sum(len(cds) for cds in filtered_contigs.values())

        total_filtered_contigs += n_contigs
        total_filtered_amino_acid_sequences += n_orfs

        if single_contig_file is not None and file_path == single_contig_file:
            total_filtered_genomes += n_contigs
        else:
            if n_contigs > 0:
                total_filtered_genomes += 1

    return total_filtered_genomes, total_filtered_contigs, total_filtered_amino_acid_sequences

def main():
    input_type = snakemake.params.input_type
    pyrodigal_dir = snakemake.params.pyrodigal_dir
    output_folder = snakemake.params.output
    bins_output_folder = snakemake.params.output_bins
    min_num_sequences = snakemake.params.min_cds
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Amino-acid sequence filtering starting...")

    input_files, bin_files, has_single_contig, has_bin_files = [], [], False, False

    if input_type == "nucl":
        input_single_contig_prots = snakemake.params.input_single_contig_prots
        logger.debug(f"Input single-contig proteins: {input_single_contig_prots}")
        input_bin_prots = snakemake.params.input_bin_prots
        logger.debug(f"Input bin proteins dir: {input_bin_prots}")
        input_files, bin_files, has_single_contig, has_bin_files = process_nucl_input(input_single_contig_prots, input_bin_prots)

    elif input_type == "prot":
        input_single_contig_prots = snakemake.params.input_single_contig_prots
        logger.debug(f"Input single-contig proteins: {input_single_contig_prots}")
        input_bin_prots = snakemake.params.input_bin_prots
        logger.debug(f"Input bin proteins: {input_bin_prots}")
        input_files, bin_files, has_single_contig, has_bin_files = process_prot_input(input_single_contig_prots, input_bin_prots)
        if has_single_contig:
            validate_prodigal_formatted_faa(input_single_contig_prots, label="single-contig")
        if has_bin_files:
            for bf in bin_files:
                validate_prodigal_formatted_faa(bf, label="bin")
    else:
        logger.error(f"Invalid input type specified: {input_type}")
        raise ValueError(f"Invalid input type specified: {input_type}")

    single_contig_file = input_files[0] if has_single_contig else None

    total_genomes, total_contigs, total_amino_acid_sequences = get_totals_before_filtering(
        input_files, single_contig_file, bin_files
    )
    logger.info(f"Total contigs before filtering: {total_contigs:,} ({total_genomes:,} genomes)")
    logger.info(f"Total amino-acid sequences before filtering: {total_amino_acid_sequences:,}")

    logger.debug(f"Output folder: {output_folder}")
    os.makedirs(output_folder, exist_ok=True)

    if has_bin_files:
        filter_and_save_bin_proteins(bin_files, bins_output_folder, min_num_sequences)

    if has_single_contig:
        filter_and_save_single_contig_proteins(single_contig_file, output_folder, min_num_sequences)

    total_filtered_genomes, total_filtered_contigs, total_filtered_amino_acid_sequences = get_totals_after_filtering(
        input_files, min_num_sequences, single_contig_file
    )

    logger.info(f"Total contigs after filtering: {total_filtered_contigs:,} ({total_filtered_genomes:,} genomes)")
    logger.info(f"Total amino-acid sequences after filtering: {total_filtered_amino_acid_sequences:,}")

    if input_type == "nucl":
        logger.debug(f"Removing original pyrodigal-gv directory: {pyrodigal_dir}")
        shutil.rmtree(pyrodigal_dir)

    logger.info("Amino-acid sequence filtering completed.")

if __name__ == "__main__":
    main()