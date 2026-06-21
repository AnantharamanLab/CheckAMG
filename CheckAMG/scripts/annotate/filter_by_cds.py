#!/usr/bin/env python3

import os
import shutil
import multiprocessing as mp
from collections import defaultdict
from functools import partial
from tqdm import tqdm

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.common.fasta_io import (
    collect_protein_inputs,
    collect_pyrodigal_outputs,
    iter_fasta_records,
    output_basename_for_bin,
)
from CheckAMG.scripts.common.validate_fasta import (
    assert_unique_fasta_names,
    validate_fasta_inputs,
)

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(
    title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Filter the input sequences by number of ORFs",
    log_path=log_file,
    append_raw=False,
)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

CHUNK_SIZE = 10000


def _contig_id_from_protein(name):
    return name.rsplit("_", 1)[0]


def _format_record(name, desc, seq):
    return (f">{name} {desc}\n" if desc else f">{name}\n") + seq + "\n"


def count_cds_and_filter_by_contig(file_path, min_cds):
    """Group amino-acid records by source contig and keep only contigs with >= min_cds ORFs."""
    contig_groups = defaultdict(list)
    for record in iter_fasta_records(file_path):
        contig_groups[_contig_id_from_protein(record.header.name)].append(record)
    return {cid: cds for cid, cds in contig_groups.items() if len(cds) >= min_cds}


def get_totals_before_filtering(single_contig_file, bin_files):
    sc_genomes = sc_contigs = sc_orfs = 0
    bin_genomes = bin_contigs = bin_orfs = 0

    if single_contig_file:
        contig_counts = defaultdict(int)
        for record in iter_fasta_records(single_contig_file):
            contig_counts[_contig_id_from_protein(record.header.name)] += 1
        sc_contigs = len(contig_counts)
        sc_genomes = sc_contigs
        sc_orfs = sum(contig_counts.values())

    for bf in bin_files:
        contig_counts = defaultdict(int)
        for record in iter_fasta_records(bf):
            contig_counts[_contig_id_from_protein(record.header.name)] += 1
        if contig_counts:
            bin_genomes += 1
            bin_contigs += len(contig_counts)
            bin_orfs += sum(contig_counts.values())

    total_genomes = sc_genomes + bin_genomes
    total_contigs = sc_contigs + bin_contigs
    total_orfs = sc_orfs + bin_orfs

    logger.info(f"Single-contig input: {sc_contigs:,} contigs ({sc_genomes:,} genomes), {sc_orfs:,} ORFs")
    logger.info(f"Bin inputs: {bin_contigs:,} contigs ({bin_genomes:,} genomes), {bin_orfs:,} ORFs")
    logger.info(
        f"Total before filtering: {total_contigs:,} contigs "
        f"({total_genomes:,} genomes), {total_orfs:,} ORFs"
    )
    return total_genomes, total_contigs, total_orfs


def process_bin_file(bin_file, output_folder, min_num_sequences):
    filtered_contigs = count_cds_and_filter_by_contig(bin_file, min_num_sequences)
    output_file = os.path.join(output_folder, output_basename_for_bin(bin_file))

    if not filtered_contigs:
        if os.path.exists(output_file):
            os.remove(output_file)
        return 0, 0

    total_orfs = 0
    buffer = []
    with open(output_file, "w", buffering=1024 * 1024) as output_handle:
        for records in filtered_contigs.values():
            for record in records:
                total_orfs += 1
                desc = getattr(record.header, "desc", "") or ""
                buffer.append(_format_record(record.header.name, desc, record.seq))
                if len(buffer) >= CHUNK_SIZE:
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
    func = partial(process_bin_file, output_folder=bins_output_folder, min_num_sequences=min_num_sequences)
    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(
            tqdm(
                pool.imap_unordered(func, bin_files),
                total=len(bin_files),
                desc="Filtering bins",
                unit="file",
            )
        )

    total_contigs = sum(r[0] for r in results)
    total_orfs = sum(r[1] for r in results)
    logger.info(f"Filtered bins: {total_contigs:,} contigs and {total_orfs:,} ORFs retained")


def _extract_and_filter_contig_group(args):
    contig_id, records, min_num_sequences = args
    return (contig_id, records) if len(records) >= min_num_sequences else None


def filter_and_save_single_contig_proteins(input_file, output_folder, min_num_sequences):
    contig_cds = defaultdict(list)
    for record in iter_fasta_records(input_file):
        contig_cds[_contig_id_from_protein(record.header.name)].append(record)

    args = [
        (
            cid,
            [(rec.header.name, getattr(rec.header, "desc", "") or "", rec.seq) for rec in records],
            min_num_sequences,
        )
        for cid, records in contig_cds.items()
    ]

    with mp.Pool(processes=mp.cpu_count()) as pool:
        filtered_results = list(
            tqdm(
                pool.imap_unordered(_extract_and_filter_contig_group, args),
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
    with open(single_contig_output_file, "w", buffering=1024 * 1024) as output_handle:
        for _, records_data in kept:
            for name, desc, seq in records_data:
                total_orfs += 1
                buffer.append(_format_record(name, desc, seq))
                if len(buffer) >= CHUNK_SIZE:
                    output_handle.write("".join(buffer))
                    buffer = []
        if buffer:
            output_handle.write("".join(buffer))

    logger.info(f"Filtered contigs: {len(kept):,} contigs and {total_orfs:,} ORFs retained")


def get_totals_after_filtering(input_files, single_contig_file, min_num_sequences):
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
        elif n_contigs > 0:
            total_filtered_genomes += 1
    return total_filtered_genomes, total_filtered_contigs, total_filtered_amino_acid_sequences


def _resolve_inputs():
    """
    Dispatch input collection & validation by `input_type`:
      - 'nucl': inputs come from upstream pyrodigal-gv (workflow-managed paths;
                may be absent if user only supplied one of contigs/bins).
                Headers are guaranteed Prodigal-formatted by construction; we
                still verify uniqueness as a safety net.
      - 'prot': inputs come directly from the user. Run the full validator
                (Prodigal format + no whitespace in seq_id + uniqueness).
    """
    input_type = snakemake.params.input_type
    input_single_contig_prots = snakemake.params.input_single_contig_prots
    input_bin_prots = snakemake.params.input_bin_prots
    logger.debug(f"Input type: {input_type}")
    logger.debug(f"Input single-contig proteins: {input_single_contig_prots}")
    logger.debug(f"Input bin proteins: {input_bin_prots}")

    if input_type == "nucl":
        input_files, bin_files, has_single, has_bins = collect_pyrodigal_outputs(
            input_single_contig_prots, input_bin_prots, logger,
        )
        for f in input_files:
            logger.debug(f"Checking for duplicate FASTA headers: {f}")
            assert_unique_fasta_names(f, logger)
    elif input_type == "prot":
        input_files, bin_files, has_single, has_bins = collect_protein_inputs(
            input_single_contig_prots, input_bin_prots, logger,
        )
        validate_fasta_inputs(input_files, input_type="prot", logger=logger, label="protein")
    else:
        logger.error(f"Invalid input type specified: {input_type}")
        raise ValueError(f"Invalid input type specified: {input_type}")

    return input_type, input_files, bin_files, has_single, has_bins


def main():
    output_folder = snakemake.params.output
    bins_output_folder = snakemake.params.output_bins
    min_num_sequences = snakemake.params.min_cds
    pyrodigal_dir = snakemake.params.pyrodigal_dir

    logger.info("Amino-acid sequence filtering starting...")

    input_type, input_files, bin_files, has_single_contig, has_bin_files = _resolve_inputs()
    single_contig_file = input_files[0] if has_single_contig else None

    total_genomes, total_contigs, total_orfs = get_totals_before_filtering(single_contig_file, bin_files)
    logger.info(f"Total contigs before filtering: {total_contigs:,} ({total_genomes:,} genomes)")
    logger.info(f"Total amino-acid sequences before filtering: {total_orfs:,}")

    logger.debug(f"Output folder: {output_folder}")
    os.makedirs(output_folder, exist_ok=True)

    if has_bin_files:
        filter_and_save_bin_proteins(bin_files, bins_output_folder, min_num_sequences)
    if has_single_contig:
        filter_and_save_single_contig_proteins(single_contig_file, output_folder, min_num_sequences)

    f_genomes, f_contigs, f_orfs = get_totals_after_filtering(
        input_files, single_contig_file, min_num_sequences,
    )
    logger.info(f"Total contigs after filtering: {f_contigs:,} ({f_genomes:,} genomes)")
    logger.info(f"Total amino-acid sequences after filtering: {f_orfs:,}")

    if input_type == "nucl":
        logger.debug(f"Removing original pyrodigal-gv directory: {pyrodigal_dir}")
        shutil.rmtree(pyrodigal_dir)

    logger.info("Amino-acid sequence filtering completed.")


if __name__ == "__main__":
    main()
