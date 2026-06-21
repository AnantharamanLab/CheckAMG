#!/usr/bin/env python3

import os
from collections import defaultdict

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.common.fasta_io import (
    collect_protein_inputs,
    collect_pyrodigal_outputs,
    iter_fasta_records
)
from CheckAMG.scripts.common.validate_fasta import validate_fasta_inputs

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Filter the query proteins", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")


def filter_split_and_track(file_path, min_cds, max_cds, max_protein_length):
    contig_cds = defaultdict(list)
    total_orfs_before = 0
    total_orfs_dropped_length = 0

    for record in iter_fasta_records(file_path):
        total_orfs_before += 1
        if max_protein_length and len(record.seq) > max_protein_length:
            total_orfs_dropped_length += 1
            continue
        contig_id = record.header.name.rsplit("_", 1)[0]
        contig_cds[contig_id].append(record)

    total_contigs_before = len(contig_cds)

    kept = {}
    dropped_min = 0
    split_events = 0

    import math

    for contig_id, records in contig_cds.items():
        n = len(records)

        if n < min_cds:
            dropped_min += 1
            continue

        if n <= max_cds:
            kept[contig_id] = records
        else:
            split_events += 1
            n_splits = math.ceil(n / max_cds)
            split_size = math.ceil(n / n_splits)

            for i in range(n_splits):
                start = i * split_size
                end = min((i + 1) * split_size, n)
                chunk = records[start:end]
                if not chunk:
                    continue
                renamed = []
                for r in chunk:
                    name = r.header.name
                    contig_part, protein_num = name.rsplit("_", 1)
                    r.header.name = f"{contig_part}_SPLIT_{i+1}_{protein_num}"
                    renamed.append(r)
                kept[f"{contig_id}_SPLIT_{i+1}"] = renamed

    stats = {
        "contigs_before": total_contigs_before,
        "contigs_after": len(kept),
        "contigs_dropped_min": dropped_min,
        "split_events": split_events,
        "orfs_before": total_orfs_before,
        "orfs_after": sum(len(v) for v in kept.values()),
        "orfs_dropped_length": total_orfs_dropped_length,
    }

    return kept, stats

def _genome_name_from_path(path):
    """Derive a genome/bin name from a FASTA path: drop `.gz` then the FASTA extension."""
    base = os.path.basename(str(path))
    if base.endswith(".gz"):
        base = base[:-3]
    stem, _ = os.path.splitext(base)
    return stem

def filter_all_to_single_output(input_files, output_faa, min_cds, max_cds, max_protein_length,
                                bin_files=None, contig_to_genome_out=None):
    global_stats = defaultdict(int)
    os.makedirs(os.path.dirname(output_faa), exist_ok=True)

    # Files that are bins map every contig they contain to the bin (genome) name;
    # contigs from the single-contig input are their own genome.
    bin_set = set(bin_files or [])

    buffer = []
    chunk_size = 10000

    contig_genome_handle = None
    seen_contigs = set()
    if contig_to_genome_out:
        os.makedirs(os.path.dirname(contig_to_genome_out), exist_ok=True)
        contig_genome_handle = open(contig_to_genome_out, "w", buffering=1024 * 1024)
        contig_genome_handle.write("Contig\tGenome\n")

    try:
        with open(output_faa, "w", buffering=1024 * 1024) as out:
            for file_path in input_files:
                kept, stats = filter_split_and_track(file_path, min_cds, max_cds, max_protein_length)

                for k, v in stats.items():
                    global_stats[k] += v

                is_bin = file_path in bin_set
                bin_genome = _genome_name_from_path(file_path) if is_bin else None

                for key, records in kept.items():
                    if contig_genome_handle is not None:
                        # `kept` keys are the original contig id, or "{contig}_SPLIT_{i}"
                        # for over-large contigs split in filter_split_and_track.
                        original_contig = key.split("_SPLIT_")[0]
                        if original_contig not in seen_contigs:
                            seen_contigs.add(original_contig)
                            genome = bin_genome if is_bin else original_contig
                            contig_genome_handle.write(f"{original_contig}\t{genome}\n")

                    for r in records:
                        desc = r.header.desc
                        header = f">{r.header.name} {desc}\n" if desc else f">{r.header.name}\n"
                        buffer.append(f"{header}{r.seq}\n")
                        if len(buffer) >= chunk_size:
                            out.write("".join(buffer))
                            buffer = []
            if buffer:
                out.write("".join(buffer))
    finally:
        if contig_genome_handle is not None:
            contig_genome_handle.close()

    return global_stats

def main():
    output_filtered_faa = snakemake.params.output_filtered_faa
    min_cds = snakemake.params.min_cds
    max_cds = snakemake.params.max_cds
    max_protein_length = snakemake.params.max_protein_length
    input_single_contig_prots = snakemake.params.input_single_contig_prots
    input_bin_prots = snakemake.params.input_bin_prots
    input_type = snakemake.params.input_type
    contig_to_genome = snakemake.params.contig_to_genome

    logger.info("Protein FASTA filtering starting...")

    if input_type == "nucl":
        input_files, bin_files, has_single, has_bins = collect_pyrodigal_outputs(
            input_single_contig_prots, input_bin_prots, logger,
        )
        assert input_files, "No valid input files found."
    elif input_type == "prot":
        input_files, bin_files, has_single, has_bins = collect_protein_inputs(
            input_single_contig_prots, input_bin_prots, logger,
        )
        assert input_files, "No valid input files found."
        validate_fasta_inputs(input_files, input_type="prot", logger=logger, label="protein")
    else:
        logger.error(f"Invalid input type specified: {input_type}")
        raise ValueError(f"Invalid input type specified: {input_type}")

    stats = filter_all_to_single_output(
        input_files,
        output_filtered_faa,
        min_cds,
        max_cds,
        max_protein_length,
        bin_files=bin_files,
        contig_to_genome_out=contig_to_genome,
    )

    logger.info(f"Contigs before filtering: {stats['contigs_before']:,}")
    logger.info(f"Proteins before filtering: {stats['orfs_before']:,}")
    
    logger.info(f"Contigs split (> {max_cds:,} proteins): {stats['split_events']:,}")
    logger.info(f"Contigs dropped (< {min_cds:,} proteins): {stats['contigs_dropped_min']:,}")
    logger.info(f"Proteins dropped (length > {max_protein_length:,} AAs): {stats['orfs_dropped_length']:,}")

    logger.info(f"Contigs after filtering: {stats['contigs_after']:,}")
    logger.info(f"Proteins after filtering: {stats['orfs_after']:,}")

    logger.info(f"Wrote contig-to-genome mapping to {contig_to_genome}")

    logger.info("Protein FASTA filtering completed.")

if __name__ == "__main__":
    main()