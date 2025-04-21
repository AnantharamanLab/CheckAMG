#!/usr/bin/env python3

import os
import sys
import resource
import platform
from Bio.SeqRecord import SeqRecord
from pyfastatools import Parser, write_fasta
from pathlib import Path
import logging
from multiprocessing import Pool
import polars as pl
import pymuscle5

os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)

def set_memory_limit(limit_in_gb):
    """
    Sets a memory limit (in GB) for this process on Linux (per-process).
    Does not cap combined usage if multiple processes run in parallel.
    """
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

print("========================================================================\n     Step 20/22: Align the filtered protein clusters with PyMuscle5     \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n     Step 20/22: Align the filtered protein clusters with PyMuscle5     \n========================================================================\n")

def load_filtered_clusters(cluster_file):
    """
    Loads the filtered protein clusters file (TSV with no header),
    returning a DataFrame with columns: protein_cluster_rep, numeric_id
    """
    return pl.read_csv(cluster_file, separator='\t', has_header=False).rename({
        "column_1": "protein_cluster_rep",
        "column_2": "numeric_id"
    })

def load_dblookup(dblookup_file):
    """
    Loads the DB lookup file (TSV with no header) mapping numeric_id -> protein_name
    """
    return pl.read_csv(dblookup_file, separator='\t', has_header=False).rename({
        "column_1": "numeric_id",
        "column_2": "protein_name"
    })

def get_clusters_to_proteins(clusters_df, dblookup_df):
    """
    Joins clusters with the DB lookup, grouping by cluster.
    Produces a DataFrame with columns: protein_cluster_rep, proteins (list of protein_names).
    """
    merged_df = clusters_df.join(dblookup_df, on="numeric_id")
    return merged_df.group_by("protein_cluster_rep").agg([
        pl.col("protein_name").alias("proteins")
    ])

def extract_sequences_from_fasta(fasta_path, protein_ids):
    """
    Loads sequences from the given FASTA into memory,
    but only for the protein_ids of interest.
    Returns a dict { protein_id: record }
    """
    sequences = {}
    for record in Parser(fasta_path):
        if record.header.name in protein_ids:
            # Remove any '*' from the Seq at load time
            record.seq = record.seq.replace("*", "")
            sequences[record.header.name] = record
    return sequences

def remove_cluster_seqs_from_dict(protein_ids, sequences_dict):
    """
    Removes sequences in 'protein_ids' from sequences_dict to free memory.
    """
    for pid in protein_ids:
        if pid in sequences_dict:
            del sequences_dict[pid]

def write_sequences_to_file(seq_list, path_out):
    """
    Writes a list of records to path_out in FASTA format.
    """
    with path_out.open('w') as outfile:
        for seq in seq_list:
            write_fasta(seq, outfile)

def align_cluster_serial(cluster_name, protein_ids, local_seq_dict, output_dir, acc_prefix, index, aligner_threads):
    """
    Performs a single-cluster alignment serially in the main process.
    If the cluster is a singleton or all sequences are identical, it writes them out without running PyMuscle5.
    Otherwise, it calls PyMuscle5.Aligner with the aligner_threads argument.
    Returns the accession basename, or None if no sequences found.
    """
    accession = f"{acc_prefix}_{index}"
    aligned_dir = Path(output_dir) / "aligned"
    aligned_dir.mkdir(parents=True, exist_ok=True)
    aligned_fasta_path = aligned_dir / f"{accession}.aligned.fasta"

    cluster_records = [local_seq_dict[pid] for pid in protein_ids if pid in local_seq_dict]
    if not cluster_records:
        logger.warning(f"No sequences found for cluster {cluster_name}. Skipping.")
        return None

    if len(cluster_records) == 1:
        write_sequences_to_file(cluster_records, aligned_fasta_path)
        return aligned_fasta_path.stem.split(".")[0]

    seq_strings = {str(rec.seq) for rec in cluster_records}
    if len(seq_strings) == 1:
        write_sequences_to_file(cluster_records, aligned_fasta_path)
        return aligned_fasta_path.stem.split(".")[0]

    sequences_for_alignment = [
        pymuscle5.Sequence(
            rec.header.name.encode(),
            rec.seq.encode()
        )
        for rec in cluster_records
        ]

    # Pass the desired number of threads to PyMuscle5's Aligner
    aligner = pymuscle5.Aligner(threads=aligner_threads)
    msa = aligner.align(sequences_for_alignment)

    with aligned_fasta_path.open('w') as out:
        for seq in msa.sequences:
            out.write(f">{seq.name.decode()}\n{seq.sequence.decode()}\n")
    return aligned_fasta_path.stem.split(".")[0]

def serial_align_clusters(truly_multi_clusters, sequences_dict, output_dir, acc_prefix, aligner_threads):
    """
    Processes each cluster sequentially (serially) using the shared sequences_dict.
    Returns a dictionary mapping cluster_name -> accession_basename.
    """
    prot_clust_to_accession = {}
    for cluster_name, prot_list, idx in truly_multi_clusters:
        basename = align_cluster_serial(cluster_name, prot_list, sequences_dict, output_dir, acc_prefix, idx, aligner_threads)
        if basename is not None:
            prot_clust_to_accession[cluster_name] = basename
        remove_cluster_seqs_from_dict(prot_list, sequences_dict)
    return prot_clust_to_accession

def main():
    cluster_file = snakemake.params.cluster_file
    dblookup_file = snakemake.params.dblookup
    fasta = snakemake.params.all_filtered_prots
    acc_prefix = snakemake.params.acc_prefix
    prot_clust_to_accession_path = snakemake.params.prot_clust_to_accession
    output_dir = snakemake.params.wdir
    threads = snakemake.threads  # Pass this to PyMuscle5's internal thread pool for each alignment.
    mem_limit = snakemake.resources.mem

    logger.info("Starting protein cluster alignment with PyMuscle5...")
    set_memory_limit(mem_limit)
    logger.debug(f"Memory limit set to {mem_limit:,} GB.")
    
    clusters_df = load_filtered_clusters(cluster_file)
    dblookup_df = load_dblookup(dblookup_file)
    num_clusters = clusters_df['protein_cluster_rep'].n_unique()
    num_proteins = dblookup_df.shape[0]
    logger.info(f"There are {num_clusters:,} clusters, covering {num_proteins:,} proteins in total.")

    cluster_info = get_clusters_to_proteins(clusters_df, dblookup_df)
    protein_ids = set()
    for row in cluster_info.iter_rows(named=True):
        protein_ids.update(row['proteins'])

    sequences_dict = extract_sequences_from_fasta(fasta, protein_ids)
    total_loaded = len(sequences_dict)
    logger.debug(f"Loaded {total_loaded:,} protein sequences in memory initially.")

    Path(output_dir, "aligned").mkdir(parents=True, exist_ok=True)

    # Phase 1: Filter out singletons or identical clusters.
    truly_multi_clusters = []
    prot_clust_to_accession = {}
    singletons_count = 0
    identical_count = 0

    cluster_index = 0
    for row in cluster_info.iter_rows(named=True):
        cluster_index += 1
        cluster_name = row['protein_cluster_rep']
        prot_list = row['proteins']
        cluster_records = [sequences_dict[pid] for pid in prot_list if pid in sequences_dict]
        if not cluster_records:
            logger.warning(f"No sequences found for cluster {cluster_name}. Skipping.")
            continue

        if len(cluster_records) == 1:
            out_path = Path(output_dir) / "aligned" / f"{acc_prefix}_{cluster_index}.aligned.fasta"
            write_sequences_to_file(cluster_records, out_path)
            prot_clust_to_accession[cluster_name] = out_path.stem.split(".")[0]
            remove_cluster_seqs_from_dict(prot_list, sequences_dict)
            singletons_count += 1
        else:
            seq_strings = {str(rec.seq) for rec in cluster_records}
            if len(seq_strings) == 1:
                out_path = Path(output_dir) / "aligned" / f"{acc_prefix}_{cluster_index}.aligned.fasta"
                write_sequences_to_file(cluster_records, out_path)
                prot_clust_to_accession[cluster_name] = out_path.stem.split(".")[0]
                remove_cluster_seqs_from_dict(prot_list, sequences_dict)
                identical_count += 1
            else:
                truly_multi_clusters.append((cluster_name, prot_list, cluster_index))

    filtered_loaded = len(sequences_dict)
    logger.debug(f"Filtered out {singletons_count:,} singleton clusters and {identical_count:,} identical clusters.")
    logger.debug(f"Remaining clusters to align: {len(truly_multi_clusters):,}.")
    logger.debug(f"Memory dictionary reduced from {total_loaded:,} to {filtered_loaded:,} sequences after filtering.")

    # Phase 2: Align the truly multi clusters serially (using PyMuscle5's internal threading)
    if truly_multi_clusters:
        partial_res = serial_align_clusters(truly_multi_clusters, sequences_dict, output_dir, acc_prefix, threads)
        prot_clust_to_accession.update(partial_res)

    # Phase 3: Check for clusters with empty or no alignment output
    unaligned_clusters = []
    for row in cluster_info.iter_rows(named=True):
        c = row['protein_cluster_rep']
        if c not in prot_clust_to_accession:
            unaligned_clusters.append(c)
        else:
            basename = prot_clust_to_accession[c]
            path_check = Path(output_dir) / "aligned" / f"{basename}.aligned.fasta"
            if not path_check.exists() or path_check.stat().st_size == 0:
                unaligned_clusters.append(c)

    if unaligned_clusters:
        logger.warning(f"The following clusters were not aligned or had empty output: {unaligned_clusters}")
        None
    else:
        logger.info("All clusters processed successfully.")

    mapping_list = [(key, value) for key, value in prot_clust_to_accession.items()]
    df = pl.DataFrame(mapping_list, schema=['protein_cluster_rep', 'accession'], orient="row")
    df.write_csv(prot_clust_to_accession_path, separator='\t')
    logger.info("Protein cluster alignment with PyMuscle5 completed.")

if __name__ == "__main__":
    main()