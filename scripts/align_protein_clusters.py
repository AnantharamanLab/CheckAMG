#!/usr/bin/env python3

import os
import sys
import math
import resource
import platform
from pyfastatools import Parser, write_fasta
from pathlib import Path
import logging
from multiprocessing import Pool
import time
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import pymuscle5
import gc; gc.collect()

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


# Module‐level globals for worker processes
sequences_dict = None
_acc_prefix = None
_output_dir = None
_aligner_threads = None
    
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
    merged_df = merged_df.group_by("protein_cluster_rep").agg([
        pl.col("protein_name").alias("proteins")
    ]).sort('protein_cluster_rep').with_row_index('index',offset=1) # Stable index: sort by cluster name
    return merged_df

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

def init_worker(seqs, acc_prefix, output_dir, aligner_threads):
    """
    Initialize globals for each worker process.

    seqs: dict of protein_id -> sequences (str) loaded from FASTA
    acc_prefix: prefix for output filenames
    output_dir: base directory where "aligned/" subfolder lives
    aligner_threads: number of threads passed to pymuscle5.Aligner (per cluster)
    """
    global sequences_dict, _acc_prefix, _output_dir, _aligner_threads, _aligner
    sequences_dict = seqs
    _acc_prefix = acc_prefix
    _output_dir = output_dir
    _aligner_threads = aligner_threads
    _aligner = pymuscle5.Aligner(threads=_aligner_threads)

def align_cluster_worker(args):
    """
    For a single cluster, filter or align:

    1. Fetch sequences from sequences_dict.
    2. If no sequences, skip.
    3. If exactly one sequence or all sequences identical:
       - Write FASTA directly to aligned/.
       - Return accession basename.
    4. Otherwise:
       - Run PyMuscle5 alignment with one thread.
       - Write MSA to aligned/.
       - Return accession basename.

    Returns (cluster_name, protein_list, accession_or_None).
    """
    cluster_name, protein_list, idx = args
    accession = f"{_acc_prefix}_{idx}"
    aligned_dir = Path(_output_dir) / "aligned"
    out_path = aligned_dir / f"{accession}.aligned.fasta"

    # Gather records
    recs = [ sequences_dict[p] for p in protein_list if p in sequences_dict ]
    if not recs:
        logger.warning(f"No sequences found for cluster {cluster_name}. Skipping.")
        return (cluster_name, protein_list, None)

    # Singleton or identical
    seqs_set = { str(r.seq) for r in recs }
    if len(recs) == 1 or len(seqs_set) == 1:
        with out_path.open('w') as fh:
            for rec in recs:
                write_fasta(rec, fh)
        return (cluster_name, protein_list, accession)

    # Multi-sequence alignment
    seq_objs = [ pymuscle5.Sequence(r.header.name.encode(), r.seq.encode()) for r in recs ]
    msa = _aligner.align(seq_objs)
    with out_path.open('w') as fh:
        for seq in msa.sequences:
            fh.write(f">{seq.name.decode()}\n{seq.sequence.decode()}\n")
    return (cluster_name, protein_list, accession)

def super_worker(group):
    """
    Aligns a list of clusters in one call to reduce IPC overhead.
    Returns a list of (cluster_name, protein_list, accession) for those succeeded.
    """
    results = []
    for args in group:
        name, prots, idx = args
        _, _, acc = align_cluster_worker((name, prots, idx))
        if acc:
            results.append((name, prots, acc))
    return results

def remove_cluster_seqs_from_dict(protein_ids, seq_dict):
    for pid in protein_ids:
        seq_dict.pop(pid, None)
        
def write_sequences_to_file(seq_list, path_out):
    """
    Writes a list of records to path_out in FASTA format.
    """
    with path_out.open('w') as outfile:
        for seq in seq_list:
            write_fasta(seq, outfile)

def main():
    cluster_file = snakemake.params.cluster_file
    dblookup_file = snakemake.params.dblookup
    fasta = snakemake.params.all_filtered_prots
    acc_prefix = snakemake.params.acc_prefix
    prot_clust_to_accession_path = snakemake.params.prot_clust_to_accession
    output_dir = snakemake.params.wdir
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem

    logger.info("Starting protein cluster alignment with PyMuscle5...")
    set_memory_limit(mem_limit)
    logger.debug(f"Memory limit set to {mem_limit:,} GB.")
    
    # Load data
    clusters_df = load_filtered_clusters(cluster_file)
    dblookup_df = load_dblookup(dblookup_file)
    num_clusters = clusters_df['protein_cluster_rep'].n_unique()
    num_proteins = dblookup_df.shape[0]
    logger.info(f"There are {num_clusters:,} clusters, covering {num_proteins:,} proteins in total.")

    cluster_info = get_clusters_to_proteins(clusters_df, dblookup_df)

    # Gather all protein IDs
    protein_ids = set()
    for prots in cluster_info['proteins'].to_list():
        protein_ids.update(prots)

    # Load sequences
    sequences_dict_local = extract_sequences_from_fasta(fasta, protein_ids)
    total_loaded = len(sequences_dict_local)
    logger.debug(f"Loaded {total_loaded:,} protein sequences in memory initially.")
    Path(output_dir, 'aligned').mkdir(parents=True, exist_ok=True)

    # Existing alignments on disk
    aligned_dir = Path(output_dir) / 'aligned'
    existing = {
        f.stem.replace(".aligned.fasta", "")
        for f in aligned_dir.glob(f"{acc_prefix}_*.aligned.fasta")
        if f.stat().st_size > 0
    }
    logger.info(f"Found {len(existing):,} existing alignments to skip.")

    prot_clust_to_accession = {}
    singletons = identical = 0

    for row in cluster_info.iter_rows(named=True):
        idx  = row['index']
        name = row['protein_cluster_rep']
        prots = row['proteins']
        acc  = f"{acc_prefix}_{idx}"
        if acc in existing:
            continue

        recs = [sequences_dict_local[p] for p in prots if p in sequences_dict_local]
        if len(recs) == 1:
            singletons += 1
        elif len(recs) > 1 and len({str(r.seq) for r in recs}) == 1:
            identical += 1
        else:
            continue

        write_sequences_to_file(recs, aligned_dir / f"{acc}.aligned.fasta")
        prot_clust_to_accession[name] = acc

    logger.info(f"Filtered out {singletons:,} singleton clusters and {identical:,} identical clusters.")

    del sequences_dict_local

    cluster_args = [
        (row['protein_cluster_rep'], row['proteins'], row['index'])
        for row in cluster_info.iter_rows(named=True)
        if f"{acc_prefix}_{row['index']}" not in existing
           and row['protein_cluster_rep'] not in prot_clust_to_accession
    ]
    total_clusters = len(cluster_args)
    total_proteins = sum(len(p) for _, p, _ in cluster_args)
    logger.info(f"There are {total_clusters:,} clusters ({total_proteins:,} proteins) to align.")

    WAVES = 10
    GROUP_SIZE = 1000
    workers = min(threads, os.cpu_count())

    wave_size = max(1, math.ceil(total_clusters / WAVES))
    processed_clusters = 0
    start_all = time.time()

    for wstart in range(0, total_clusters, wave_size):
        wave_args = cluster_args[wstart : wstart + wave_size]

        # gather seqs only for this wave
        wave_ids = {pid for _, plist, _ in wave_args for pid in plist}
        seqs_wave = extract_sequences_from_fasta(fasta, wave_ids)

        groups = [wave_args[i:i+GROUP_SIZE] for i in range(0, len(wave_args), GROUP_SIZE)]
        init_args = (seqs_wave, acc_prefix, output_dir, 1)

        with Pool(workers, initializer=init_worker, initargs=init_args) as pool:
            for results in pool.imap(super_worker, groups):
                for name, _, acc in results:
                    prot_clust_to_accession[name] = acc
                processed_clusters += len(results)
                if processed_clusters % 10000 == 0:
                    rate = processed_clusters / (time.time() - start_all)
                    eta  = (total_clusters - processed_clusters) / rate / 3600
                    logger.info(f"{processed_clusters:,}/{total_clusters:,} clusters aligned "
                                f"(η≈{eta:3.1f} h)")

        # clear per‑wave seqs after the pool exits
        del seqs_wave
        gc.collect()

    # Check for clusters with empty or no alignment output
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