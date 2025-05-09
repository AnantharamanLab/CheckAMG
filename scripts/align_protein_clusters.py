#!/usr/bin/env python3

import os
import sys
import math
import resource
import platform
from pyfastatools import Parser, write_fasta
from pyfastx import Fasta
from collections import namedtuple
from pathlib import Path
import logging
from multiprocessing import Pool
from functools import lru_cache
import time
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import pymuscle5
import pyfamsa
import gc

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

print("========================================================================\n            Step 20/22: Align the filtered protein clusters             \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n            Step 20/22: Align the filtered protein clusters             \n========================================================================\n")


# Module‐level globals for worker processes
sequences_dict = None
_acc_prefix    = None
_output_dir    = None
_aligner_fast  = None
_aligner_full  = None 
_FASTA_IDX     = None
LARGE_CLUSTER_SIZE = 128
Rec = namedtuple("Rec", ["header", "seq"])
class Header:
    """Minimal FASTA header object compatible with write_fasta()"""
    def __init__(self, name, desc=""):
        self.name = name
        self.desc = desc
    
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

def extract_sequences_indexed(ids):
    """
    Retrieve amino acid sequences from pyfastx index for a set of protein IDs.
    Removes terminal stop codons if present.
    Returns a dict of {protein_id: Rec}
    """
    d = {}
    for pid in ids:
        try:
            seq = _FASTA_IDX[pid].seq
            d[pid] = Rec(Header(pid), seq.replace("*", ""))
        except KeyError:
            continue
    return d

def write_sequences_to_file(seq_list, path_out):
    """
    Write a list of Rec sequences to a FASTA file.
    """
    with path_out.open('w') as outfile:
        for seq in seq_list:
            write_fasta(seq, outfile)
            
@lru_cache(maxsize=100_000)
def _mk_seq(name: bytes, seq: bytes):
    """
    Create a PyMuscle5 Sequence object from name and sequence.
    Cached to avoid recomputing for identical sequences.
    """
    return pymuscle5.Sequence(name, seq)

def init_worker(acc_prefix, output_dir, aligner_threads):
    """
    Initialize aligners and global paths for multiprocessing workers.
    """
    global _acc_prefix, _output_dir, _aligner_fast, _aligner_full, _famsa
    _acc_prefix  = acc_prefix
    _output_dir  = output_dir
    _aligner_fast = pymuscle5.Aligner(threads=1, refine_iterations=0)
    _aligner_full = pymuscle5.Aligner(threads=1)
    _famsa = pyfamsa.Aligner(threads=aligner_threads)

def align_cluster_worker(args):
    """
    Aligns one protein cluster:
    - Skip singleton or identical sequences
    - Use pyFAMSA for large clusters
    - Use PyMuscle5 for smaller clusters (faster or refined depending on size)
    Returns (cluster_name, accession).
    """
    cluster_name, protein_list, idx = args
    acc   = f"{_acc_prefix}_{idx}"
    out_f = Path(_output_dir, "aligned", f"{acc}.aligned.fasta")

    recs = [sequences_dict[p] for p in protein_list if p in sequences_dict]
    if not recs:
        return None

    if len(recs) == 1 or len({r.seq for r in recs}) == 1:
        with out_f.open("w") as fh:
            for r in recs: write_fasta(r, fh)
        return (cluster_name, acc)
    
    # Use pyFAMSA for large clusters
    if len(recs) > LARGE_CLUSTER_SIZE:
        seqs = [pyfamsa.Sequence(r.header.name.encode(), r.seq.encode()) for r in recs]
        msa  = _famsa.align(seqs) # thread‑safe, no files
        with out_f.open("w") as fh:
            for s in msa: # PyFAMSA MSA iterable
                fh.write(f">{s.id.decode()}\n{s.sequence.decode()}\n")
        return (cluster_name, acc)
    
    # Use PyMuscle5 for remaining clusters
    seq_objs = [_mk_seq(r.header.name.encode(), r.seq.encode()) for r in recs]
    # aligner  = _aligner_fast if len(recs) <= 8 else _aligner_full
    aligner  = _aligner_full
    msa = aligner.align(seq_objs)
    with out_f.open("w") as fh:
        for s in msa.sequences:
            fh.write(f">{s.name.decode()}\n{s.sequence.decode()}\n")
    return (cluster_name, acc)

def super_worker(group):
    """
    Batch-process a group of cluster alignment tasks to reduce IPC cost.
    """
    return [align_cluster_worker(a) for a in group]

def _fmt_eta(hours_remaining: float) -> str:
    """Return ETA as dd:hh:mm:ss (some may be 0)."""
    total_seconds = int(round(hours_remaining * 3600))
    days,  sec = divmod(total_seconds, 86_400) # 86400 s = 24 h
    hours, sec = divmod(sec, 3_600)
    minutes, sec = divmod(sec, 60)

    return f"{days}d:{hours:02d}h:{minutes:02d}m:{sec}s"
        
def main():
    cluster_file = snakemake.params.cluster_file
    dblookup_file = snakemake.params.dblookup
    input_prot_dir = snakemake.params.input_prot_dir
    confidence_levels = snakemake.params.confidence_levels
    acc_prefix = snakemake.params.acc_prefix
    prot_clust_to_accession_path = snakemake.params.prot_clust_to_accession
    output_dir = snakemake.params.wdir
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem

    logger.info("Starting protein cluster alignment with PyMuscle5...")
    set_memory_limit(mem_limit)
    logger.debug(f"Memory limit set to {mem_limit:,} GB.")
    
    # Build paths for each confidence level file
    input_faa_paths = [os.path.join(input_prot_dir, f"{level}_confidence_viral.faa") for level in confidence_levels]
    logger.debug(f"Input protein files: {input_faa_paths}")
    assert len(input_faa_paths) > 0, "No input protein files found. Please check the input directory."

    if len(input_faa_paths) == 1:
        fasta = input_faa_paths[0]
        logger.debug(f"Using {fasta} as the input protein file.")
    elif len(input_faa_paths) == 2:
        fasta = os.path.join(input_prot_dir, f"{"_".join(confidence_levels)}_confidence_viral.faa")
        unique_records = {}
        for faa_path in input_faa_paths:
            # Parse each fasta file
            for record in Parser(faa_path):
                # Use the record's id as a key for uniqueness;
                # if a record with the same ID is already seen, it is skipped.
                unique_records[record.header.name] = record
        # Write the union (non-redundant records) to the output fasta file.
        with open(fasta, "w") as outfile:
            for record in unique_records.values():
                write_fasta(record, outfile)
        logger.debug(f"Combined protein file created: {fasta}")
        logger.debug(f"Using {fasta} as the input protein file.")
    elif len(input_faa_paths) == 3:
        fasta = os.path.join(input_prot_dir, "filtered_proteins.faa")
        logger.debug(f"Using {fasta} as the input protein file.")
    else:
        fasta = None
        logger.error("More than 3 or less than 1 confidence levels provided. Please check the input files.")
        raise ValueError("More than 3 or less than 1 confidence levels provided. Please check the input files.")

    # build pyfastx index once
    global _FASTA_IDX
    _FASTA_IDX = Fasta(str(fasta), build_index=True)
    
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
    sequences_dict_local = extract_sequences_indexed(protein_ids)
    total_loaded = len(sequences_dict_local)
    logger.debug(f"Loaded {total_loaded:,} protein sequences in memory initially.")
    Path(output_dir, 'aligned').mkdir(parents=True, exist_ok=True)

    # Existing alignments on disk
    aligned_dir = Path(output_dir) / 'aligned'
    existing = {
        f.stem.replace(".aligned", "")
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
    workers = min(threads, os.cpu_count())
    GB_PER_WAVE = 3
    est_bytes = total_proteins * 250 # rough upper bound
    max_in_ram = (GB_PER_WAVE * workers) * (1 << 30)
    WAVES = max(1, math.ceil(est_bytes / max_in_ram))
    GROUP_SIZE = max( min(1_000, # never exceed 1 000
                        math.ceil(total_clusters / (workers * 20))),
                    10) # but at least 10

    wave_size  = math.ceil(total_clusters / WAVES)
    logger.debug(f"Auto‑tuned WAVES={WAVES}, GROUP_SIZE={GROUP_SIZE}, wave_size={wave_size}")

    start_all = time.time()
    finished_total = 0
    if total_clusters >= 100_000:
        LOG_EVERY = 10_000
        next_log = 10_000
    elif total_clusters >= 10_000:
        LOG_EVERY = 1_000
        next_log = 1_000
    elif total_clusters >= 1_000:
        LOG_EVERY = 100
        next_log = 100
    else:
        LOG_EVERY = 10
        next_log = 10

    for wstart in range(0, total_clusters, wave_size):
        wave      = cluster_args[wstart : wstart + wave_size]

        # pull only the sequences required for this wave into shared memory
        wave_ids = {pid for _, plist, _ in wave for pid in plist}
        global sequences_dict
        sequences_dict = extract_sequences_indexed(wave_ids)

        # slice the wave into smaller groups for the worker pool
        groups = [wave[i : i + GROUP_SIZE] for i in range(0, len(wave), GROUP_SIZE)]

        with Pool(workers, initializer=init_worker,
                initargs=(acc_prefix, output_dir, 1)) as pool:

            for grp_idx, results in enumerate(pool.imap_unordered(super_worker, groups)):
                finished_total += len(results) # true completions

                for res in results: # store accessions
                    if res:
                        name, acc = res
                        prot_clust_to_accession[name] = acc

                while finished_total >= next_log:
                    rate = finished_total / (time.time() - start_all)
                    hrs = (total_clusters - finished_total) / rate / 3600
                    logger.info(f"{finished_total:,}/{total_clusters:,} clusters aligned (estimated time remaining: {_fmt_eta(hrs)})...")
                    next_log += LOG_EVERY
                    
        sequences_dict = None
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