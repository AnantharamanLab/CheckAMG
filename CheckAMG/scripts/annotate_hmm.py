#!/usr/bin/env python3

import os
import sys
import resource
import platform
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pathlib import Path
from pyhmmer import easel, plan7, hmmer
import load_prot_paths
import uuid
from multiprocessing import Pool
from datetime import datetime
from pyfastatools import Parser, write_fasta
import math
from math import ceil
from collections import defaultdict
from tqdm import tqdm
from functools import partial
from random import shuffle

# Global cache of pre-loaded HMM models, shared by forked workers via copy-on-write
HMM_MODELS = {}

# Global caches for thresholds
KEGG_THRESHOLDS = {}

# Load KEGG thresholds
KEGG_THRESHOLDS_PATH = snakemake.params.kegg_cutoff_file
if Path(KEGG_THRESHOLDS_PATH).exists():
    df = pl.read_csv(KEGG_THRESHOLDS_PATH)
    KEGG_THRESHOLDS = dict(zip(df["id"].to_list(), df["threshold"].to_list()))


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

print("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================\n")

def assign_db(db_path):
    if "KEGG" in str(db_path) or "kegg" in str(db_path) or "kofam" in str(db_path):
        return "KEGG"
    elif "FOAM" in str(db_path) or "foam" in str(db_path):
        return "FOAM"
    elif "Pfam" in str(db_path) or "pfam" in str(db_path):
        return "Pfam"
    elif "dbcan" in str(db_path) or "dbCAN" in str(db_path) or "dbCan" in str(db_path):
        return "dbCAN"
    elif "METABOLIC_custom" in str(db_path) or "metabolic_custom" in str(db_path):
        return "METABOLIC"
    elif "VOG" in str(db_path) or "vog" in str(db_path):
        return "VOG"
    elif "eggNOG" in str(db_path) or "eggnog" in str(db_path):
        return "eggNOG"
    elif "PHROG" in str(db_path) or "phrog" in str(db_path):
        return "PHROG"
    elif "user_custom" in str(db_path):
        return "user_custom"
    else:
        return None 

def extract_query_info(hits, db_path):
    if "Pfam" in str(db_path) or "pfam" in str(db_path):
        hmm_id = hits.query.accession.decode()
    elif "FOAM" in str(db_path) or "foam" in str(db_path):
        hmm_id = hits.query.accession.decode()
    elif "eggNOG" in str(db_path) or "eggnog" in str(db_path):
        hmm_id = hits.query.name.decode().split(".")[0]
    else:
        query_name = hits.query.name.decode()
        if ".wlink.txt.mafft" in query_name:
            hmm_id = query_name.split(".")[1]
        else:
            hmm_id = query_name.replace("_alignment", "").replace(".mafft", "").replace(".txt", "").replace(".hmm", "").replace("_protein.alignment", "")
    return hmm_id

def aggregate_sequences(prots):
    all_sequences = []
    for fasta_file in prots:
        all_sequences.extend(Parser(fasta_file).all())
    return all_sequences

def split_aggregated_sequences(all_sequences, chunk_size):
    for i in range(0, len(all_sequences), chunk_size):
        yield all_sequences[i:i + chunk_size]

def determine_chunk_size(n_sequences, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8, n_processes=1):
    total_bytes = n_sequences * est_bytes_per_seq
    allowed_bytes = max_chunk_fraction * mem_limit * (1024**3) / n_processes
    n_chunks = max(1, math.ceil(total_bytes / allowed_bytes), n_processes)
    return math.ceil(n_sequences / n_chunks)

def filter_hmm_results(tsv_path, hmm_path, out_path):
    db = assign_db(hmm_path)
    results = {}

    # Load hits and keep only the best one per sequence
    with open(tsv_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 7:
                continue
            sequence, hmm_id, evalue, score, length, start, end = parts
            evalue = float(evalue)
            score = float(score)
            start = int(start)
            end = int(end)
            key = sequence

            current_best = results.get(key)
            new_hit = (hmm_id, score, int(length), int(start), int(end), evalue)

            if current_best is None:
                results[key] = new_hit
            else:
                # Keep hit with better (lower) evalue, or higher score if evalue is the same
                if new_hit[5] < current_best[5] or (new_hit[5] == current_best[5] and new_hit[1] > current_best[1]):
                    results[key] = new_hit

    # Write best hit per sequence
    with open(out_path, 'w') as out:
        out.write("hmm_id\tsequence\tscore\tcoverage\tdb\n")
        for seq, hit in results.items():
            hmm_id, score, length, start, end, _ = hit
            coverage = (end - start + 1) / length
            out.write(f"{hmm_id}\t{seq}\t{score:.6f}\t{coverage:.3f}\t{db}\n")

def _hmm_worker(args):
    return hmmsearch_worker(*args)

def get_kegg_threshold(hmm_id):
    return KEGG_THRESHOLDS.get(hmm_id, None)

def hmmsearch_worker(key, seq_path, db_str, seq_lengths, out_dir, min_coverage, min_score, min_bitscore_fraction, evalue, cpus):
    outfile = Path(out_dir) / f"{key}_search.tsv"
    errfile = Path(out_dir) / f"{key}_search.err"
    alphabet = easel.Alphabet.amino()
    hmm_list = HMM_MODELS[db_str]
    db = assign_db(db_str)

    with open(outfile, 'w') as out, open(errfile, 'w') as err, \
         easel.SequenceFile(seq_path, digital=True, alphabet=alphabet) as seqs:
        for hits in hmmer.hmmsearch(queries=hmm_list, sequences=seqs, E=0.1, cpus=cpus): # Use a permissive evalue to ensure reproducibility when chunk size changes due to different memory limits
            hmm = hits.query  # plan7.HMM object
            hmm_id = extract_query_info(hits, db_str)
            for hit in hits:
                hit_name = hit.name.decode()
                for dom in hit.domains.included:
                    a = dom.alignment
                    alignment_len = a.target_to - a.target_from + 1
                    coverage = alignment_len / seq_lengths[hit_name]

                    # Apply GA/TC cutoffs where available
                    if db == "Pfam" and hmm.cutoffs.gathering is not None:
                        if dom.score < hmm.cutoffs.gathering1: # use sequence-level GA, not domain GA
                            continue
                        
                    elif db == "KEGG":
                        kegg_thresh = get_kegg_threshold(hmm_id)
                        if kegg_thresh is not None and dom.score < kegg_thresh:
                            # Apply a heuristic like `anvi-run-kegg-kofams` from Anvi'o does, since KEGG thresholds are sometimes too strict
                            if hit.evalue > evalue or dom.score < min_bitscore_fraction * kegg_thresh:
                                continue
                        elif kegg_thresh is None and (dom.score < min_score or coverage < min_coverage):
                            continue
                        
                    elif db == "FOAM" and hmm.cutoffs.trusted is not None:
                        if dom.score < hmm.cutoffs.trusted1: # use sequence-level TC, not domain TC
                            # Apply the same heuristic as KEGG, since FOAM comes from KEGG
                            if hit.evalue > evalue or dom.score < min_bitscore_fraction * hmm.cutoffs.trusted1:
                                continue
                            
                    elif db == "METABOLIC" and hmm.cutoffs.gathering is not None:
                        if dom.score < hmm.cutoffs.gathering1: # use sequence-level GA, not domain GA
                            continue
                        
                    else:
                        # fallback filtering
                        if dom.score < min_score or coverage < min_coverage:
                            continue

                    out.write(f"{hit_name}\t{hmm_id}\t{hit.evalue:.1E}\t{dom.score:.6f}\t{hit.length}\t{a.target_from}\t{a.target_to}\n")

    return str(outfile)
        
def main():
    protein_dir = snakemake.params.protein_dir
    wdir = snakemake.params.wdir
    hmm_vscores = snakemake.params.hmm_vscores
    cov_fraction = snakemake.params.cov_fraction
    db_dir = snakemake.params.db_dir
    output = Path(snakemake.params.vscores)
    all_hmm_results = Path(snakemake.params.all_hmm_results)
    num_threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    minscore = snakemake.params.min_bitscore
    min_bitscore_fraction = snakemake.params.min_bitscore_fraction_heuristic
    evalue = snakemake.params.max_evalue

    logger.info("Protein HMM alignments starting...")

    run_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    tmp_dir = Path(wdir) / f"hmmsearch_tmp_{run_id}"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    aggregated = aggregate_sequences(load_prot_paths.load_prots(protein_dir))
    seq_lengths = defaultdict(int)
    seq_lengths.update({rec.header.name: len(rec.seq) for rec in aggregated})
    logger.debug(f"Memory limit set to {mem_limit} GB")
    set_memory_limit(mem_limit)

    priority_order = ["KEGG", "FOAM", "PHROG", "VOG", "Pfam", "eggNOG", "dbCAN", "METABOLIC", "user_custom"]
    hmm_paths = sorted([Path(db_dir) / f for f in os.listdir(db_dir) if f.endswith(('.H3M', '.h3m'))],
                    key=lambda x: priority_order.index(assign_db(x)) if assign_db(x) in priority_order else float('inf'))
    
    
    for db in hmm_paths:
        db_str = str(db)
        HMM_MODELS[db_str] = list(plan7.HMMFile(db_str))

    aggregated = aggregate_sequences(load_prot_paths.load_prots(protein_dir))
    seq_lengths = {rec.header.name: len(rec.seq) for rec in aggregated}

    db_counts = {db: sum(1 for _ in plan7.HMMFile(str(db))) for db in hmm_paths}
    total_models = sum(db_counts.values())
    OVERSUB = 10
    total_tasks = num_threads * OVERSUB
    jobs = []
    N = len(aggregated)

    for db in hmm_paths:
        db_str = str(db)
        tasks = max(1, round(total_tasks * db_counts[db] / total_models))
        cs = ceil(N / tasks)
        for i in range(0, N, cs):
            chunk = aggregated[i:i+cs]
            chunk_file = Path(tmp_dir) / f"chunk_{db.stem}_{i//cs}.faa"
            with open(chunk_file, 'w') as f:
                for rec in chunk:
                    write_fasta(rec, f)
            jobs.append((
                f"{db.stem}_{i//cs}", chunk_file, db_str,
                seq_lengths, tmp_dir, cov_fraction, minscore, min_bitscore_fraction, evalue, 1
            ))
    shuffle(jobs) # Shuffle jobs so the big databases don't always run all at first
    
    logger.info(f"Running HMMsearch with {num_threads} maximum jobs in parallel...")
    logger.debug(f"Using a minimum bit score of {minscore} and a minimum coverage fraction of {cov_fraction} for fallback filtering when database-provided cutoffs are not available.")
    logger.debug(f"Using a minimum bitscore fraction of {min_bitscore_fraction} and maximum E-value of {evalue} for heuristic filtering of HMM hits for KEGG and FOAM HMMs.")
    # Run HMMsearch in parallel
    result_paths = []
    with Pool(processes=num_threads) as pool:
        for res in tqdm(pool.imap_unordered(_hmm_worker, jobs), total=len(jobs), desc="HMMsearches", unit="chunk"):
            result_paths.append(res)

    # Combine and filter result files
    logger.info("Combining and filtering HMMscan results...")
    filtered_paths = []
    for result_path in result_paths:
        hmm_path = result_path.split("_chunk")[0] + ".h3m"
        filtered_path = result_path.replace("_search.tsv", "_filtered.tsv")
        logger.debug(f"Filtering results from {result_path} using {hmm_path} to {filtered_path}")
        filter_hmm_results(result_path, hmm_path, filtered_path)
        filtered_paths.append(filtered_path)

    schema = {
        "hmm_id": pl.Utf8,
        "sequence": pl.Utf8,
        "score": pl.Float64,
        "hmm_coverage": pl.Float64,
        "db": pl.Utf8
    }

    dfs = []
    for p in filtered_paths:
        try:
            df = pl.read_csv(p, separator="\t", schema=schema, ignore_errors=True)
            dfs.append(df)
        except Exception as e:
            logger.warning(f"Failed to load {p}: {e}")

    combined_df = pl.concat(dfs)
    combined_df.write_csv(all_hmm_results, separator="\t")

    vscores_df = pl.read_csv(hmm_vscores, schema_overrides={"id": pl.Utf8, "V-score": pl.Float64, "VL-score": pl.Float64, "db": pl.Categorical, "name": pl.Utf8})
    
    # Normalize Pfam IDs for join (strip version suffix)
    vscores_df = vscores_df.with_columns([
        pl.when(pl.col("db") == "Pfam")
        .then(pl.col("id").str.replace(r"\.\d+$", ""))
        .otherwise(pl.col("id")).alias("id_norm")
    ])
    combined_df = combined_df.with_columns([
        pl.when(pl.col("db") == "Pfam")
        .then(pl.col("hmm_id").str.replace(r"\.\d+$", ""))
        .otherwise(pl.col("hmm_id")).alias("hmm_id_norm")
    ])

    # Join on normalized id columns
    merged_df = combined_df.rename({"hmm_id": "id"}).join(
        vscores_df, left_on="hmm_id_norm", right_on="id_norm", how="left"
    )

    # Keep original columns for downstream logic
    merged_df = merged_df.with_columns([
        pl.col("id").alias("hmm_id")
    ])
    merged_df = merged_df.filter(pl.col("V-score").is_not_null())
    cols_to_drop = ["name", "db_right", "id", "id_norm", "hmm_id_norm"]
    for col in cols_to_drop:
        if col in merged_df.columns:
            merged_df = merged_df.drop(col)
    merged_df = merged_df.sort(["sequence", "score", "V-score", "db"])
    merged_df.write_csv(output, separator="\t")
    
    for f in tmp_dir.iterdir():
        f.unlink()
    tmp_dir.rmdir()

    logger.info("Protein HMM alignments completed.")

if __name__ == "__main__":
    main()