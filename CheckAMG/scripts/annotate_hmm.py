#!/usr/bin/env python3

import os
import sys
import resource
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pathlib import Path
from pyhmmer import easel, plan7, hmmer
import uuid
from datetime import datetime
from pyfastatools import Parser, write_fasta
import math
from tqdm import tqdm

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

def aggregate_sequences(protein_dir):
    all_sequences = []
    protein_dir = Path(protein_dir)
    for fasta_file in protein_dir.rglob("*"):
        if fasta_file.suffix.lower() in (".faa", ".fasta"):
            all_sequences.extend(Parser(str(fasta_file)).all())
    return all_sequences

def split_aggregated_sequences(all_sequences, chunk_size):
    for i in range(0, len(all_sequences), chunk_size):
        yield all_sequences[i:i + chunk_size]

def determine_chunk_size(n_sequences, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8):
    total_bytes = n_sequences * est_bytes_per_seq
    allowed_bytes = max_chunk_fraction * mem_limit * (1024**3)
    n_chunks = max(1, math.ceil(total_bytes / allowed_bytes))
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

def get_kegg_threshold(hmm_id):
    return KEGG_THRESHOLDS.get(hmm_id, None)

def hmmsearch_serial(batch_key, batch_fasta, db_path, seq_lengths, out_dir, min_coverage, min_score, min_bitscore_fraction, evalue, cpus):
    outfile = Path(out_dir) / f"{batch_key}_search.tsv"
    alphabet = easel.Alphabet.amino()
    hmm_list = list(plan7.HMMFile(db_path))
    db = assign_db(db_path)
    with open(outfile, 'w') as out, easel.SequenceFile(batch_fasta, digital=True, alphabet=alphabet) as seqs:
        for hits in hmmer.hmmsearch(queries=hmm_list, sequences=seqs, E=0.1, cpus=cpus):
            hmm = hits.query
            hmm_id = extract_query_info(hits, db_path)
            for hit in hits:
                hit_name = hit.name.decode()
                for dom in hit.domains.included:
                    a = dom.alignment
                    alignment_len = a.target_to - a.target_from + 1
                    coverage = alignment_len / seq_lengths[hit_name]
                    # Apply GA/TC cutoffs where available
                    if db == "Pfam" and hmm.cutoffs.gathering is not None:
                        if dom.score < hmm.cutoffs.gathering1:
                            continue
                    elif db == "KEGG":
                        kegg_thresh = get_kegg_threshold(hmm_id)
                        if kegg_thresh is not None and dom.score < kegg_thresh:
                            if hit.evalue > evalue or dom.score < min_bitscore_fraction * kegg_thresh:
                                continue
                        elif kegg_thresh is None and (dom.score < min_score or coverage < min_coverage):
                            continue
                    elif db == "FOAM" and hmm.cutoffs.trusted is not None:
                        if dom.score < hmm.cutoffs.trusted1:
                            if hit.evalue > evalue or dom.score < min_bitscore_fraction * hmm.cutoffs.trusted1:
                                continue
                    elif db == "METABOLIC" and hmm.cutoffs.gathering is not None:
                        if dom.score < hmm.cutoffs.gathering1:
                            continue
                    else:
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

    logger.info("Protein HMM searches starting...")
    run_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    tmp_dir = Path(wdir) / f"hmmsearch_tmp_{run_id}"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    aggregated = aggregate_sequences(protein_dir)
    seq_lengths = {rec.header.name: len(rec.seq) for rec in aggregated}

    set_memory_limit(mem_limit)
    logger.debug(f"Memory limit set to {mem_limit} GB")

    priority_order = ["KEGG", "FOAM", "PHROG", "VOG", "Pfam", "eggNOG", "dbCAN", "METABOLIC", "user_custom"]
    hmm_paths = sorted([Path(db_dir) / f for f in os.listdir(db_dir) if f.endswith((".H3M", ".h3m"))],
                    key=lambda x: priority_order.index(assign_db(x)) if assign_db(x) in priority_order else float('inf'))

    N = len(aggregated)
    chunk_size = determine_chunk_size(N, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8)
    batch_files = []
    batch_keys = []
    for idx, batch in enumerate(split_aggregated_sequences(aggregated, chunk_size)):
        batch_key = f"seqbatch_{idx}"
        batch_fasta = Path(tmp_dir) / f"{batch_key}.faa"
        with open(batch_fasta, 'w') as f:
            for rec in batch:
                write_fasta(rec, f)
        batch_files.append(str(batch_fasta))
        batch_keys.append(batch_key)
    logger.info(f"Splitting {N:,} sequences into {len(batch_keys):,} batches of {chunk_size:,}")
    
    result_paths = []
    for db_path in hmm_paths:
        db_name = assign_db(db_path)
        logger.info(f"Running hmmsearches against {db_name} profile HMMs...")
        if len(batch_keys) > 1:
            for batch_key, batch_fasta in tqdm(zip(batch_keys, batch_files), total=len(batch_keys), desc=f"HMMsearches ({db_name})", unit="batch"):
                out_path = hmmsearch_serial(
                    batch_key=f"{db_path.stem}_{batch_key}",
                    batch_fasta=batch_fasta,
                    db_path=str(db_path),
                    seq_lengths=seq_lengths,
                    out_dir=tmp_dir,
                    min_coverage=cov_fraction,
                    min_score=minscore,
                    min_bitscore_fraction=min_bitscore_fraction,
                    evalue=evalue,
                    cpus=num_threads
                )
                result_paths.append((out_path, db_path))
        else:
            for batch_key, batch_fasta in zip(batch_keys, batch_files):
                out_path = hmmsearch_serial(
                    batch_key=f"{db_path.stem}_{batch_key}",
                    batch_fasta=batch_fasta,
                    db_path=str(db_path),
                    seq_lengths=seq_lengths,
                    out_dir=tmp_dir,
                    min_coverage=cov_fraction,
                    min_score=minscore,
                    min_bitscore_fraction=min_bitscore_fraction,
                    evalue=evalue,
                    cpus=num_threads
                )
                result_paths.append((out_path, db_path))

    logger.info("Combining and filtering HMMsearch results...")
    filtered_paths = []
    for result_path, db_path in result_paths:
        filtered_path = result_path.replace("_search.tsv", "_filtered.tsv")
        logger.debug(f"Filtering results from {result_path} using {db_path} to {filtered_path}")
        filter_hmm_results(result_path, db_path, filtered_path)
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
    merged_df = combined_df.rename({"hmm_id": "id"}).join(
        vscores_df, left_on="hmm_id_norm", right_on="id_norm", how="left"
    )
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

    logger.info("Protein HMM searches completed.")

if __name__ == "__main__":
    main()