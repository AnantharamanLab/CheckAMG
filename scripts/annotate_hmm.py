#!/usr/bin/env python3

import os
import sys
import resource
import platform
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import collections
import math
import pyhmmer
from pyhmmer import easel, plan7
from pyhmmer.errors import AllocationError, EaselError
from pyfastatools import Parser, write_fasta
import load_prot_paths
import uuid
import gc
from datetime import datetime

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
        
if snakemake.params.build_or_annotate =="build":
    print("========================================================================\n                Step 5/22: Assign functions to proteins                 \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n                Step 5/22: Assign functions to proteins                 \n========================================================================\n")
elif snakemake.params.build_or_annotate == "annotate":
    print("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================\n")

def load_hmms(hmmdb_path):
    """Load HMM profiles from a given HMM database file."""
    with plan7.HMMFile(hmmdb_path) as hmm_file:
        hmms = list(hmm_file)
    return hmms

def filter_results(results):
    """Filter the search results to keep only the best scoring hit for each sequence within each database."""
    best_results = {}
    for result in results:
        key = (result.sequence, result.db_path)
        if key not in best_results or result.score > best_results[key].score:
            best_results[key] = result
    return list(best_results.values())

def get_hmm_coverage(domain):
    """Calculate the alignment coverage for a given domain."""
    n_aligned_positions = domain.alignment.hmm_to - domain.alignment.hmm_from + 1
    return n_aligned_positions / domain.alignment.hmm_length

def extract_query_info(hits, db_path):
    """Extract query information depending on the database."""
    if "Pfam" in db_path:
        hmm_id = hits.query_accession.decode()
    elif "eggNOG" in db_path:
        hmm_id = hits.query_name.decode().split(".")[0]
    else:
        query_name = hits.query_name.decode()
        if ".wlink.txt.mafft" in query_name:
            hmm_id = query_name.split(".")[1]
        else:
            hmm_id = query_name.replace("_alignment", "").replace(".mafft", "").replace(".txt", "").replace(".hmm", "").replace("_protein.alignment", "")
    return hmm_id

def run_hmmsearch_on_chunk(sequence_chunk, hmm_list, db_path, e_value_threshold, num_cpus, cov_fraction, output_path):
    Result = collections.namedtuple("Result", ["hmm_id", "sequence", "score", "db_path", "hmm_coverage"])
    results = []

    # Create the amino acid alphabet
    aa_alphabet = easel.Alphabet.amino()

    # Read sequences
    with easel.SequenceFile(sequence_chunk, format="fasta", digital=True, alphabet=aa_alphabet) as seqs_file:
        proteins = list(seqs_file)

    # Search each sequence using hmm profiles
    for hits in pyhmmer.hmmsearch(queries=hmm_list, sequences=proteins, cpus=num_cpus, E=e_value_threshold):
        for hit in hits:
            for domain in hit.domains.included:
                hmm_coverage = get_hmm_coverage(domain)
                if hmm_coverage >= cov_fraction:
                    hmm_id = extract_query_info(hits, db_path)
                    results.append(Result(hmm_id, hit.name.decode(), hit.score, db_path, hmm_coverage))
    filtered_results = filter_results(results)

    # Write the results to a temporary file
    with open(output_path, "w") as output_file:
        for result in filtered_results:
            output_file.write(f"{result.hmm_id}\t{result.sequence}\t{result.score}\t{result.db_path}\t{result.hmm_coverage}\n")

    # Clean up memory
    del results, filtered_results, proteins
    gc.collect()
        
def aggregate_sequences(prots):
    """Aggregate all sequences from all fasta files into a single list."""
    all_sequences = []
    for fasta_file in prots:
        sequences = list(Parser(fasta_file).all())
        all_sequences.extend(sequences)
    return all_sequences

def split_aggregated_sequences(all_sequences, chunk_size):
    """Split aggregated sequences into evenly sized chunks."""
    for i in range(0, len(all_sequences), chunk_size):
        yield all_sequences[i:i + chunk_size]

def determine_chunk_size(n_sequences, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8):
    """
    Determine chunk size such that the total memory used per chunk is
    no more than max_chunk_fraction (e.g., 80%) of the maximum memory.
    
    Parameters:
    - n_sequences: Total number of sequences.
    - mem_limit: Maximum allowable memory in GB.
    - est_bytes_per_seq: Estimated memory used per sequence in bytes.
    - max_chunk_fraction: Fraction of mem_limit to target per chunk.
    
    Returns:
    - The number of sequences per chunk.
    """
    # Total estimated memory needed for all sequences in bytes.
    total_bytes = n_sequences * est_bytes_per_seq
    
    # Allowed memory per chunk in bytes.
    allowed_bytes = max_chunk_fraction * mem_limit * (1024**3)
    
    # Determine the number of chunks required so that each chunk's memory usage does not exceed the allowed_bytes.
    n_chunks = math.ceil(total_bytes / allowed_bytes)
    n_chunks = max(1, n_chunks) # Ensure at least one chunk.
    
    # Compute the chunk size by splitting the total sequences into n_chunks.
    chunk_size = math.ceil(n_sequences / n_chunks)
    return chunk_size

def process_database(hmmdb_path, aggregated_sequences, output_dir, e_value_threshold, num_cpus, cov_fraction, chunk_size, run_id):
    hmm_list = load_hmms(hmmdb_path)

    chunk_id = 0
    for chunk in split_aggregated_sequences(aggregated_sequences, chunk_size):
        chunk_file = os.path.join(output_dir, f"{run_id}_chunk_{chunk_id}.fasta")
        with open(chunk_file, "w") as f:
            for record in chunk:
                write_fasta(record, f)

        output_path = os.path.join(output_dir, f"{run_id}_results_{os.path.basename(hmmdb_path)}_{chunk_id}.tsv")
        run_hmmsearch_on_chunk(chunk_file, hmm_list, hmmdb_path, e_value_threshold, num_cpus, cov_fraction, output_path)

        # Remove chunk file to save space
        os.remove(chunk_file)
        del chunk
        gc.collect()

        chunk_id += 1

    del hmm_list
    gc.collect()

def combine_results(output_dir, run_id):
    combined_results = []

    for file_name in os.listdir(output_dir):
        if file_name.startswith(f"{run_id}_results_") and file_name.endswith(".tsv"):
            file_path = os.path.join(output_dir, file_name)
            if os.path.getsize(file_path) > 0:  # Skip empty files
                chunk_df = pl.read_csv(file_path, separator="\t", has_header=False)
                chunk_df.columns = ["hmm_id", "sequence", "score", "db_path", "hmm_coverage"]
                combined_results.append(chunk_df)
                del chunk_df
                # gc.collect()

    if combined_results:
        combined_df = pl.concat(combined_results)
        del combined_results
        gc.collect()
        return combined_df
    else:
        logging.error("No HMMsearch results found; all CSV files are empty.")
        raise ValueError("No valid data found in any result files.")

def assign_db(db_path):
    """Assign the database label based on db_path."""
    if "KEGG" in db_path or "kegg" in db_path or "kofam" in db_path:
        return "KEGG"
    elif "Pfam" in db_path or "pfam" in db_path:
        return "Pfam"
    elif "dbcan" in db_path or "dbCAN" in db_path or "dbCan" in db_path:
        return "dbCAN"
    elif "METABOLIC_custom" in db_path or "metabolic_custom" in db_path:
        return "METABOLIC"
    elif "VOG" in db_path or "vog" in db_path:
        return "VOG"
    elif "eggNOG" in db_path or "eggnog" in db_path:
        return "eggNOG"
    elif "PHROG" in db_path or "phrog" in db_path:
        return "PHROG"
    elif "user_custom" in db_path:
        return "user_custom"
    else:
        return None 

def main():
    protein_dir = snakemake.params.protein_dir
    wdir = snakemake.params.wdir 
    hmm_vscores = snakemake.params.hmm_vscores
    cov_fraction = snakemake.params.cov_fraction
    db_dir = snakemake.params.db_dir
    output = snakemake.params.vscores
    all_hmm_results = snakemake.params.all_hmm_results
    num_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem

    logger.info("Protein HMM alignments starting...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    set_memory_limit(mem_limit)
            
    prots = load_prot_paths.load_prots(protein_dir)
    
    # Aggregate sequences from all files
    aggregated_sequences = aggregate_sequences(prots)
    
    # Determine chunk size
    n_sequences = len(aggregated_sequences)
    chunk_size = determine_chunk_size(n_sequences, mem_limit)
    logger.debug(f"Chunk size: {chunk_size:,} sequences per chunk, for {n_sequences//chunk_size} chunks.")
    
    run_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    tmp_dir = os.path.join(wdir, f"pyhmmer_tmp_{run_id}")
    os.makedirs(tmp_dir, exist_ok=True)
    priority_order = ["KEGG", "PHROG", "VOG", "Pfam", "eggNOG", "dbCAN", "METABOLIC", "user_custom"]

    # Load database paths
    hmmdbs = [os.path.join(db_dir, db) for db in os.listdir(db_dir) if db.endswith('.H3M') or db.endswith('.h3m')]

    # Sort databases based on the priority order
    hmmdbs_sorted = sorted(hmmdbs, key=lambda x: priority_order.index(assign_db(x)) if assign_db(x) in priority_order else float('inf'))

    logger.info(f"Loading HMM profiles.")

    for hmmdb_path in hmmdbs_sorted:
        db_label = assign_db(hmmdb_path)
        logger.info(f"Running HMMsearch against database: {db_label}")
        try:
            process_database(hmmdb_path, aggregated_sequences, tmp_dir, 1E-10, num_cpus, cov_fraction, chunk_size, run_id)
        except (AllocationError, EaselError, SystemError, MemoryError) as e:
            logging.error(f"Memory limit exceeded. Please increase the memory limit or reduce the number of input sequences. Error: {e}")
            raise

    logger.info("Combining results from all databases.")
    combined_df = combine_results(tmp_dir, run_id)

    # Format the final DataFrame
    combined_df = combined_df.rename({"hmm_id": "id"})
    combined_df = combined_df.with_columns(pl.col("db_path").map_elements(assign_db, return_dtype=pl.Utf8).alias("db"))
    combined_df = combined_df.drop('db_path')
    combined_df = combined_df.rename({"id": "hmm_id"})
    combined_df = combined_df.sort(['sequence', 'score', 'db', 'hmm_id'])

    # Write all hmm results to file
    combined_df.write_csv(all_hmm_results, separator="\t")

    # Load V-scores CSV
    vscores_df = pl.read_csv(hmm_vscores, schema_overrides={"id": pl.Utf8, "V-score": pl.Float64, "VL-score": pl.Float64, "db": pl.Categorical, "name": pl.Utf8})

    # Merge with V-scores and filter
    merged_df = combined_df.join(vscores_df, left_on='hmm_id', right_on='id', how='left').filter(pl.col("V-score").is_not_null())
    merged_df = merged_df.drop('name', 'db_right')
    merged_df = merged_df.sort(['sequence', 'score', 'V-score', 'db'])
    merged_df.write_csv(output, separator="\t")
    
    # Cleanup: remove temporary files if needed
    for file_name in os.listdir(tmp_dir):
        file_path = os.path.join(tmp_dir, file_name)
        os.remove(file_path)
    os.rmdir(tmp_dir)
    
    logger.info("Protein HMM alignments completed.")

if __name__ == "__main__":
    main()
