#!/usr/bin/env python3

import os
import sys
import resource
import platform
import logging
import pyhmmer
from pyfastatools import Parser
from concurrent.futures import ThreadPoolExecutor, as_completed
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl

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

print("========================================================================\n Step 21/22: Build HMM profiles from the aux protein cluster alignments\n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n Step 21/22: Build HMM profiles from the aux protein cluster alignments\n========================================================================\n")

def read_prot_clust_to_accession(file_path):
    """
    Reads the protein cluster to accession mapping file.
    """
    return pl.read_csv(file_path, separator='\t')

def build_hmm_profile(aligned_fasta, output_dir, accession):
    """Builds and saves an HMM profile for a given aligned fasta file."""
    cluster_id = os.path.splitext(os.path.basename(aligned_fasta))[0].replace(".aligned", "")
    output_path = os.path.join(output_dir, f"{accession}.hmm")

    try:
        # Initialize builder and background for each thread
        alphabet = pyhmmer.easel.Alphabet.amino()
        builder = pyhmmer.plan7.Builder(alphabet)
        background = pyhmmer.plan7.Background(alphabet)

        with pyhmmer.easel.MSAFile(aligned_fasta, digital=True, alphabet=alphabet) as msa_file:
            msa = msa_file.read()
            msa.name = msa.accession = accession.encode()
            hmm, _, _ = builder.build_msa(msa, background)
            with open(output_path, 'wb') as out:
                hmm.write(out)
    except pyhmmer.errors.UnexpectedError as e:
        logger.error(f"UnexpectedError with {aligned_fasta}: {e}")
        raise e
    except ValueError as e:
        logger.error(f"ValueError with {aligned_fasta}: {e}")
        raise e
    except Exception as e:
        logger.error(f"Error building HMM for {aligned_fasta}: {e}")
        raise e

    return cluster_id, accession

def extract_sequences_from_aligned_fasta(aligned_fasta):
    """
    Extracts protein sequences and their IDs from an aligned FASTA file.
    """
    sequences = []
    for record in Parser(aligned_fasta):
        sequences.append((record.header.name, str(record.seq)))
    return sequences

def parallel_build_hmm_profiles(prot_clust_to_accession, aligned_dir, output_dir, threads):
    output_data = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_cluster = {}

        for row in prot_clust_to_accession.iter_rows(named=True):
            protein_cluster = row['protein_cluster_rep']
            accession = row['accession']
            aligned_fasta = os.path.join(aligned_dir, f"{accession}.aligned.fasta")

            if os.path.exists(aligned_fasta):
                # Pass aligned_fasta, output_dir, and accession directly
                future = executor.submit(
                    build_hmm_profile,
                    aligned_fasta=aligned_fasta,
                    output_dir=output_dir,
                    accession=accession
                )
                future_to_cluster[future] = (protein_cluster, accession, aligned_fasta)

        for future in as_completed(future_to_cluster):
            protein_cluster, accession, aligned_fasta = future_to_cluster[future]
            try:
                cluster_id, hmm_accession = future.result()
                sequences = extract_sequences_from_aligned_fasta(aligned_fasta)
                for protein_id, sequence in sequences:
                    output_data.append((protein_id, protein_cluster, hmm_accession))
            except Exception as exc:
                logger.error(f"Error building HMM for {protein_cluster}: {exc}")
                raise exc

    return output_data

def main():
    aligned_dir = snakemake.params.alignments_dir
    output_dir = snakemake.params.outdir
    prot_clust_to_accession_path = snakemake.params.prot_clust_to_accession
    output_table_path = snakemake.params.output_table
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Building HMM profiles from auxiliary protein alignments...")
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Load pyhmmer alphabet, builder, and background
    alphabet = pyhmmer.easel.Alphabet.amino()
    builder = pyhmmer.plan7.Builder(alphabet)
    background = pyhmmer.plan7.Background(alphabet)
    
    # Load the protein cluster to accession mapping
    prot_clust_to_accession = read_prot_clust_to_accession(prot_clust_to_accession_path)

    # Build HMM profiles in parallel
    output_data = parallel_build_hmm_profiles(
        prot_clust_to_accession=prot_clust_to_accession,
        aligned_dir=aligned_dir,
        output_dir=output_dir,
        threads=threads
    )

    # Create a DataFrame for the output table
    output_df = pl.DataFrame(output_data, schema=['Protein', 'protein_cluster_rep', 'hmm'], orient="row")

    # Add a new column with the numeric part extracted from the `hmm` column
    output_df = output_df.with_columns(
        output_df["hmm"]
        .str.extract(r"phaux_(\d+)", 1) # Extract the numeric part from `hmm`
        .cast(pl.Int64) # Convert to integer for sorting
        .alias("hmm_numeric")
    )

    # Sort by the numeric column then drop it
    output_df = output_df.sort("hmm_numeric").drop("hmm_numeric")

    # Calculate HMM families statistics:
    # Group by protein_cluster_rep and count the number of proteins for each family.
    grouped = output_df.group_by("protein_cluster_rep").len()
    singleton_count = grouped.filter(pl.col("len") == 1).height
    multi_count = grouped.filter(pl.col("len") > 1).height
    total_families = grouped.height
    logger.info(f"Total HMM profiles: {total_families:,}")
    logger.info(f"Singleton families: {singleton_count:,}")
    logger.info(f"Multi-sequence families: {multi_count:,}")

    # Save the output table to a file
    output_df.write_csv(output_table_path, separator='\t')
    
    logger.info("Building HMM profiles and output table completed.")
    logger.info(f"Total HMM profiles built: {len(output_df['hmm'].unique()):,}")
    
if __name__ == "__main__":
    main()

