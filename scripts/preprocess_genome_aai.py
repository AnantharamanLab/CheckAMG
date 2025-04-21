#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pyfastatools import Parser
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

print("========================================================================\n      Step 11/22: Preprocess alignment output for AAI calculations      \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n      Step 11/22: Preprocess alignment output for AAI calculations      \n========================================================================\n")

def parse_faa_file(faa_file, is_vMAG=False):
    """
    Parses an amino-acid fasta file to obtain protein names and their genomes.
    """ 
    try:
        for header in Parser(faa_file).all_prodigal_headers():
            seq_id, scaffold, gene_number = str(header.name()), str(header.scaffold), str(header.id)
            if is_vMAG:
                genome = os.path.splitext(os.path.basename(faa_file))[0]
            else:
                genome = scaffold
            yield seq_id, gene_number, scaffold, genome
    except RuntimeError as e:
        logger.error(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
        raise RuntimeError(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
    
def gene_to_genome(input_faa_files):
    """
    Obtains gene-to-genome mapping from prodigal-formatted input
    faa(s)

    input_files: path(s) to amino-acid fasta files to be parsed.
    
    Returns a Polars DataFrame with genes mapped to the genomes
    that encode them.
    """
    # Accumulate all faa data
    faa_data = {"protein": [], "genome": [], "scaffold": [], "gene_number": []}

    # Assumes the first file is not vMAG and subsequent files are
    for i, file in enumerate(input_faa_files):
        is_vMAG = i > 0 # First file is not vMAG, others are
        for seq_id, gene_number, scaffold, genome in parse_faa_file(file, is_vMAG = is_vMAG): 
            faa_data["protein"].append(seq_id)
            faa_data["genome"].append(genome)
            faa_data["scaffold"].append(scaffold)
            faa_data["gene_number"].append(int(gene_number))

    # Convert the list of dictionaries (faa_data) to a Polars DataFrame
    gene_to_genome = pl.DataFrame(faa_data).sort(["genome", "scaffold", "gene_number"])
    
    return gene_to_genome

def assign_genome_ids(gene_map: pl.DataFrame) -> pl.DataFrame:
    """
    Takes a gene_map DataFrame with columns ["protein", "genome"].
    Returns a new DataFrame with columns:
        ["protein", "genome", "genome_id"] 
    where genome_id is an integer code for each distinct genome.
    """
    # Extract the distinct genomes
    unique_genomes = gene_map.select("genome").unique().sort("genome")

    # Create an integer ID for each genome by enumerating
    genome_ids = unique_genomes.with_row_index("genome_id")

    # Join back to the original gene_map, so each "genome" now has a "genome_id"
    gene_map_with_ids = (
        gene_map.join(genome_ids, on="genome", how="left")
    )
    return gene_map_with_ids

def assign_protein_ids(gene_map_with_ids: pl.DataFrame) -> pl.DataFrame:
    """
    Assigns integer IDs for each protein. This is especially useful if
    'protein' can be extremely long or repeated millions of times.
    Produces new columns ["protein_id"].
    """
    unique_ptns = gene_map_with_ids.select("protein").unique().sort("protein")
    ptn_ids = unique_ptns.with_row_index("protein_id")

    # Join to attach "protein_id" to each row
    gene_map_final = gene_map_with_ids.join(ptn_ids, on="protein", how="left")
    return gene_map_final

def replace_protein_names_in_alignments(alignment_file, protein_to_id, out_file, chunk_size=10_000_000):
    """
    Reads the alignment file in chunks (each containing chunk_size lines),
    replaces the protein names in the first two columns with their corresponding
    protein IDs using the protein_to_id dictionary, and writes the modified lines
    to out_file.

    Parameters:
      alignment_file: path to the input alignment file.
      protein_to_id: dictionary mapping protein names to protein IDs.
      out_file: path to write the transformed file.
      chunk_size: number of lines to process at once.
    """
    with open(alignment_file, "r") as infile, open(out_file, "w", newline='') as outfile:
        chunk = []
        for line in infile:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                processed_chunk = []
                for l in chunk:
                    parts = l.rstrip("\n").split("\t")
                    if len(parts) >= 2:
                        # Replace the first two columns with their protein IDs
                        parts[0] = str(protein_to_id.get(parts[0], parts[0]))
                        parts[1] = str(protein_to_id.get(parts[1], parts[1]))
                    processed_chunk.append("\t".join(parts))
                # Write the processed chunk at once
                outfile.write("\n".join(processed_chunk) + "\n")
                chunk = [] # Clear the chunk list
        # Process any remaining lines in the last chunk
        if chunk:
            processed_chunk = []
            for l in chunk:
                parts = l.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    parts[0] = str(protein_to_id.get(parts[0], parts[0]))
                    parts[1] = str(protein_to_id.get(parts[1], parts[1]))
                processed_chunk.append("\t".join(parts))
            outfile.write("\n".join(processed_chunk) + "\n")

def main():
    input_prot_subdir = snakemake.params.input_prot_subdir
    wdir = snakemake.params.wdir
    out_g2g = snakemake.params.gene_map_file
    alignment_file = snakemake.params.search_output
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    logger.info("Preparing inputs for AAI calculations...")
    
    if os.path.exists(os.path.join(input_prot_subdir, "vMAG_proteins")) and any(f.endswith('.faa') for f in os.listdir(os.path.join(input_prot_subdir, "vMAG_proteins")) if os.path.isfile
(os.path.join(os.path.join(input_prot_subdir, "vMAG_proteins"), f))):
        input_files = [os.path.join(input_prot_subdir, "single_contig_proteins.faa")] + [os.path.join(input_prot_subdir, "vMAG_proteins", faa) for faa in os.listdir(os.path.join(input_prot_subdir, "vMAG_proteins")) if faa.endswith(".faa")]
    else:
        input_files = [os.path.join(input_prot_subdir, "single_contig_proteins.faa")]
        
    gene_map = gene_to_genome(input_files).sort(["genome", "scaffold", "gene_number"])
    gene_map = assign_genome_ids(gene_map).sort(["genome", "scaffold", "gene_number"])
    gene_map = assign_protein_ids(gene_map).sort("protein_id").select(["protein", "genome", "protein_id", "genome_id"])
    
    logger.debug(f"Gene-to-genome mapping: {gene_map}")
    logger.debug(f"Number of genomes: {gene_map.select('genome').n_unique()}")
    logger.debug(f"Number of proteins: {gene_map.select('protein').n_unique()}")
    
    os.makedirs(wdir, exist_ok=True)
    gene_map.write_csv(out_g2g, separator="\t")
    logger.info("Gene-to-genome mapping completed and written to disk.")
    
    # Convert gene_map to a dictionary for quick lookup and free memory.
    protein_to_id = dict(zip(gene_map["protein"].to_list(), gene_map["protein_id"].to_list()))
    del gene_map
    gc.collect()
    
    logger.info("Replacing protein names with protein ids in the alignment file...")
    alignment_out = str(alignment_file).replace(".tsv", ".processed.tsv")
    replace_protein_names_in_alignments(alignment_file, protein_to_id, alignment_out)
    try:
        assert (os.path.exists(alignment_out) and os.path.getsize(alignment_out) > 0), f"Alignment file not found or is empty at {alignment_out}."
    except AssertionError as e:
        logger.error(e)
        raise e
    
    logger.info(f"Removing the original alignment file at {alignment_file}.")
    os.remove(alignment_file)
    
    logger.info("Preparing for AAI calculations completed.")
    
if __name__ == "__main__":
    main()