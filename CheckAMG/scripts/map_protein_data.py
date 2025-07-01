#!/usr/bin/env python3

import os
import sys
import resource
import platform
from pathlib import Path
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pyfastatools import Parser

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

print("========================================================================\n            Step 6/11: Map gene- and genome-level metadata             \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n            Step 6/11: Map gene- and genome-level metadata             \n========================================================================\n")

def parse_faa_file(faa_file, is_vMAG=False):
    """
    Parses an amino-acid fasta file to obtain protein names and their genomes.
    """ 
    try:
        for header in Parser(faa_file).all_prodigal_headers():
            seq_id, scaffold, gene_number, start, stop, frame = str(header.name()), str(header.scaffold), str(header.id), str(header.start), str(header.end), str(header.strand.value)
            if is_vMAG:
                genome = os.path.splitext(os.path.basename(faa_file))[0]
            else:
                genome = scaffold
            yield genome, scaffold, seq_id, gene_number, start, stop, frame
    except RuntimeError as e:
        logger.error(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
        raise RuntimeError(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
    
def process_data(single_contig_prots, vmag_prots, ranks):
    """
    Obtains gene-level data from prodigal-formatted input faa(s)
    
    Returns a Polars DataFrame with gene- and genome-level data.
    """

    # Accumulate all faa data in chunks to avoid large memory usage
    faa_data = {
        "genome": [],
        "contig": [],
        "protein": [],
        "gene_number": [],
        "contig_pos_start": [],
        "contig_pos_end": [],
        "frame": [],
        "is_vMAG": []
    }

    if vmag_prots:
        for file in vmag_prots:
            for genome, scaffold, seq_id, gene_number, start, stop, frame in parse_faa_file(file, is_vMAG=True):
                faa_data["genome"].append(genome)
                faa_data["contig"].append(scaffold)
                faa_data["protein"].append(seq_id)
                faa_data["gene_number"].append(int(gene_number))
                faa_data["contig_pos_start"].append(int(start))
                faa_data["contig_pos_end"].append(int(stop))
                faa_data["frame"].append(int(frame))
                faa_data["is_vMAG"].append("true")
    
    if single_contig_prots:
        for genome, scaffold, seq_id, gene_number, start, stop, frame in parse_faa_file(single_contig_prots, is_vMAG=False):
            faa_data["genome"].append(genome)
            faa_data["contig"].append(scaffold)
            faa_data["protein"].append(seq_id)
            faa_data["gene_number"].append(int(gene_number))
            faa_data["contig_pos_start"].append(int(start))
            faa_data["contig_pos_end"].append(int(stop))
            faa_data["frame"].append(int(frame))
            faa_data["is_vMAG"].append("false")

    # Convert the list of dictionaries (faa_data) to a Polars DataFrame
    faa_dataframe = pl.DataFrame(faa_data)
    
    return faa_dataframe

def main():
    all_genes = snakemake.params.gene_index
    ranks = snakemake.params.cluster_taxa_levels
    out_parent = snakemake.params.out_parent
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    protein_dir = snakemake.params.protein_dir
    vmag_proteins_subdir = snakemake.params.vmag_proteins_subdir
    
    if not os.path.exists(out_parent):
        os.makedirs(out_parent, exist_ok=True)
    
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        single_contig_prots = os.path.join(protein_dir, 'single_contig_proteins.faa')
    else:
        single_contig_prots = None
    
    if os.path.exists(vmag_proteins_subdir) and os.path.isdir(vmag_proteins_subdir):
        vmag_prots = [os.path.join(vmag_proteins_subdir, f) for f in os.listdir(vmag_proteins_subdir) if f.endswith('.faa')]
    else:
        vmag_prots = None
    
    logger.info("Starting genome and gene data mapping...")
    processed_data = process_data(single_contig_prots, vmag_prots, ranks)
    
    processed_data.write_csv(all_genes, separator="\t")
    logger.info("Genome and gene data mapping completed.")

if __name__ == "__main__":
    main()
