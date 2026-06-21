#!/usr/bin/env python3

import os
from pathlib import Path
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pyfastatools import Parser

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Map gene- and genome-level metadata", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

def parse_faa_file(faa_file, is_bin=False):
    """
    Parses an amino-acid fasta file to obtain protein names and their genomes.
    """ 
    try:
        for header in Parser(faa_file).all_prodigal_headers():
            seq_id, scaffold, gene_number, start, stop, frame = str(header.name()), str(header.scaffold), str(header.id), str(header.start), str(header.end), str(header.strand.value)
            if is_bin:
                genome = os.path.splitext(os.path.basename(faa_file))[0]
            else:
                genome = scaffold
            yield genome, scaffold, seq_id, gene_number, start, stop, frame
    except RuntimeError as e:
        logger.error(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
        raise RuntimeError(f"Error parsing .faa file {faa_file}: {e}. Headers are not in prodigal format.")
    
def process_data(single_contig_prots, bin_prots):
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
        "is_bin": []
    }

    if bin_prots:
        for file in bin_prots:
            for genome, scaffold, seq_id, gene_number, start, stop, frame in parse_faa_file(file, is_bin=True):
                faa_data["genome"].append(genome)
                faa_data["contig"].append(scaffold)
                faa_data["protein"].append(seq_id)
                faa_data["gene_number"].append(int(gene_number))
                faa_data["contig_pos_start"].append(int(start))
                faa_data["contig_pos_end"].append(int(stop))
                faa_data["frame"].append(int(frame))
                faa_data["is_bin"].append("true")
    
    if single_contig_prots:
        for genome, scaffold, seq_id, gene_number, start, stop, frame in parse_faa_file(single_contig_prots, is_bin=False):
            faa_data["genome"].append(genome)
            faa_data["contig"].append(scaffold)
            faa_data["protein"].append(seq_id)
            faa_data["gene_number"].append(int(gene_number))
            faa_data["contig_pos_start"].append(int(start))
            faa_data["contig_pos_end"].append(int(stop))
            faa_data["frame"].append(int(frame))
            faa_data["is_bin"].append("false")

    # Convert the list of dictionaries (faa_data) to a Polars DataFrame
    faa_dataframe = pl.DataFrame(faa_data)
    
    return faa_dataframe

def main():
    out_path = snakemake.params.gene_index
    out_parent = snakemake.params.out_parent
    save_to_parquet = snakemake.params.save_to_parquet
    
    protein_dir = snakemake.params.protein_dir
    bin_proteins_subdir = snakemake.params.bin_proteins_subdir
    
    if not os.path.exists(out_parent):
        os.makedirs(out_parent, exist_ok=True)
    
    if os.path.exists(os.path.join(protein_dir, 'single_contig_proteins.faa')):
        single_contig_prots = os.path.join(protein_dir, 'single_contig_proteins.faa')
    else:
        single_contig_prots = None
    
    if os.path.exists(bin_proteins_subdir) and os.path.isdir(bin_proteins_subdir):
        bin_prots = [os.path.join(bin_proteins_subdir, f) for f in os.listdir(bin_proteins_subdir) if f.endswith('.faa')]
    else:
        bin_prots = None
    
    logger.info("Starting genome and gene data mapping...")
    processed_data = process_data(single_contig_prots, bin_prots)
    
    if save_to_parquet:
        out_path_parquet = as_parquet_path_if_enabled(out_path, save_to_parquet)
        processed_data.write_parquet(out_path_parquet)
    else:
        processed_data.write_csv(out_path, separator="\t")
    logger.info("Genome and gene data mapping completed.")

if __name__ == "__main__":
    main()
