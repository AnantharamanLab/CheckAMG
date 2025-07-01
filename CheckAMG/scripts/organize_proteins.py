#!/usr/bin/env python3

import os
import sys
import load_prot_paths
from pyfastatools import Parser, write_fasta
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from concurrent.futures import ThreadPoolExecutor, as_completed

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
        logging.FileHandler(log_file, mode='a'),  # Append mode for log file
        logging.StreamHandler(sys.stdout)  # Print to console
    ]
)
logger = logging.getLogger()

print("========================================================================\n   Step 10/11: Organize proteins into putative AMGs, APGs, and AReGs    \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n   Step 10/11: Organize proteins into putative AMGs, APGs, and AReGs    \n========================================================================\n")

def load_protein_sequences(prot_paths):
    """Loads protein sequences from provided FASTA files into a dictionary."""
    prot_records = {}
    for prot_path in prot_paths:
        for record in Parser(prot_path):
            prot_records[record.header.name] = record
    logger.debug(f"Loaded {len(prot_records):,} protein sequences.")
    return prot_records

def adjust_sequence_header(record, protein):
    """Adjusts the header of a sequence record based on available annotations."""
    if 'top_hit_hmm_id' in protein and 'top_hit_description' in protein:
        if protein['top_hit_db'] != "PHROG":
            # Format the header without brackets
            hmm_id = protein['top_hit_hmm_id'][0] if isinstance(protein['top_hit_hmm_id'], list) else protein['top_hit_hmm_id']
            description = protein['top_hit_description'][0] if isinstance(protein['top_hit_description'], list) else protein['top_hit_description']
            record.header.desc = f'"{hmm_id} {description}"' if hmm_id and description else f'"{hmm_id or description}"'
        else:
            record.header.desc = determine_phrog_description(protein)
    record.header.name = record.header.name.split(" ")[0] # Keep only the sequence name

def determine_phrog_description(protein):
    """Determines PHROG-based description for a record."""
    db_fields = [('dbCAN_hmm_id', 'dbCAN_Description'), ('KEGG_hmm_id', 'KEGG_Description'), ('FOAM_hmm_id', 'FOAM_Description'), ('METABOLIC_hmm_id', 'METABOLIC_Description'), ('Pfam_hmm_id', 'Pfam_Description')]
    for hmm_id, desc in db_fields:
        if protein[hmm_id] and protein[desc]:
            return f'"{protein[hmm_id]} {protein[desc]}"'
    return ''

def write_fasta_str(record):
    return f">{record.header.name} {record.header.desc}\n{str(record.seq)}\n"

def organize_proteins(category_table_path, category, all_genes_df):
    """
    Organizes proteins based on annotations.
    This function reads metabolic/regulatory/physiology gene data and gene annotations from specified file paths.
    It then categorizes the genes into three broad sets based on their characteristics: 
    1. avgs_high: High-confidence AVGs (classified as "high" by `viral_origin_confidence`).
    2. avgs_medium: Medium-confidence AVGs (classified as "medium" by `viral_origin_confidence`).
    3. avgs_low: Low-confidence AVGs (classified as "low" by `viral_origin_confidence`).
    4. avgs_all: Any AVG classified as an auxiliary metabolic/regulatory/physiology gene, regardless of genomic context.
    
    Parameters:
    category_table_path (str): The file path to the metabolic/regulatory/physiology genes curated table in TSV format.
    
    Returns:
    dict: A dictionary containing categorized auxiliary genes and AMGs/AReGs/APGs.
    """
    
    category_genes_df = pl.read_csv(category_table_path, separator='\t')
    all_category_genes = set(category_genes_df['Protein'].to_list())
    acronym = {"metabolic": "AMG", "regulatory": "AReG", "physiology": "APG"}.get(category, "Unknown")
    logger.info(f"There are a total of {len(all_category_genes):,} putative {acronym}s.")
    
    # Extract relevant gene context information
    all_genes_info = (
        all_genes_df
        .select(["Protein", "Circular_Contig", "Viral_Origin_Confidence", "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes"])
        .to_dicts()
    )
    
    # Classify proteins using viral_origin_confidence function
    avgs_high = set()
    avgs_medium = set()
    avgs_low = set()

    for gene in all_genes_info:
        protein = gene["Protein"]
        confidence = gene["Viral_Origin_Confidence"]

        if protein in all_category_genes:
            if confidence == "high":
                avgs_high.add(protein)
            elif confidence == "medium":
                avgs_medium.add(protein)
            else:
                avgs_low.add(protein)

    logger.info(f"There are {len(avgs_high):,} {acronym}s classified as HIGH confidence viral origin.")
    logger.info(f"There are {len(avgs_medium):,} {acronym}s classified as MEDIUM confidence viral origin.")
    logger.info(f"There are {len(avgs_low):,} {acronym}s classified as LOW confidence viral origin.")

    return {
        "avgs_high": avgs_high, 
        "avgs_medium": avgs_medium, 
        "avgs_low": avgs_low, 
        "avgs_all": all_category_genes
    }

def write_organized_files(organized_dict, category, category_table_path, prot_records, output_dir):
    """
    Writes organized protein sequences to separate FASTA files based on the provided categories.
    
    Parameters:
    organized_dict (dict): A dictionary where keys are category names and values are sets of protein names.
    prot_records (dict): A dictionary of protein records where keys are protein names and values are sequence strings.
    output_dir (str): The directory where the output FASTA files will be saved.
    
    Returns:
    None
    """
    
    # Read category annotation table and convert to lookup dictionary
    category_genes_df = pl.read_csv(category_table_path, separator='\t')
    category_genes_lookup = {
        row["Protein"]: row for row in category_genes_df.iter_rows(named=True)
    }

    # Get acronym prefix for output files
    acronym = {"metabolic": "AMG", "regulatory": "AReG", "physiology": "APG"}.get(category, "Unknown")
    
    # Define output filenames by confidence level
    filename_dict = {
        "avgs_high" : f"{acronym}s_high_confidence.faa",
        "avgs_medium" : f"{acronym}s_medium_confidence.faa",
        "avgs_low" : f"{acronym}s_low_confidence.faa",
        "avgs_all" : f"{acronym}s_all.faa"
    }
    
    # Create output directory if it doesn't exist
    output_subdir = os.path.join(output_dir, f"faa_{category}")
    os.makedirs(output_subdir, exist_ok=True)

    def write_fasta_file(key, protein_names):
        output_fasta = os.path.join(output_subdir, filename_dict[key])
        fasta_lines = []

        for protein_name in protein_names:
            if protein_name in prot_records:
                protein_row = category_genes_lookup.get(protein_name)
                if protein_row:
                    adjust_sequence_header(prot_records[protein_name], protein_row)
                fasta_lines.append(write_fasta_str(prot_records[protein_name]))
            else:
                logger.warning(f"Protein '{protein_name}' not found in the loaded protein records.")

        with open(output_fasta, "w") as out_f:
            out_f.write("".join(fasta_lines))
    
    # Write each confidence level file in parallel
    with ThreadPoolExecutor() as executor:
        futures = [
            executor.submit(write_fasta_file, key, protein_names)
            for key, protein_names in organized_dict.items()
        ]
        for f in as_completed(futures):
            f.result()

def main():
    metab_table_path = snakemake.params.metabolism_table
    phys_table_path = snakemake.params.physiology_table
    reg_table_path = snakemake.params.regulation_table
    all_genes_path = snakemake.params.all_genes_annotated
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)    

    logger.info("Organizing proteins based on annotations and writing AMG/AReG/APG classifications...")
    
    input_prots_dir = snakemake.params.protein_dir
    fasta_outdir = snakemake.params.aux_fasta_dir
    prot_paths = load_prot_paths.load_prots(input_prots_dir)
    prot_records = load_protein_sequences(prot_paths)
    
    all_genes_df = pl.read_csv(all_genes_path, separator='\t')
    
    for category in ["metabolic", "regulatory", "physiology"]:
        if category == "metabolic":
            organized_dict = organize_proteins(metab_table_path, category, all_genes_df)
            write_organized_files(organized_dict, category, metab_table_path, prot_records, fasta_outdir)
        elif category == "regulatory":
            organized_dict = organize_proteins(reg_table_path, category, all_genes_df)
            write_organized_files(organized_dict, category, reg_table_path, prot_records, fasta_outdir)
        elif category == "physiology":
            organized_dict = organize_proteins(phys_table_path, category, all_genes_df)
            write_organized_files(organized_dict, category, phys_table_path, prot_records, fasta_outdir)
            
    logger.debug(f"Results were written to {fasta_outdir}.")
    logger.info(f"Organization completed.")
    
if __name__ == "__main__":
    main()
