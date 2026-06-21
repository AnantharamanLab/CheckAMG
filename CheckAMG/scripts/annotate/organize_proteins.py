#!/usr/bin/env python3

import os
from pyfastatools import Parser, write_fasta
from pathlib import Path
import re
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from concurrent.futures import ThreadPoolExecutor, as_completed

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled
from CheckAMG.scripts.annotate import utils

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Organize AMGs, APGs, and AReGs proteins", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")


def load_protein_sequences(prot_paths):
    """Loads protein sequences from provided FASTA files into a dictionary."""
    prot_records = {}
    for prot_path in prot_paths:
        for record in Parser(prot_path):
            prot_records[record.header.name] = record
    logger.debug(f"Loaded {len(prot_records):,} protein sequences.")
    return prot_records

_EC_NUM_RE = re.compile(r'(?i)(?:\bEC[: ]\s*)?(?:\d+|\-)\.(?:\d+|\-)\.(?:\d+|\-)\.(?:\d+|\-)')

def _normalize_function_text(s):
    if s is None:
        return ""
    if isinstance(s, list):
        s = s[0] if s else ""
    s = str(s)
    s = s.replace("\t", " ").replace("\n", " ").replace("\r", " ")
    s = re.sub(r"\s+", " ", s).strip()
    s = s.strip('"').strip()
    s = re.sub(r"\s+", " ", s).strip()
    return s

def _truncate_function_text(s, max_len=150):
    s = _normalize_function_text(s)
    if not s:
        return ""
    if len(s) <= max_len:
        return s

    extra = 12 # allow a small spillover to avoid broken delimiters

    def _extend_tail(end_idx: int) -> int:
        i = end_idx
        # include immediate closers/punctuation if present right after the cut
        while i < len(s) and i < (max_len + extra):
            ch = s[i]
            if ch in ")]};,": # keep common closers and separators
                i += 1
                continue
            break
        return i

    best_end = None
    for m in _EC_NUM_RE.finditer(s):
        end = m.end()
        if end <= max_len:
            best_end = end

    if best_end is not None and best_end > 0:
        end2 = _extend_tail(best_end)
        out = s[:end2].rstrip(" ")
        return out.rstrip(" ;,")

    # fallback: cut on last space, then try to close an unclosed '(' or '[' if nearby
    cut = s.rfind(" ", 0, max_len + 1)
    if cut <= 0:
        cut = max_len
    out = s[:cut].rstrip(" ")

    # avoid leaving a dangling opener if the closer is very close
    if out.count("(") > out.count(")"):
        j = s.find(")", cut)
        if j != -1 and j < (max_len + extra):
            out = s[:j + 1].rstrip(" ")
    if out.count("[") > out.count("]"):
        j = s.find("]", cut)
        if j != -1 and j < (max_len + extra):
            out = s[:j + 1].rstrip(" ")

    out = out.rstrip(" ;,")
    return out

def adjust_sequence_header(record, protein):
    """Adjusts the header string fields based on available annotations."""
    seq_name = record.header.name.split(" ")[0]  # Keep only the sequence name
    func = protein.get("Function", "")
    func = _truncate_function_text(func, max_len=150)
    desc = f'"{func}"' if func else ""
    return seq_name, desc

def write_fasta_str(name, desc, seq):
    if desc:
        return f">{name} {desc}\n{seq}\n"
    return f">{name}\n{seq}\n"

def organize_proteins(category_table_path, category, all_genes_df, save_to_parquet=False):
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

    if save_to_parquet:
        category_genes_df = pl.read_parquet(category_table_path)
    else:
        category_genes_df = pl.read_csv(category_table_path, separator='\t')
    all_category_genes = set(category_genes_df['Protein'].to_list())
    acronym = {"metabolic": "AMG", "regulatory": "AReG", "physiology": "APG"}.get(category, "Unknown")
    logger.info(f"There are a total of {len(all_category_genes):,} predicted {acronym}s.")

    # Extract relevant gene context information
    all_genes_info = (
        all_genes_df
        .select(["Protein", "Viral_Origin_Confidence", "Viral_Flanking_Genes_Left_Dist", "Viral_Flanking_Genes_Right_Dist", "MGE_Flanking_Genes_Left_Dist", "MGE_Flanking_Genes_Right_Dist"])
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

def write_organized_files(organized_dict, category, category_table_path, prot_records, output_dir, save_to_parquet=False):
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
    if save_to_parquet:
        category_genes_df = pl.read_parquet(category_table_path)
    else:
        category_genes_df = pl.read_csv(category_table_path, separator='\t')
    category_genes_lookup = {
        row["Protein"]: row for row in category_genes_df.iter_rows(named=True)
    }

    # Get acronym prefix for output files
    acronym = {"metabolic": "AMG", "regulatory": "AReG", "physiology": "APG"}.get(category, "Unknown")

    # Define output filenames by confidence level
    filename_dict = {
        "avgs_high": f"{acronym}s_high_confidence.faa",
        "avgs_medium": f"{acronym}s_medium_confidence.faa",
        "avgs_low": f"{acronym}s_low_confidence.faa",
        "avgs_all": f"{acronym}s_all.faa"
    }

    # Create output directory if it doesn't exist
    output_subdir = os.path.join(output_dir, f"faa_{category}")
    os.makedirs(output_subdir, exist_ok=True)

    def write_fasta_file(key, protein_names):
        output_fasta = os.path.join(output_subdir, filename_dict[key])
        fasta_lines = []

        for protein_name in protein_names:
            record = prot_records.get(protein_name)
            if record is None:
                logger.warning(f"Protein '{protein_name}' not found in the loaded protein records.")
                continue

            protein_row = category_genes_lookup.get(protein_name, {})
            name, desc = adjust_sequence_header(record, protein_row)
            fasta_lines.append(write_fasta_str(name, desc, str(record.seq)))

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
    save_to_parquet = snakemake.params.save_to_parquet
    metab_table_path = as_parquet_path_if_enabled(snakemake.params.metabolism_table, save_to_parquet)
    phys_table_path = as_parquet_path_if_enabled(snakemake.params.physiology_table, save_to_parquet)
    reg_table_path = as_parquet_path_if_enabled(snakemake.params.regulation_table, save_to_parquet)
    all_genes_path = as_parquet_path_if_enabled(snakemake.params.all_genes_annotated, save_to_parquet)
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Organizing proteins based on annotations and writing AMG/AReG/APG classifications...")

    input_prots_dir = snakemake.params.protein_dir
    bin_proteins_subdir = snakemake.params.bin_proteins_subdir
    fasta_outdir = snakemake.params.aux_fasta_dir
    prot_paths = utils.load_prots(input_prots_dir, bin_proteins_subdir)
    prot_records = load_protein_sequences(prot_paths)

    if save_to_parquet:
        all_genes_df = pl.read_parquet(all_genes_path)
    else:
        all_genes_df = pl.read_csv(all_genes_path, separator='\t')

    for category in ["metabolic", "physiology", "regulatory"]:
        if category == "metabolic":
            organized_dict = organize_proteins(metab_table_path, category, all_genes_df, save_to_parquet)
            write_organized_files(organized_dict, category, metab_table_path, prot_records, fasta_outdir, save_to_parquet)
        elif category == "regulatory":
            organized_dict = organize_proteins(reg_table_path, category, all_genes_df, save_to_parquet)
            write_organized_files(organized_dict, category, reg_table_path, prot_records, fasta_outdir, save_to_parquet)
        elif category == "physiology":
            organized_dict = organize_proteins(phys_table_path, category, all_genes_df, save_to_parquet)
            write_organized_files(organized_dict, category, phys_table_path, prot_records, fasta_outdir, save_to_parquet)

    logger.debug(f"Results were written to {fasta_outdir}.")
    logger.info(f"Organization completed.")

if __name__ == "__main__":
    main()
