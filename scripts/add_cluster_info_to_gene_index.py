#!/usr/bin/env python3

import os
import sys
import resource
import platform
from pathlib import Path
import logging
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

print("========================================================================\n"
      "      Step 15/22: Map genome and protein clusters to protein info       \n"
      "========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n"
              "      Step 15/22: Map genome and protein clusters to protein info       \n"
              "========================================================================\n")

def merge_data(protein_clusters_tsv, gene_index_tsv):
    protein_clusters = pl.read_csv(protein_clusters_tsv, separator='\t', infer_schema_length=0)
    gene_index = pl.read_csv(gene_index_tsv, separator='\t', infer_schema_length=0)
    merged_df = protein_clusters.join(gene_index, on="protein", how="left")
    for col in merged_df.columns:
        if col.endswith("_right"):
            merged_df.drop_in_place(col)
    return merged_df
    
def main():
    gene_index_tsv = snakemake.params.gene_index_annotated
    protein_clusters_tsv = snakemake.params.protein_clusters_tsv
    out_parent = snakemake.params.out_parent
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    if not os.path.exists(out_parent):
        os.makedirs(out_parent, exist_ok=True)
    
    logger.info("Starting to add genome and protein cluster assignments to the annotated gene index...")
    processed_data = merge_data(protein_clusters_tsv, gene_index_tsv)
    processed_data.write_csv(gene_index_tsv, separator="\t")
    logger.info("Adding cluster assignments completed.")

if __name__ == "__main__":
    main()
