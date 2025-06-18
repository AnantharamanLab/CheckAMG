#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
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

print("========================================================================\n      Step 11/11: Write the final summarized auxiliary gene table       \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n      Step 11/11: Writing the final summarized auxiliary gene table       \n========================================================================\n")

# Function to determine Viral Origin Confidence
def viral_origin_confidence(circular, viral_window, viral_flank_up, viral_flank_down, mge_flank):
    # This can certainly be simplified and c4ompacted, but for clarity, I like it as is
    confidence_score = 0
    
    # 1) Being flanked by viral genes raises confidence
    # 1a) If both viral flanks are present, confidence is raised
    # 1b) If only one viral flank is present but the contig is circular, confidence is still raised
    # 1c) If neither viral flank is present, confidence is lowered
    if viral_flank_up:
        confidence_score += 1
    if viral_flank_down:
        confidence_score += 1
    if circular:
        confidence_score += 1
    # 2) Being in a viral window raises confidence,
    #    Not being in a viral window lowers confidence
    if viral_window:
        confidence_score += 1
    # 3) Being flanked by MGE genes lowers confidence,
    #    not being flanked by MGE genes raises confidence
    if not mge_flank:
        confidence_score += 1
    
    if confidence_score >= 4:
        return "high"
    elif confidence_score < 4 and confidence_score > 1:
        return "medium"
    elif confidence_score <= 1:
        return "low"
    else:
        logger.error(f"Unexpected confidence score: {confidence_score}. This should not happen.")
        raise ValueError(f"Unexpected confidence score: {confidence_score}. This should not happen.")
    
# Function to classify proteins based on their presence in metabolic, physiology, and regulatory tables
def classify_proteins(final_df, metabolism_df, physiology_df, regulatory_df):
    logger.debug(f"Columns in final_df before classification: {final_df.columns}")
    logger.debug(f"Columns in metabolism_df: {metabolism_df.columns}")
    logger.debug(f"Columns in physiology_df: {physiology_df.columns}")
    logger.debug(f"Columns in regulatory_df: {regulatory_df.columns}")
    
    # Ensure column names are consistent
    required_col = "Protein"
    for df_name, df in zip(["metabolism_df", "physiology_df", "regulatory_df"], [metabolism_df, physiology_df, regulatory_df]):
        if required_col not in df.columns:
            logger.error(f"Column '{required_col}' not found in {df_name}")
            raise ValueError(f"Column '{required_col}' not found in {df_name}")

    # Extract unique protein lists for fast lookup
    metabolic_proteins = set(metabolism_df["Protein"].to_list())
    physiology_proteins = set(physiology_df["Protein"].to_list())
    regulatory_proteins = set(regulatory_df["Protein"].to_list())

    # Assign classifications
    def classify(protein):
        if protein in metabolic_proteins:
            return "metabolic"
        elif protein in physiology_proteins:
            return "physiological"
        elif protein in regulatory_proteins:
            return "regulatory"
        else:
            return "unclassified"

    final_df = final_df.with_columns(
        pl.col("Protein").map_elements(classify, return_dtype=pl.Utf8).alias("classification"),
        pl.struct(["Circular_Contig", "Virus_Like_Window", "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes"])
        .map_elements(lambda x: viral_origin_confidence(x["Circular_Contig"], x["Virus_Like_Window"], x["Viral_Flanking_Genes_Upstream"], x["Viral_Flanking_Genes_Downstream"], x["MGE_Flanking_Genes"]),
                    return_dtype=pl.Utf8)
        .alias("Confidence")
    ).select(
        [
            "Protein",
            "Contig",
            "Genome",
            "classification",
            "Confidence",
            # "Circular_Contig", "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes", ## Debugging
            "KEGG_hmm_id",
            "KEGG_Description",
            "FOAM_hmm_id",
            "FOAM_Description",
            "Pfam_hmm_id",
            "Pfam_Description",
            "dbCAN_hmm_id",
            "dbCAN_Description",
            "METABOLIC_hmm_id",
            "METABOLIC_Description",
            "PHROG_hmm_id",
            "PHROG_Description",
            "top_hit_hmm_id",
            "top_hit_description",
            "top_hit_db",
        ]
    ).rename(
        {
            "classification": "Protein Classification",
            "Confidence": "Protein Viral Origin Confidence",
            "KEGG_hmm_id": "KEGG KO",
            "KEGG_Description": "KEGG KO Name",
            "FOAM_hmm_id": "FOAM ID",
            "FOAM_Description": "FOAM Annotation",
            "Pfam_hmm_id": "Pfam Accession",
            "Pfam_Description": "Pfam Name",
            "dbCAN_hmm_id": "CAZy Family",
            "dbCAN_Description": "CAZy Activities",
            "METABOLIC_hmm_id": "METABOLIC db ID",
            "METABOLIC_Description": "METABOLIC Annotation",
            "PHROG_hmm_id": "PHROG Number",
            "PHROG_Description": "PHROG Annotation",
            "top_hit_hmm_id": "Best Scoring HMM",
            "top_hit_description": "Best Scoring HMM Annotation",
            "top_hit_db": "Best Scoring HMM Origin",
        }
    )
    # Free memory
    del metabolism_df, physiology_df, regulatory_df
    gc.collect()

    return final_df

# Function to join the auxiliary status, hmm dataframe, all genes dataframe, and classification for CheckAMG annotate mode
def merge_dataframes(all_genes_df, metabolism_df, physiology_df, regulatory_df):
    logger.debug(f"Columns in all_genes_df: {all_genes_df.columns}")
    logger.debug(f"all_genes_df shape: {all_genes_df.shape}")
    return classify_proteins(all_genes_df, metabolism_df, physiology_df, regulatory_df)

def main():
    all_genes_path = snakemake.params.all_genes_annotated
    gene_index_path = snakemake.params.gene_index
    metabolism_path = snakemake.params.metabolism_table
    physiology_path = snakemake.params.physiology_table
    regulatory_path = snakemake.params.regulation_table
    final_table_path = snakemake.params.final_table
    mem_limit = snakemake.resources.mem
    threads = snakemake.threads
    set_memory_limit(mem_limit)
        
    all_genes_df = pl.read_csv(all_genes_path, separator='\t')
    gene_index_df = pl.read_csv(gene_index_path, separator='\t')
    gene_index_df = gene_index_df.select(["protein"] + [
        col for col in gene_index_df.columns if col.endswith("_protein_cluster") or col.endswith("_genome_cluster")
        ]).rename({"protein": "Protein"})
    metabolism_df = pl.read_csv(metabolism_path, separator='\t')
    physiology_df = pl.read_csv(physiology_path, separator='\t')
    regulatory_df = pl.read_csv(regulatory_path, separator='\t')
    
    logger.info(f"Generating the final table with proteins, annotations, and classifications...")
    
    final_df = merge_dataframes(all_genes_df, metabolism_df, physiology_df, regulatory_df)
    
    # Sort and write the final table output
    confidence_order = {"high": 0, "medium": 1, "low": 2}
    final_df = final_df.with_columns(
        pl.col("Protein Viral Origin Confidence").replace(confidence_order).cast(pl.Int32).alias("Protein Viral Origin Confidence_sort")
    )
    classification_order = {"metabolic": 0, "physiological": 1, "regulatory": 2, "unclassified": 3}
    final_df = final_df.with_columns(
        pl.col("Protein Classification").replace(classification_order).cast(pl.Int32).alias("Protein Classification_sort")
    ).sort(["Protein Classification_sort", "Protein Viral Origin Confidence_sort", "Protein"]).drop(["Protein Classification_sort", "Protein Viral Origin Confidence_sort"])
    final_df.write_csv(final_table_path, separator='\t')
    logger.info(f"Final table written to {final_table_path}")

    # Log classification summary
    logger.debug(f"Columns in final_df after classification: {final_df.columns}")
    logger.debug(f"Protein Classification value counts:\n{final_df['Protein Classification'].value_counts()}")
    logger.debug(f"Protein Viral Origin Confidence value counts:\n{final_df['Protein Viral Origin Confidence'].value_counts()}")
    
if __name__ == "__main__":
    main()