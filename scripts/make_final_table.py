#!/usr/bin/env python3

import os
import sys
import logging
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import duckdb
import re
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

if snakemake.params.build_or_annotate =="build":
    print("========================================================================\n     Step 22/22: Writing the final summarized auxiliary gene table      \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n     Step 22/22: Writing the final summarized auxiliary gene table      \n========================================================================\n")
elif snakemake.params.build_or_annotate == "annotate":
    print("========================================================================\n     Step 11/11: Writing the final summarized auxiliary gene table      \n========================================================================")
    with open(log_file, "a") as log:
        log.write("========================================================================\n     Step 11/11: Writing the final summarized auxiliary gene table      \n========================================================================\n")

def build_protein_metrics(aux_scores_df: pl.DataFrame, hmm_df: pl.DataFrame) -> pl.DataFrame:
    """
    Builds a wide table (one row per protein) with, for each rank:
      - {rank}-Level Quantile (from quantile_weighted),
      - {rank}-Level Auxiliary Score,
      - {rank} Designation (derived from the quantile).
    The designation is made as follows (assuming quantile_weighted is already an integer 1-10):
        1           -> "strongly auxiliary"
        2–3         -> "weakly auxiliary"
        4–7         -> "conditional"
        8–9         -> "weakly core"
        10          -> "strongly core"
    Then the table is joined to hmm_df (on protein/Protein) so that each protein row includes its HMM.
    """
    rank_dfs = []
    
    for rank in aux_scores_df["rank"].unique().to_list():
        rank_df = (
            aux_scores_df
            .filter(pl.col("rank") == rank)
            .select(["protein", "quantile_weighted", "aux_score_weighted"])
            .with_columns(
                pl.when(pl.col("quantile_weighted") == 1).then(pl.lit("strongly auxiliary"))
                  .when(pl.col("quantile_weighted").is_in([2, 3])).then(pl.lit("weakly auxiliary"))
                  .when(pl.col("quantile_weighted").is_in([4, 5, 6, 7])).then(pl.lit("conditional"))
                  .when(pl.col("quantile_weighted").is_in([8, 9])).then(pl.lit("weakly core"))
                  .when(pl.col("quantile_weighted") == 10).then(pl.lit("strongly core"))
                  .otherwise(pl.lit("conditional"))
                  .alias("Designation")
            )
            .rename({
                "quantile_weighted": f"{rank}-Level Quantile",
                "aux_score_weighted": f"{rank}-Level Auxiliary Score"
            })
        )
        rank_dfs.append(rank_df)
    
    merged_df = rank_dfs[0]
    for i in range(1, len(rank_dfs)):
        merged_df = merged_df.join(rank_dfs[i], on="protein", how="left")
    
    # Join with HMM info. (Assumes hmm_df has column "Protein" and "hmm")
    merged_df = merged_df.join(hmm_df, left_on="protein", right_on="Protein", how="right")
    return merged_df

def build_family_aux_status(aux_scores_df: pl.DataFrame, hmm_df: pl.DataFrame) -> pl.DataFrame:
    """
    Aggregates at the HMM (family) level.
    For each (protein, rank), we first compute an integer quantile score (q_int) 
    (i.e. rounded quantile_weighted). Then, joining with hmm_df so that each protein
    is linked to its hmm, we group by hmm and compute the average q_int across all
    occurrences (across ranks). Finally, we assign a family-level final designation using:
        avg  [1, 2)   -> "strongly auxiliary"
        avg  [2, 4)   -> "weakly auxiliary"
        avg  [4, 7)   -> "conditional"
        avg  [7, 9)   -> "weakly core"
        else [9, 10]  -> "strongly core"
    """
    df_flag = aux_scores_df.with_columns(
        pl.col("quantile_weighted").round(0).cast(pl.Int32).alias("q_int")
    )
    
    # Join with HMM table so that each protein gets its hmm
    joined = (
        df_flag.select(["protein", "rank", "q_int"])
        .join(
            hmm_df.select(["Protein", "hmm"]),
            left_on="protein", right_on="Protein", how="inner"
        )
    )
    
    hmm_aggregated = (
        joined
        .group_by("hmm")
        .agg(
            pl.mean("q_int").alias("avg_quantile")
        )
    )
    
    family_status = (
        hmm_aggregated
        .with_columns(
            pl.when(pl.col("avg_quantile") < 2).then(pl.lit("strongly auxiliary"))
              .when(pl.col("avg_quantile") < 4).then(pl.lit("weakly auxiliary"))
              .when(pl.col("avg_quantile") < 7).then(pl.lit("conditional"))
              .when(pl.col("avg_quantile") < 9).then(pl.lit("weakly core"))
              .otherwise(pl.lit("strongly core"))
              .alias("auxiliary_status")
        )
        .select(["hmm", "auxiliary_status", "avg_quantile"])
    )
    return family_status

# Function to determine Viral Origin Confidence
def viral_origin_confidence(circular, viral_window, viral_flank_up, viral_flank_down, mge_flank):
    # This can certainly be simplified and compacted, but for clarity, I like it as is
    confidence_score = 0
    
    # 1) Being flanked by viral genes raises confidence
    # 1a) If both viral flanks are present, confidence is raised
    # 1b) If only one viral flank is present but the contig is circular, confidence is still raised
    # 1c) If neither viral flank is present, confidence is lowered
    if viral_flank_up and viral_flank_down:
        confidence_score += 1
    elif (viral_flank_up or viral_flank_down) and circular:
        confidence_score += 1
    else:
        confidence_score += 0
    # 2) Being in a viral window raises confidence,
    #    Not being in a viral window lowers confidence
    if viral_window:
        confidence_score += 1
    else:
        confidence_score += 0
    # 3) Being flanked by MGE genes lowers confidence,
    #    not being flanked by MGE genes raises confidence
    if not mge_flank:
        confidence_score += 1
    else:
        confidence_score += 0
    
    if confidence_score == 3:
        return "high"
    elif confidence_score == 2:
        return "medium"
    elif confidence_score <= 1:
        return "low"
    else:
        logger.error(f"Unexpected confidence score: {confidence_score}. This should not happen.")
        raise ValueError(f"Unexpected confidence score: {confidence_score}. This should not happen.")
    
# Function to classify proteins based on their presence in metabolic, physiology, and regulatory tables
def classify_proteins(final_df, metabolism_df, physiology_df, regulatory_df, build_or_annotate):
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

    # Apply classification
    if build_or_annotate == "build":
        membership_df = final_df.group_by("hmm").agg(pl.n_unique("Protein").alias("hmm_membership"))
        final_df = final_df.join(membership_df, on="hmm", how="left")
        del membership_df
        
        rank_cols = [col for col in final_df.columns if col.endswith("_genome_cluster")] +\
                    [col for col in final_df.columns if col.endswith("_protein_cluster")] +\
                    [col for col in final_df.columns if col.endswith("Auxiliary Score")] +\
                    [col for col in final_df.columns if col.endswith("Quantile")]
        rank_cols = sorted(rank_cols)
        selected_cols = [
                "Protein",
                "Contig",
                "Genome",
                "hmm",
                "hmm_membership",
                "auxiliary_status",
                "classification",
                "Confidence",
                # "Circular_Contig", "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes", ## Debugging
        ] + rank_cols + [
                "KEGG_hmm_id",
                "KEGG_Description",
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
        
        final_df = final_df.with_columns(
            pl.col("Protein").map_elements(classify, return_dtype=pl.Utf8).alias("classification"),
            pl.struct(["Circular_Contig", "Virus_Like_Window", "Viral_Flanking_Genes_Upstream", "Viral_Flanking_Genes_Downstream", "MGE_Flanking_Genes"])
            .map_elements(lambda x: viral_origin_confidence(x["Circular_Contig"], x["Virus_Like_Window"], x["Viral_Flanking_Genes_Upstream"], x["Viral_Flanking_Genes_Downstream"], x["MGE_Flanking_Genes"]),
                        return_dtype=pl.Utf8)
            .alias("Confidence")
        ).select(
            selected_cols
        ).rename(
            {
                "hmm": "Family",
                "hmm_membership": "Size",
                "auxiliary_status": "Auxiliary Status",
                "classification": "Protein Classification",
                "Confidence": "Protein Viral Origin Confidence",
                "KEGG_hmm_id": "KEGG KO",
                "KEGG_Description": "KEGG KO Name",
                "Pfam_hmm_id": "Pfam Accession",
                "Pfam_Description": "Pfam Name",
                "dbCAN_hmm_id": "CAZy Family",
                "dbCAN_Description": "CAZy Activities",
                "METABOLIC_hmm_id": "METABOLIC db ID",
                "METABOLIC_Description": "METABOLIC Annotation",
                "PHROG_hmm_id": "PHROG Number",
                "PHROG_Description": "PHROG Annotation",
                "top_hit_hmm_id": "Top Hit HMM",
                "top_hit_description": "Top Hit Annotation",
                "top_hit_db": "Top Hit HMM Origin",
            }
        )
        
        rank_rename = {col: col.replace("class", "Protein Class").replace("family", "Protein Family").replace("genus", "Protein Genus").replace("species", "Protein Species").replace("_protein_cluster", " Protein Cluster").replace("_genome_cluster", " Genome Cluster") for col in rank_cols}
        final_df = final_df.rename(rank_rename)
        
    elif build_or_annotate == "annotate":
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
    else:
        logging.error(f"Invalid option for build_or_annotate: {build_or_annotate}")
        raise ValueError(f"Invalid option for build_or_annotate: {build_or_annotate}")

    # Free memory
    del metabolism_df, physiology_df, regulatory_df
    gc.collect()

    return final_df

# Function to join the auxiliary status, hmm dataframe, all genes dataframe, and classification for CheckAMG build mode
def merge_dataframes_build(aux_metric_df, hmm_df, all_genes_df, gene_index_df, metabolism_df, physiology_df, regulatory_df, aux_status_df):
    # Log available columns before joining
    logger.debug(f"Columns in aux_metric_df: {aux_metric_df.columns}")
    dups_in_am = (
        aux_metric_df
        .group_by("Protein")
        .agg(pl.count("Protein").alias("count"))
        .filter(pl.col("count") > 1)
    )
    logger.debug(f"Duplicates in aux_metric_df: {dups_in_am}")
    logger.debug(f"Columns in aux_status_df: {aux_status_df.columns}")
    dups_in_as = (
        aux_status_df
        .group_by("hmm")
        .agg(pl.count("hmm").alias("count"))
        .filter(pl.col("count") > 1)
    )
    logger.debug(f"Duplicates in aux_status_df: {dups_in_as}")
    logger.debug(f"Columns in hmm_df: {hmm_df.columns}")
    dups_in_hmm = (
        hmm_df
        .group_by("Protein")
        .agg(pl.count("Protein").alias("count"))
        .filter(pl.col("count") > 1)
    )
    logger.debug(f"Duplicate proteins in hmm_df: {dups_in_hmm}")
    logger.debug(f"Columns in all_genes_df: {all_genes_df.columns}")
    dups_in_ag = (
        all_genes_df
        .group_by("Protein")
        .agg(pl.count("Protein").alias("count"))
        .filter(pl.col("count") > 1)
    )
    logger.debug(f"Duplicates in all_genes_df: {dups_in_ag}")
    logger.debug(f"Columns in gene_index_df: {gene_index_df.columns}")
    dups_in_gi = (
        gene_index_df
        .group_by("Protein")
        .agg(pl.count("Protein").alias("count"))
        .filter(pl.col("count") > 1)
    )
    logger.debug(f"Duplicates in gene_index_df: {dups_in_gi}")

    # Join the auxiliary status and auxiliary metrics dfs
    merged_df = (
        aux_metric_df.lazy()
        .join(aux_status_df.lazy(), on="hmm", how="full")
        .lazy()
    )
    # Join with the other data
    merged_df = (
        merged_df
        .join(hmm_df.lazy(), on="hmm", how="left")
        .join(all_genes_df.lazy(), on="Protein", how="left")
        .join(gene_index_df.lazy(), on="Protein", how="left")
        .collect()
    )

    # Log result
    logger.debug(f"Merged dataframe shape: {merged_df.shape}")
    logger.debug(f"Columns in merged_df: {merged_df.columns}")
    
    merged_df = classify_proteins(merged_df, metabolism_df, physiology_df, regulatory_df, "build").unique()
    
    # Clean up memory
    del aux_metric_df, hmm_df, all_genes_df, gene_index_df, metabolism_df, physiology_df, regulatory_df
    gc.collect()

    return merged_df

# Function to join the auxiliary status, hmm dataframe, all genes dataframe, and classification for CheckAMG annotate mode
def merge_dataframes_genome(all_genes_df, metabolism_df, physiology_df, regulatory_df):
    logger.debug(f"Columns in all_genes_df: {all_genes_df.columns}")
    logger.debug(f"all_genes_df shape: {all_genes_df.shape}")
    return classify_proteins(all_genes_df, metabolism_df, physiology_df, regulatory_df, "annotate")

# Choose a consensus classification for a protein family
def consensus_classification(class_list: list[str]) -> str:
    """
    Return one of {"metabolic", "physiological", "regulatory", "unclassified"}
    following these rules:
    - if the family has size 1, then the lone protein's classification
    - otherwise, tally a weighted vote:
        - metabolic / physiological / regulatory = 1.0
        - unclassified (or null/""/NA)           = 0.4
    - class with highest score wins; ties, then first in priority list.
    """
    if len(class_list) == 1:
        return class_list[0]

    weights = {"metabolic": 1.0,
               "physiological": 1.0,
               "regulatory": 1.0,
               "unclassified": 0.4,
               None: 0.4, "": 0.4}
    priority = ["metabolic", "physiological", "regulatory", "unclassified"]

    score = {}
    for c in class_list:
        c_norm = c if c in weights else "unclassified"
        score[c_norm] = score.get(c_norm, 0) + weights[c_norm]

    best = max(score.items(), key=lambda kv: (kv[1], -priority.index(kv[0])))
    return best[0]

# Choose a consensus functional annotation for a protein family
def consensus_annotation(annot_score_list: list[dict]) -> str:
    """
    Determine the consensus functional annotation for a protein family 
    using weighted support from individual protein annotations.

    Rules:
    - Each non-empty annotation contributes support weighted by its bitscore.
    - Annotations are penalized (x0.4) if they contain:
        - the whole word "hypothetical" (case-insensitive),
        - the whole word "unknown" (case-insensitive),
        - or the substring "DUF" (must be uppercase).
    - The annotation with the highest cumulative weight is chosen.
    - If no valid annotations exist or the top one is penalized, 
      the result defaults to "Unknown function".
    """
    votes: dict[str, float] = {}
    word_pattern = re.compile(r"\b(hypothetical|unknown)\b", flags=re.IGNORECASE)
    duf_pattern = re.compile(r"DUF")

    for item in annot_score_list:
        annotation = (item.get("annotation") or "").strip()
        score = item.get("score", 0.0)
        if not annotation or not isinstance(score, (int, float)):
            continue
        weight = score
        if word_pattern.search(annotation) or duf_pattern.search(annotation):
            weight *= 0.4
        votes[annotation] = votes.get(annotation, 0.0) + weight

    if not votes:
        return "Unknown function"

    winner, _ = max(votes.items(), key=lambda kv: kv[1])
    if word_pattern.search(winner) or duf_pattern.search(winner):
        return "Unknown function"
    return winner

def main():
    build_or_annotate = snakemake.params.build_or_annotate
    all_genes_path = snakemake.params.all_genes_annotated
    gene_index_path = snakemake.params.gene_index
    metabolism_path = snakemake.params.metabolism_table
    physiology_path = snakemake.params.physiology_table
    regulatory_path = snakemake.params.regulation_table
    final_table_path = snakemake.params.final_table
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
        
    all_genes_df = pl.read_csv(all_genes_path, separator='\t')
    gene_index_df = pl.read_csv(gene_index_path, separator='\t')
    gene_index_df = gene_index_df.select(["protein"] + [
        col for col in gene_index_df.columns if col.endswith("_protein_cluster") or col.endswith("_genome_cluster")
        ]).rename({"protein": "Protein"})
    metabolism_df = pl.read_csv(metabolism_path, separator='\t')
    physiology_df = pl.read_csv(physiology_path, separator='\t')
    regulatory_df = pl.read_csv(regulatory_path, separator='\t')
    
    if build_or_annotate == "build":
        logger.info(f"Generating the final table with proteins, families, auxiliary status, annotations, and classifications...")
        
        aux_scores_path = snakemake.params.aux_scores
        hmm_path = snakemake.params.hmm_table
        hmm_results = snakemake.params.hmm_results
        acc_prefix = snakemake.params.acc_prefix

        aux_scores_df = pl.read_csv(aux_scores_path, separator='\t')
        hmm_df = pl.read_csv(hmm_path, separator='\t')

        # 1) Build HMM-level metrics
        aux_metric_df = build_protein_metrics(aux_scores_df, hmm_df)
        
        # 2) Build HMM-level auxiliary status
        family_status_df = build_family_aux_status(aux_scores_df, hmm_df)

        # 3) Merge with everything else
        final_df = merge_dataframes_build(
            aux_metric_df, hmm_df, all_genes_df, gene_index_df,
            metabolism_df, physiology_df, regulatory_df, family_status_df
        )
        
        # Sort and write the final table output
        classification_order = {"metabolic": 0, "physiological": 1, "regulatory": 2, "unclassified": 3}
        final_df = final_df.with_columns(
            pl.col("Protein Classification").replace(classification_order).cast(pl.Int32).alias("Protein Classification_sort")
        )
        
        status_order = {
            "strongly auxiliary": 0,
            "weakly auxiliary": 1,
            "conditional": 2,
            "weakly core": 3,
            "strongly core": 4
        }
        final_df = final_df.with_columns(
            pl.col("Auxiliary Status").replace(status_order).cast(pl.Int32).alias("Auxiliary Status_sort")
        )
        
        confidence_order = {"high": 0, "medium": 1, "low": 2}
        final_df = final_df.with_columns(
            pl.col("Protein Viral Origin Confidence").replace(confidence_order).cast(pl.Int32).alias("Protein Viral Origin Confidence_sort")
        )
        
        final_df = final_df.sort(
            ["Protein Classification_sort",
             "Auxiliary Status_sort",
             "Size",
             "Family",
             "Protein Viral Origin Confidence_sort",
             "Protein"
             ], descending=[False, False, True, False, False, False])
        final_df = final_df.drop(
            ["Protein Classification_sort",
             "Auxiliary Status_sort",
             "Protein Viral Origin Confidence_sort"
             ])
        
        final_df.write_csv(final_table_path.replace(".tsv", "_detailed.tsv"), separator='\t')
        final_df_simple = final_df.drop(
            [col for col in final_df.columns if \
                col.endswith("Quantile") or \
                    col.endswith("Auxiliary Score") or \
                        col.endswith("Cluster")] + \
                            ["Size"]
        )
        final_df_simple.write_csv(final_table_path, separator='\t')

        # Count proteins per final Auxiliary Status
        # Define the order for auxiliary status (lowest means strongest auxiliary)
        aux_order = ["strongly auxiliary", "weakly auxiliary", "conditional", "weakly core", "strongly core"]
        # Define the order for protein classification categories
        class_order = ["metabolic", "physiological", "regulatory", "unclassified"]
        protein_designation_counts = final_df.group_by("Auxiliary Status").agg(
            pl.col("Protein").n_unique().alias("num_proteins")
        )
        # Overall Protein Counts by Auxiliary Status
        protein_counts_list = protein_designation_counts.to_dicts()
        protein_counts_list = sorted(
            protein_counts_list,
            key=lambda d: aux_order.index(d["Auxiliary Status"]) if d["Auxiliary Status"] in aux_order else len(aux_order)
        )
        overall_protein_report = [f"{d['Auxiliary Status'].title()}: {d['num_proteins']:,}" for d in protein_counts_list]
        logger.info(f"PROTEIN counts by auxiliary status:")
        for report in overall_protein_report:
            logger.info(report)

        # Family (HMM) Counts by Auxiliary Status
        family_designation_counts = family_status_df.group_by("auxiliary_status").agg(
            pl.col("hmm").n_unique().alias("num_families")
        )
        family_counts_list = family_designation_counts.to_dicts()
        family_counts_list = sorted(
            family_counts_list,
            key=lambda d: aux_order.index(d["auxiliary_status"]) if d["auxiliary_status"] in aux_order else len(aux_order)
        )
        overall_family_report = [f"{d['auxiliary_status'].title()}: {d['num_families']:,}" for d in family_counts_list]
        logger.info(f"FAMILY (HMM) counts by auxiliary status:")
        for report in overall_family_report:
            logger.info(report)

        # Protein Counts by Classification and Auxiliary Status (separate string per classification)
        protein_by_classification = final_df.group_by(["Protein Classification", "Auxiliary Status"]).agg(
            pl.col("Protein").n_unique().alias("num_proteins")
        )
        protein_classification_list = protein_by_classification.to_dicts()
        from collections import defaultdict
        class_grouped = defaultdict(list)
        for entry in protein_classification_list:
            class_grouped[entry["Protein Classification"]].append(entry)

        for class_cat in class_order:
            if class_cat in class_grouped:
                entries = sorted(
                    class_grouped[class_cat],
                    key=lambda d: aux_order.index(d["Auxiliary Status"]) if d["Auxiliary Status"] in aux_order else len(aux_order)
                )
                reports = [f"{d['Auxiliary Status'].title()}: {d['num_proteins']:,}" for d in entries]
                logger.info(f"{class_cat.upper()} protein counts by auxiliary status:")
                for report in reports:
                    logger.info(f"{report}")
        logger.info(f"Protein-level table written to {final_table_path}")
        
        # Write the final table to DuckDB
        duck_path = final_table_path.replace(".tsv", ".duckdb")
        con = duckdb.connect(duck_path)
        con.execute(f"DROP TABLE IF EXISTS {acc_prefix}")
        con.register(f"df", final_df)
        con.execute(f"CREATE TABLE {acc_prefix} AS SELECT * FROM df")
        con.close()
        logger.info(f"Protein-level DuckDB database written to {duck_path}")
        
        # Make a family-level table
        hmm_results_df = pl.read_csv(hmm_results, separator="\t").rename({"sequence": "Protein"})
        anno_cols = [
            ("KEGG KO Name", "KEGG"),
            ("Pfam Name", "Pfam"),
            ("PHROG Annotation", "PHROG"),
            ("METABOLIC Annotation", "METABOLIC"),
            ("CAZy Activities", "dbCAN"),
            ("Top Hit Annotation", None) # Top hit db is stored in another column
        ]

        # Collect all annotations + their associated db
        annot_long_list = []
        for col, db in anno_cols:
            if col not in final_df.columns:
                continue
            if db is not None:
                temp = final_df.select(["Protein", "Family", col]).rename({col: "annotation"}).with_columns(
                    pl.lit(db).alias("db")
                )
            else:
                # For "Top Hit Annotation", get db from "Top Hit HMM Origin"
                temp = final_df.select(["Protein", "Family", col, "Top Hit HMM Origin"]).rename({
                    col: "annotation",
                    "Top Hit HMM Origin": "db"
                })
            annot_long_list.append(temp)

        # Combine all annotation columns into a single long-format table
        annot_long_df = pl.concat(annot_long_list).filter(pl.col("annotation").is_not_null())

        # Join with bitscores
        bitscore_annot = annot_long_df.join(hmm_results_df, on=["Protein", "db"], how="inner")

        # Group by Family, collect annotations + score
        family_annos_weighted = (
            bitscore_annot
            .group_by("Family")
            .agg(
                pl.struct(["annotation", "score"]).alias("annotation_score_list")
            )
        )

        # Extract just the column names
        anno_colnames = [col for col, _ in anno_cols]

        family_level_df = final_df.with_columns(
            pl.struct(anno_colnames).alias("anno_struct")
        )

        # Group into Family-level summary
        family_level_df = (
            family_level_df
            .group_by("Family", maintain_order=True)
            .agg([
                pl.len().alias("Size"),
                pl.col("Protein Classification").alias("class_list"),
                pl.col("anno_struct").alias("anno_structs"),
            ])
        )

        family_level_df = family_level_df.join(family_annos_weighted, on="Family", how="left")
        
        # Convert to Python for applying consensus logic
        df_dicts = []
        for row in family_level_df.iter_rows(named=True):
            classification = consensus_classification(row["class_list"])
            if row["annotation_score_list"] is None:
                predicted_function = "Unknown function"
            else:
                predicted_function = consensus_annotation(row["annotation_score_list"])
            df_dicts.append({
                "Family": row["Family"],
                "Size": row["Size"],
                "Classification": classification,
                "Predicted Function": predicted_function,
            })

        # Convert back to Polars
        family_level_df = pl.DataFrame(df_dicts)
        family_level_df = family_level_df.join(family_status_df.rename({"hmm": "Family", "auxiliary_status": "Auxiliary Status"}).drop("avg_quantile"), on="Family", how="left")
        # Add auxiliary status and sort by custom logic
        classification_order = {"metabolic": 0, "physiological": 1, "regulatory": 2, "unclassified": 3}
        aux_order = {
            "strongly auxiliary": 0,
            "weakly auxiliary": 1,
            "conditional": 2,
            "weakly core": 3,
            "strongly core": 4
        }

        family_level_df = family_level_df.with_columns([
            pl.col("Classification").replace(classification_order).cast(pl.Int32).alias("Classification_sort"),
            pl.col("Auxiliary Status").replace(aux_order).cast(pl.Int32).alias("Auxiliary_Status_sort")
        ])

        family_level_df = family_level_df.sort(
            by=["Classification_sort", "Size", "Auxiliary_Status_sort"],
            descending=[False, True, False]
        ).drop(["Classification_sort", "Auxiliary_Status_sort"]).select(
            ["Family", "Size", "Auxiliary Status", "Classification", "Predicted Function"]
        )

        family_level_path = final_table_path.replace(".tsv", "_protein_families.tsv")
        family_level_df.write_csv(family_level_path, separator="\t")

        logger.info(f"Family-level table written to {family_level_path}")
        
        # Write the family-level table to DuckDB
        duck_fam_path = family_level_path.replace(".tsv", ".duckdb")
        con = duckdb.connect(duck_fam_path)
        con.execute(f"DROP TABLE IF EXISTS {acc_prefix}_protein_families")
        con.register(f"df", family_level_df)
        con.execute(f"CREATE TABLE {acc_prefix}_protein_families AS SELECT * FROM df")
        con.close()
        logger.info(f"Family-level DuckDB database written to {duck_fam_path}")
            
    elif build_or_annotate == "annotate":
        logger.info(f"Generating the final table with proteins, annotations, and classifications...")
        
        final_df = merge_dataframes_genome(all_genes_df, metabolism_df, physiology_df, regulatory_df)
        
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
    else:
        logging.error(f"Invalid option for build_or_annotate: {build_or_annotate}")
        raise ValueError(f"Invalid option for build_or_annotate: {build_or_annotate}")
    
    # Log classification summary
    logger.debug(f"Columns in final_df after classification: {final_df.columns}")
    logger.debug(f"Protein Classification value counts:\n{final_df['Protein Classification'].value_counts()}")
    logger.debug(f"Protein Viral Origin Confidence value counts:\n{final_df['Protein Viral Origin Confidence'].value_counts()}")
    if build_or_annotate == "build":
        logger.debug(f"Auxiliary Status value counts:\n{final_df['Auxiliary Status'].value_counts()}")
        # Log summary of families with multiple unique protein classifications
        family_classification = final_df.group_by("Family").agg(
            pl.col("Protein Classification").n_unique().alias("unique_classifications")
        )
        multiple_classifications = family_classification.filter(pl.col("unique_classifications") > 1)
        if multiple_classifications.height > 0:
            logger.warning(f"There are {multiple_classifications.height} families with >1 unique Protein Classifications.")
            classifications = final_df.join(multiple_classifications, on="Family", how="inner")\
                .filter(pl.col("unique_classifications") > 1).sort(["Family", "Protein Classification"])
            multiple_output = final_table_path.replace('.tsv', '_conflicting_classifications.tsv')
            logger.warning(f"Writing families with multiple unique Protein Classifications to {multiple_output}")
            classifications.write_csv(multiple_output, separator='\t')
            logger.debug(f"Families with multiple unique Protein Classifications:\n{multiple_classifications}")
            pl.Config.set_tbl_cols(20)
            logger.debug(f"Classifications: {classifications.drop([col for col in classifications.columns if \
                col.endswith("Quantile") or \
                    col.endswith("Auxiliary Score") or \
                        col.endswith("Cluster")] + \
                            ["Size", "Contig", "Genome", "unique_classifications"])}")
        else:
            logger.debug("No families with multiple unique Protein Classifications found.")
    
if __name__ == "__main__":
    main()