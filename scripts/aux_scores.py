#!/usr/bin/env python3
import os
import sys
import resource
import platform
import logging
import numba
import numpy as np
from numpy.typing import NDArray
from scipy.sparse import csr_matrix
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
logging.getLogger("numba").setLevel(logging.INFO)

print("========================================================================\n       Step 16/22: Calculate auxiliary scores of protein clusters       \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n       Step 16/22: Calculate auxiliary scores of protein clusters       \n========================================================================\n")

# Constants
EPSILON = 1e-16 # Small constant to avoid log(0)

# Types
FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]

def entropy(distr: FloatArray, axis: int | None = None) -> float | FloatArray:
    """Calculate entropy for dense arrays."""
    summand = distr * np.log2(distr + EPSILON)
    if axis is None:
        result = max(-np.sum(summand), 0.0)
    else:
        result = -np.sum(summand, axis=axis)
        result = np.maximum(result, 0.0)
    return result  # type: ignore

def sparse_entropy(sparse_matrix: csr_matrix, axis: int = 1) -> np.ndarray:
    """Calculate entropy along the given axis of a sparse matrix."""
    log_sparse_matrix = sparse_matrix.copy()
    log_sparse_matrix.data = np.log2(np.asarray(log_sparse_matrix.data) + EPSILON)
    summand = sparse_matrix.multiply(log_sparse_matrix)
    entropy_values = -summand.sum(axis=axis).A1
    return np.maximum(entropy_values, 0.0)

def uniform_distr(n: int) -> FloatArray:
    return np.ones(n) / n

def calc_background_entropy(n_genomes: int) -> float:
    """Calculate the entropy of a uniform distribution as the background entropy."""
    background_distr = uniform_distr(n_genomes)
    return entropy(background_distr)

def create_membership_array(n_genomes: int, genome_sizes: IntArray, protein_cluster_label: IntArray) -> csr_matrix:
    if len(genome_sizes) != n_genomes:
        raise ValueError("Length of genome_sizes must be equal to n_genomes")
    if len(protein_cluster_label) != genome_sizes.sum():
        raise ValueError("Length of protein_cluster_label must be equal to the sum of genome_sizes")
    unique_labels, pseudo_labels = np.unique(protein_cluster_label, return_inverse=True)
    n_clusters = len(unique_labels)
    ptn_label = np.arange(n_genomes).repeat(genome_sizes)
    membership = csr_matrix((np.ones(len(pseudo_labels)), (pseudo_labels, ptn_label)), shape=(n_clusters, n_genomes))
    return membership, unique_labels

def process_genome_cluster(genome_cluster: str, rank: str, rank_df: pl.DataFrame) -> pl.DataFrame | None:
    local_df = rank_df.filter(pl.col(f"{rank}_genome_cluster") == genome_cluster)
    local_df = local_df.sort("genome")
    unique_genomes = local_df.select("genome").unique().to_numpy()
    logger.debug(f"[{rank}] Genome cluster {genome_cluster} - unique genomes (sorted): {unique_genomes}")

    genome_size_df = local_df.group_by("genome").agg(pl.count("protein").alias("genome_size")).sort("genome")
    genome_sizes = genome_size_df["genome_size"].to_numpy()
    logger.debug(f"[{rank}] Genome sizes for cluster {genome_cluster}: {genome_sizes}")

    local_df = local_df.filter(pl.col(f"{rank}_protein_cluster").is_not_null())
    if local_df.is_empty():
        logger.warning(f"[{rank}] Genome cluster {genome_cluster} has no valid protein cluster labels; skipping.")
        return None

    protein_cluster_labels = local_df.select(f"{rank}_protein_cluster").to_numpy()
    logger.debug(f"[{rank}] Protein cluster labels for cluster {genome_cluster}: {protein_cluster_labels}")

    n_genomes = len(genome_sizes)
    expected_proteins = genome_sizes.sum()
    actual_proteins = len(protein_cluster_labels)
    if expected_proteins != actual_proteins:
        logger.warning(
            f"[{rank}] Skipping genome cluster {genome_cluster} due to mismatch: expected proteins ({expected_proteins}) "
            f"vs protein_cluster_labels count ({actual_proteins})."
        )
        return None

    membership, unique_labels = create_membership_array(n_genomes, genome_sizes.flatten(), protein_cluster_labels.flatten())
    logger.debug(f"[{rank}] Membership matrix shape for cluster {genome_cluster}: {membership.shape}")
    membership_row_sums = membership.sum(axis=1).A1
    logger.debug(f"[{rank}] Membership matrix row sums for cluster {genome_cluster}: {membership_row_sums}")

    membership_sum = membership.sum(axis=1).A1
    pc_distrs = membership.multiply(1 / membership_sum[:, None])
    pc_entropy = sparse_entropy(pc_distrs, axis=1)

    weights = membership.multiply(1 / genome_sizes.flatten())
    weights = weights.sum(axis=1).A1
    weights /= weights.max()

    background_entropy = calc_background_entropy(n_genomes) + EPSILON
    gain = background_entropy - pc_entropy
    gain_ratio = gain / background_entropy
    aux_scores = 1.0 - gain_ratio

    weighted_gain = background_entropy - (weights * pc_entropy)
    weighted_gain_ratio = weighted_gain / background_entropy
    aux_scores_weighted = 1.0 - weighted_gain_ratio

    pseudo_label_map = {i: label for i, label in enumerate(unique_labels)}
    result_df = pl.DataFrame({
        f"{rank}_protein_cluster": [pseudo_label_map[i] for i in range(len(aux_scores))],
        f"protein_cluster_aux_score_{rank}": aux_scores,
        f"protein_cluster_aux_score_{rank}_weighted": aux_scores_weighted,
    })
    logger.debug(f"[{rank}] Processed genome cluster {genome_cluster}: result_df shape: {result_df.shape}")
    return result_df

def calculate_entropy_based_aux_scores(df: pl.DataFrame, exclude_singletons: bool, ranks: list[str], threads: int) -> pl.DataFrame:
    # Use lazy evaluation to delay and fuse operations.
    lazy_df = df.lazy().unique()
    df = lazy_df.collect()
    df = df.select([
        'genome',
        'protein',
        *[f"{rank}_genome_cluster" for rank in ranks],
        *[f"{rank}_protein_cluster" for rank in ranks],
    ])
    all_data = []
    logger.debug(f"All genes dataframe after unique(): {df.shape}")
    logger.debug(f"Columns in the dataframe: {df.columns}")

    for rank in sorted(ranks):
        logger.debug(f"Processing dataframe for rank: {rank}")
        
        # Lazily select only the relevant columns
        lazy_rank = df.lazy().select([
            'genome',
            'protein',
            f"{rank}_genome_cluster",
            f"{rank}_protein_cluster",
        ])

        if exclude_singletons:
            logger.info(f"Excluding singleton genome clusters at {rank}-level...")
            genome_cluster_size_df = (
                lazy_rank
                .unique(subset=[f"{rank}_genome_cluster", "genome"])
                .group_by(f"{rank}_genome_cluster")
                .agg(pl.col("genome").n_unique().alias("genome_cluster_size"))
            )

            # Filter out singleton genome clusters
            lazy_rank = lazy_rank.join(genome_cluster_size_df, on=f"{rank}_genome_cluster", how="left")
            lazy_rank = lazy_rank.filter(pl.col("genome_cluster_size") > 1)

        # Collect the result of all fused operations into memory
        rank_df = lazy_rank.collect()

        if exclude_singletons:
            # Continue with already eager `rank_df` for protein cluster filtering
            singleton_gc = genome_cluster_size_df.filter(pl.col("genome_cluster_size") == 1).collect()
            sc_list = singleton_gc[f"{rank}_genome_cluster"].to_numpy()
            n_singleton_genomes = df.filter(pl.col(f"{rank}_genome_cluster").is_in(sc_list)).select("genome").n_unique()
            n_singleton_genome_prots = df.filter(pl.col(f"{rank}_genome_cluster").is_in(sc_list)).select("protein").n_unique()
            logger.info(f"Filtered {n_singleton_genome_prots:,} proteins from {n_singleton_genomes:,} singleton genomes at the {rank}-level.")

            logger.info(f"Excluding singleton protein clusters at {rank}-level...")
            protein_cluster_size_df = (
                rank_df.group_by(f"{rank}_protein_cluster")
                .agg(pl.col("protein").n_unique().alias("protein_cluster_size"))
            )
            sp_ids = protein_cluster_size_df.filter(pl.col("protein_cluster_size") == 1)[f"{rank}_protein_cluster"].to_numpy()
            rank_df = rank_df.filter(~pl.col(f"{rank}_protein_cluster").is_in(sp_ids))
            n_singleton_protein_prots = df.filter(pl.col(f"{rank}_protein_cluster").is_in(sp_ids)).select("protein").n_unique()
            logger.info(f"Filtered {n_singleton_protein_prots:,} singleton proteins at the {rank}-level.")
            logger.info(f"Number of remaining proteins at the {rank}-level: {rank_df.select('protein').n_unique():,}")

        # Apply uniqueness and not-null filters
        rank_df = rank_df.unique().filter(pl.col(f"{rank}_protein_cluster").is_not_null())
        genome_clusters = rank_df.select(f"{rank}_genome_cluster").unique().to_numpy()
        genome_cluster_aux_scores_dfs = []
        with ThreadPoolExecutor(max_workers=threads) as executor:
            future_to_cluster = {executor.submit(process_genome_cluster, gc, rank, rank_df): gc for gc in genome_clusters}
            for future in as_completed(future_to_cluster):
                result_df = future.result()
                if result_df is not None:
                    genome_cluster_aux_scores_dfs.append(result_df)
        aux_scores_df = pl.concat(genome_cluster_aux_scores_dfs)
        rank_df = rank_df.join(aux_scores_df, on=f"{rank}_protein_cluster", how='left')
        rank_df = rank_df.filter(pl.col("protein").is_not_null())
        rank_df = rank_df.rename({
            f"{rank}_genome_cluster": "genome_cluster",
            f"{rank}_protein_cluster": "protein_cluster",
            f"protein_cluster_aux_score_{rank}": "aux_score",
            f"protein_cluster_aux_score_{rank}_weighted": "aux_score_weighted"
        })
        rank_df = rank_df.with_columns(pl.lit(rank).alias('rank')).unique()
        all_data.append(rank_df)
    combined_df = pl.concat(all_data)
    combined_df = combined_df.select([
        "protein", "rank", "genome_cluster", "protein_cluster", "aux_score", "aux_score_weighted"
    ]).sort("protein", "rank").unique()
    return combined_df

def assign_quantile_bins(scores: np.ndarray, quantiles: np.ndarray) -> np.ndarray:
    bins = np.digitize(scores, quantiles, right=False)
    bins[bins > 10] = 10
    return bins

@numba.njit
def make_unique_quantiles(quantiles: np.ndarray) -> np.ndarray:
    jitter_amount = EPSILON
    unique_quantiles = np.copy(quantiles)
    for i in range(1, len(unique_quantiles)):
        if unique_quantiles[i] <= unique_quantiles[i - 1]:
            unique_quantiles[i] = unique_quantiles[i - 1] + jitter_amount
    return unique_quantiles

def calculate_quantiles_per_rank(df: pl.DataFrame, ranks: list[str]) -> pl.DataFrame:
    all_data = []
    for rank in ranks:
        rank_df = df.filter(pl.col("rank") == rank)
        if rank_df.is_empty():
            all_data.append(rank_df)
            continue
        aux_scores = rank_df["aux_score"].to_numpy()
        aux_scores_weighted = rank_df["aux_score_weighted"].to_numpy()
        quantiles = np.quantile(aux_scores, np.linspace(0, 1, 11))
        quantiles_weighted = np.quantile(aux_scores_weighted, np.linspace(0, 1, 11))
        quantiles = make_unique_quantiles(quantiles)
        quantiles_weighted = make_unique_quantiles(quantiles_weighted)
        bins_all = assign_quantile_bins(aux_scores, quantiles)
        bins_weighted_all = assign_quantile_bins(aux_scores_weighted, quantiles_weighted)
        rank_df = rank_df.with_columns([
            pl.Series("quantile", bins_all),
            pl.Series("quantile_weighted", bins_weighted_all)
        ])
        all_data.append(rank_df)
    combined_df = pl.concat(all_data, how="vertical") if all_data else pl.DataFrame()
    return combined_df

def main():
    outdir_path = snakemake.params.outdir
    all_genes_path = snakemake.params.gene_index
    cluster_ranks = snakemake.params.cluster_taxa_levels
    aux_scores_path = snakemake.params.aux_scores
    exclude_singletons = snakemake.params.exclude_singletons
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    logger.info("Calculating auxiliary scores of input proteins...")
    logger.debug(f"Maximum memory allowed to be allocated: {mem_limit} GB")
    set_memory_limit(mem_limit)
    
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)
    
    # Read input in lazy mode and then collect unique rows
    lazy_all_genes = pl.scan_csv(all_genes_path, separator='\t')
    all_genes = lazy_all_genes.unique().collect()
    
    aux_scores_df = calculate_entropy_based_aux_scores(all_genes, exclude_singletons, cluster_ranks, threads)
    aux_scores_df_with_quantiles = calculate_quantiles_per_rank(aux_scores_df, cluster_ranks)
    aux_scores_df_with_quantiles.write_csv(aux_scores_path, separator='\t')
    logger.info("Calculation of auxiliary scores completed.")

if __name__ == "__main__":
    main()
