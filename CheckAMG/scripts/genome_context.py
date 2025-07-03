import os
import sys
import resource
import platform
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
import numpy as np
import pandas as pd
from joblib import load
from numba import njit
import concurrent.futures
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import logging

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
logging.getLogger("numba").setLevel(logging.INFO)

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

print("========================================================================\n         Step 8/11: Analyze the genomic context of annotations          \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n         Step 8/11: Analyze the genomic context of annotations          \n========================================================================\n")

def calculate_gene_lengths(data):
    """
    Calculate gene lengths and protein lengths in amino acids.
    """
    data = data.with_columns([
        (pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1).alias('gene_length_bases'),
        ((pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1) / 3).cast(pl.Int32).alias('prot_length_AAs')
    ])
    return data

def calculate_contig_statistics(data, circular_contigs):
    """
    Calculate contig average V-scores/VL-scores and assign a circular_contig flag.
    """
    stats = data.group_by("contig", maintain_order=True).agg([
        pl.col("KEGG_V-score").mean().alias("contig_avg_KEGG_V-score"),
        pl.col("Pfam_V-score").mean().alias("contig_avg_Pfam_V-score"),
        pl.col("PHROG_V-score").mean().alias("contig_avg_PHROG_V-score"),
        pl.col("KEGG_VL-score").mean().alias("contig_avg_KEGG_VL-score"),
        pl.col("Pfam_VL-score").mean().alias("contig_avg_Pfam_VL-score"),
        pl.col("PHROG_VL-score").mean().alias("contig_avg_PHROG_VL-score")
    ])
    result = data.join(stats, on="contig")
    result = result.with_columns(pl.col("contig").is_in(circular_contigs).alias("circular_contig"))
    return result

@njit
def window_avg(scores, lengths, window_size, minimum_percentage):
    """
    Two-pointer method to calculate average V/VL-scores within a variable-length window.
    """
    n = len(lengths)
    out = np.full(n, np.nan, dtype=np.float64)
    prefix_len = np.zeros(n+1, dtype=np.float64)
    prefix_score = np.zeros(n+1, dtype=np.float64)
    prefix_valid_len = np.zeros(n+1, dtype=np.float64)
    prefix_count = np.zeros(n+1, dtype=np.float64)

    for i in range(n):
        prefix_len[i+1] = prefix_len[i] + lengths[i]
        if not np.isnan(scores[i]):
            prefix_score[i+1] = prefix_score[i] + scores[i]
            prefix_valid_len[i+1] = prefix_valid_len[i] + lengths[i]
            prefix_count[i+1] = prefix_count[i] + 1
        else:
            prefix_score[i+1] = prefix_score[i]
            prefix_valid_len[i+1] = prefix_valid_len[i]
            prefix_count[i+1] = prefix_count[i]

    left_ptr = 0
    right_ptr = 0
    for i in range(n):
        while prefix_len[i] - prefix_len[left_ptr] > window_size:
            left_ptr += 1
        while right_ptr + 1 < n and prefix_len[right_ptr+1] - prefix_len[i+1] < window_size:
            right_ptr += 1
        total_len = prefix_len[right_ptr+1] - prefix_len[left_ptr]
        if total_len == 0:
            out[i] = np.nan
            continue
        valid_len = prefix_valid_len[right_ptr+1] - prefix_valid_len[left_ptr]
        pct_valid = 100.0 * valid_len / total_len
        if pct_valid >= minimum_percentage:
            sum_scores = prefix_score[right_ptr+1] - prefix_score[left_ptr]
            count_valid = prefix_count[right_ptr+1] - prefix_count[left_ptr]
            if count_valid > 0:
                out[i] = sum_scores / count_valid
            else:
                out[i] = np.nan
        else:
            out[i] = np.nan
    return out

def process_window_statistics_for_contig(contig, data, window_size, minimum_percentage):
    """
    Calculate window averages for KEGG_VL-score, Pfam_VL-score, PHROG_VL-score,
    KEGG_V-score, Pfam_V-score, PHROG_V-score using a two-pointer approach for
    a single contig.
    """
    df = data.filter(pl.col("contig") == contig)
    lengths = df["gene_length_bases"].to_numpy()
    kegg_vl = df["KEGG_VL-score"].to_numpy()
    pfam_vl = df["Pfam_VL-score"].to_numpy()
    phrog_vl = df["PHROG_VL-score"].to_numpy()
    kegg_v = df["KEGG_V-score"].to_numpy()
    pfam_v = df["Pfam_V-score"].to_numpy()
    phrog_v = df["PHROG_V-score"].to_numpy()

    if len(lengths) == 0:
        return df

    kegg_vl_out = window_avg(kegg_vl, lengths, window_size, minimum_percentage)
    pfam_vl_out = window_avg(pfam_vl, lengths, window_size, minimum_percentage)
    phrog_vl_out = window_avg(phrog_vl, lengths, window_size, minimum_percentage)
    kegg_v_out = window_avg(kegg_v, lengths, window_size, minimum_percentage)
    pfam_v_out = window_avg(pfam_v, lengths, window_size, minimum_percentage)
    phrog_v_out = window_avg(phrog_v, lengths, window_size, minimum_percentage)

    df = df.with_columns([
        pl.Series("window_avg_KEGG_VL-score", kegg_vl_out),
        pl.Series("window_avg_Pfam_VL-score", pfam_vl_out),
        pl.Series("window_avg_PHROG_VL-score", phrog_vl_out),
        pl.Series("window_avg_KEGG_V-score", kegg_v_out),
        pl.Series("window_avg_Pfam_V-score", pfam_v_out),
        pl.Series("window_avg_PHROG_V-score", phrog_v_out)
    ])
    return df

def calculate_window_statistics(data, window_size, minimum_percentage, n_cpus):
    """
    Calculate window averages for the entire dataset by processing each contig in parallel.
    """
    data = data.sort(["contig", "gene_number"])
    contigs = data["contig"].unique().to_list()

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(
                process_window_statistics_for_contig,
                contig, data, window_size, minimum_percentage
            )
            for contig in contigs
        ]
        results = [f.result() for f in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Calculating sliding-window averages", unit="contig")]
    return pl.concat(results, how="vertical")

def prepare_lgbm_features(df, feature_names):
    # Ensure required features are present, fill missing
    features = {}
    for col in feature_names:
        if col in df.columns:
            features[col] = df[col].to_numpy()
        else:
            # fill missing with 0.0 for floats, 0 for ints/bools
            features[col] = np.zeros(len(df), dtype=float)
    X = pd.DataFrame({c: features[c] for c in feature_names})

    # Identify columns
    numeric_cols = [c for c in X.columns if X[c].dtype != bool and not ('flank' in c.lower() or 'circular' in c)]
    boolean_cols = [c for c in X.columns if X[c].dtype == bool or ('flank' in c.lower() or 'circular' in c)]

    # Group by contig and impute per group
    for col in numeric_cols:
        X[col] = X.groupby(df['contig'])[col].transform(lambda x: x.fillna(x.median()))
    for col in boolean_cols:
        # Convert to int first if necessary
        if X[col].dtype != int:
            X[col] = X[col].astype(float)
        X[col] = X.groupby(df['contig'])[col].transform(lambda x: x.fillna(x.mode().iloc[0] if not x.mode().empty else 0))

    # Ensure order of columns
    X = X[feature_names]
    return X

def viral_origin_confidence_lgbm(df, lgbm_model, thresholds, feature_names):
    """
    df: pl.DataFrame (must have columns matching those used in model training)
    lgbm_model: fitted sklearn Light GM or CalibratedClassifierCV
    thresholds: dict, { 'high': {'threshold': float, ...}, 'medium': {'threshold': float, ...}, 'low': {'threshold': float, ...} }
    feature_names: list of feature columns (ordered)
    imputer_num: fitted SimpleImputer for numeric cols
    imputer_bool: fitted SimpleImputer for bool cols
    Returns  polars DataFrame with added column 'Viral_Origin_Confidence' (high/medium/low)
    """
    # Convert polars to pandas for sklearn
    df_pd = df.to_pandas()

    # Prepare features
    X = prepare_lgbm_features(df_pd, feature_names)

    # Predict proba
    y_proba = lgbm_model.predict_proba(X)[:, 1]
    # Assign confidence
    conf = np.full(y_proba.shape, 'low', dtype=object)
    conf[y_proba >= thresholds['medium']['threshold']] = 'medium'
    conf[y_proba >= thresholds['high']['threshold']] = 'high'
    # Add to polars DataFrame
    df = df.with_columns([
        pl.Series('LGBM_viral_prob', y_proba),
        pl.Series('Viral_Origin_Confidence', conf)
    ])
    return df

@njit
def flank_two_pointer_vscores(lengths, scores, max_flank_length, min_vscore):
    """
    Two-pointer approach to compute separate upstream and downstream flags for v-scores.
    Returns two boolean arrays (left/out/up and right/down) for each gene.
    """
    n = len(lengths)
    left_out = np.zeros(n, dtype=np.bool_)
    right_out = np.zeros(n, dtype=np.bool_)
    prefix_len = np.zeros(n+1, dtype=np.float64)
    prefix_meet = np.zeros(n+1, dtype=np.int64)
    for i in range(n):
        prefix_len[i+1] = prefix_len[i] + lengths[i]
        if not np.isnan(scores[i]) and scores[i] >= min_vscore:
            prefix_meet[i+1] = prefix_meet[i] + 1
        else:
            prefix_meet[i+1] = prefix_meet[i]
    left_ptr = 0
    right_ptr = 0
    for i in range(n):
        # Check upstream flank (genes before i)
        while prefix_len[i] - prefix_len[left_ptr] > max_flank_length:
            left_ptr += 1
        left_has = (prefix_meet[i] - prefix_meet[left_ptr]) > 0
        left_out[i] = left_has
        # Check downstream flank (genes after i)
        while right_ptr + 1 < n and prefix_len[right_ptr+1] - prefix_len[i+1] <= max_flank_length:
            right_ptr += 1
        right_has = (prefix_meet[right_ptr+1] - prefix_meet[i+1]) > 0
        right_out[i] = right_has
    return left_out, right_out

def verify_flanking_vscores(contig_data, minimum_vscore, max_flank_length):
    """
    Verify flanking v-score values by checking upstream and downstream separately.
    Returns the contig_data joined with new columns:
      KEGG_verified_flank_up, KEGG_verified_flank_down,
      Pfam_verified_flank_up, Pfam_verified_flank_down,
      PHROG_verified_flank_up, PHROG_verified_flank_down.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_scores = contig_data['KEGG_V-score'].to_numpy()
    pfam_scores = contig_data['Pfam_V-score'].to_numpy()
    phrog_scores = contig_data['PHROG_V-score'].to_numpy()

    kegg_left, kegg_right = flank_two_pointer_vscores(lengths, kegg_scores, max_flank_length, minimum_vscore)
    pfam_left, pfam_right = flank_two_pointer_vscores(lengths, pfam_scores, max_flank_length, minimum_vscore)
    phrog_left, phrog_right = flank_two_pointer_vscores(lengths, phrog_scores, max_flank_length, minimum_vscore)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str]*len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_verified_flank_up": kegg_left,
        "KEGG_verified_flank_down": kegg_right,
        "Pfam_verified_flank_up": pfam_left,
        "Pfam_verified_flank_down": pfam_right,
        "PHROG_verified_flank_up": phrog_left,
        "PHROG_verified_flank_down": phrog_right
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

def create_in_set_array(hmm_ids, valid_set):
    """
    Returns a NumPy int array with 1 if hmm_ids[i] is in valid_set, else 0.
    This is done in Python space to avoid string membership checks in nopython.
    """
    arr = np.zeros(len(hmm_ids), dtype=np.int64)
    for i, val in enumerate(hmm_ids):
        # If val is None or not in set, it remains 0, otherwise 1
        if val is not None and val in valid_set:
            arr[i] = 1
    return arr

@njit
def flank_two_pointer_integers(lengths, in_set, max_flank_length):
    """
    Two-pointer approach over an integer array in_set (0 or 1) to compute separate upstream and downstream flags.
    Returns two boolean arrays for the left and right flanks.
    """
    n = len(lengths)
    left_out = np.zeros(n, dtype=np.bool_)
    right_out = np.zeros(n, dtype=np.bool_)
    prefix_len = np.zeros(n+1, dtype=np.float64)
    prefix_inset = np.zeros(n+1, dtype=np.int64)
    for i in range(n):
        prefix_len[i+1] = prefix_len[i] + lengths[i]
        prefix_inset[i+1] = prefix_inset[i] + in_set[i]
    left_ptr = 0
    right_ptr = 0
    for i in range(n):
        while prefix_len[i] - prefix_len[left_ptr] > max_flank_length:
            left_ptr += 1
        left_has = (prefix_inset[i] - prefix_inset[left_ptr]) > 0
        left_out[i] = left_has
        while right_ptr + 1 < n and prefix_len[right_ptr+1] - prefix_len[i+1] <= max_flank_length:
            right_ptr += 1
        right_has = (prefix_inset[right_ptr+1] - prefix_inset[i+1]) > 0
        right_out[i] = right_has
    return left_out, right_out

def verify_flanking_hallmark(contig_data, hallmark_accessions, max_flank_length):
    """
    Verify flanking hallmark genes by checking both upstream and downstream separately.
    Returns the contig_data joined with new columns:
      KEGG_verified_flank_up, KEGG_verified_flank_down,
      Pfam_verified_flank_up, Pfam_verified_flank_down,
      PHROG_verified_flank_up, PHROG_verified_flank_down.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()

    kegg_arr = create_in_set_array(kegg_hmm, hallmark_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, hallmark_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, hallmark_accessions)

    kegg_left, kegg_right = flank_two_pointer_integers(lengths, kegg_arr, max_flank_length)
    pfam_left, pfam_right = flank_two_pointer_integers(lengths, pfam_arr, max_flank_length)
    phrog_left, phrog_right = flank_two_pointer_integers(lengths, phrog_arr, max_flank_length)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str]*len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_verified_flank_up": kegg_left,
        "KEGG_verified_flank_down": kegg_right,
        "Pfam_verified_flank_up": pfam_left,
        "Pfam_verified_flank_down": pfam_right,
        "PHROG_verified_flank_up": phrog_left,
        "PHROG_verified_flank_down": phrog_right
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

def check_flanking_insertions(contig_data, mobile_accessions, max_flank_length):
    """
    Check for genes that have known mobile-element proteins in their left or right flank.
    Only one flank is needed to be positive (require_both_flanks=False).
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()

    # Create integer indicator arrays for mobile gene membership.
    kegg_arr = create_in_set_array(kegg_hmm, mobile_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, mobile_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, mobile_accessions)

    # Get left and right flank indicators (each is an array of booleans)
    kegg_left, kegg_right = flank_two_pointer_integers(lengths, kegg_arr, max_flank_length)
    pfam_left, pfam_right = flank_two_pointer_integers(lengths, pfam_arr, max_flank_length)
    phrog_left, phrog_right = flank_two_pointer_integers(lengths, phrog_arr, max_flank_length)

    # Combine left and right flags with a logical OR for each gene.
    kegg_combined = np.logical_or(kegg_left, kegg_right)
    pfam_combined = np.logical_or(pfam_left, pfam_right)
    phrog_combined = np.logical_or(phrog_left, phrog_right)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str] * len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_MGE_flank": kegg_combined,
        "Pfam_MGE_flank": pfam_combined,
        "PHROG_MGE_flank": phrog_combined
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

def process_genomes(data, circular_contigs, minimum_percentage,
                    window_size, max_flank_length, minimum_vscore,
                    lgbm_model, thresholds, feature_names,
                    use_hallmark=False,
                    hallmark_accessions=None, mobile_accessions=None,
                    n_cpus=1):
    logger.debug(f"Calculating lengths for {data.shape[0]:,} genes.")
    logger.debug(f"Data before calculating gene lengths: {data.head()}")
    data = calculate_gene_lengths(data)

    logger.info("Calculating contig statistics.")
    logger.debug(f"Data before calculating contig statistics: {data.head()}")
    data = calculate_contig_statistics(data, circular_contigs)

    logger.info("Calculating window statistics.")
    logger.debug(f"Data before calculating window statistics: {data.head()}")
    logger.debug(f"Column dtypes before conversion: {data.schema}")
    score_columns = [
        "KEGG_V-score","KEGG_VL-score","Pfam_V-score","Pfam_VL-score","PHROG_V-score","PHROG_VL-score",
        "contig_avg_KEGG_V-score","contig_avg_Pfam_V-score","contig_avg_PHROG_V-score",
        "contig_avg_KEGG_VL-score","contig_avg_Pfam_VL-score", "contig_avg_PHROG_VL-score"
    ]
    for col in score_columns:
        if col in data.columns:
            data = data.with_columns(pl.col(col).cast(pl.Float64, strict=False))
    logger.debug(f"Column dtypes after conversion: {data.schema}")

    # Parallel window statistics calculated per contig.
    data = calculate_window_statistics(data, window_size, minimum_percentage, n_cpus)
    data = data.unique()

    # Parallel verification of flanking regions by partitioning by contig.
    if use_hallmark and hallmark_accessions is not None:
        logger.info("Verifying flanking hallmark genes.")
        logger.debug(f"Data before verifying flanking hallmark genes: {data.head()}")
        contig_dfs = data.partition_by("contig")
        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            futures = [
                executor.submit(verify_flanking_hallmark, df, hallmark_accessions, max_flank_length)
                for df in contig_dfs
            ]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanks for viral hallmarks", unit="contig")]
        data = pl.concat(results, how="vertical")
    else:
        logger.info("Verifying flanking V-scores.")
        logger.debug(f"Data before verifying flanking V-scores: {data.head()}")
        contig_dfs = data.partition_by("contig")
        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            futures = [
                executor.submit(verify_flanking_vscores, df, minimum_vscore, max_flank_length)
                for df in contig_dfs
            ]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc=f"Checking flanks for V-score={minimum_vscore}", unit="contig")]
        data = pl.concat(results, how="vertical")

    logger.info("Checking for genes in potential mobile genetic element regions.")
    logger.debug(f"Data before checking for flanking insertions: {data.head()}")
    contig_dfs = data.partition_by("contig")
    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(check_flanking_insertions, df, mobile_accessions, max_flank_length)
            for df in contig_dfs
        ]
        results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanks for mobile genes", unit="contig")]
    data = pl.concat(results, how="vertical")
    
    logger.info("Assigning viral origin confidence using LightGBM with genome context features.")
    logger.debug(f"Data before assigning viral origin confidence: {data.head()}")
    contig_dfs = data.partition_by("contig")
    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(viral_origin_confidence_lgbm, df, lgbm_model, thresholds, feature_names)
            for df in contig_dfs
        ]
        results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Fitting models", unit="contig")]
    data = pl.concat(results, how="vertical")

    data = data.unique().sort(["genome", "contig", "gene_number"])
    return data

def main():
    input_file = snakemake.params.gene_index_annotated
    output_file = snakemake.params.context_table
    circular_contigs_file = snakemake.params.circular_contigs
    minimum_percentage = snakemake.params.annotation_percent_threshold
    window_size = snakemake.params.window_size
    minimum_vscore = snakemake.params.minimum_flank_vscore
    max_flank_length = snakemake.params.max_flank_length
    lgbm_model = load(snakemake.params.lgbm_model)
    feature_names = list(load(snakemake.params.feature_names))
    thresholds = load(snakemake.params.thresholds)
    outparent = snakemake.params.outparent
    n_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Starting genome context analysis...")
    os.makedirs(outparent, exist_ok=True)

    if not os.path.exists(circular_contigs_file) or os.path.getsize(circular_contigs_file) == 0:
        circular_contigs = set()
        logger.warning("No results found for checking circular contigs. All values for 'circular_contig' will be False.")
        logger.debug(f"Reading input file: {input_file}")
    else:
        logger.debug(f"Reading input files: {input_file} and {circular_contigs_file}")
        circular_contigs = set(pl.read_csv(circular_contigs_file, separator='\t')['contig'].to_list())

    data = pl.read_csv(input_file, separator='\t')
    logger.debug(f"Loaded data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")
    data = data.sort(["contig", "gene_number"]).unique()
    logger.debug(f"Unique data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")
    logger.debug(f"Data before processing: {data.head()}")

    use_hallmark = snakemake.params.use_hallmark
    hallmark_path = snakemake.params.hallmark_path
    hallmark_ids = None
    if use_hallmark:
        logger.debug(f"Reading hallmark file: {hallmark_path}")
        hallmark_data = pl.read_csv(hallmark_path)
        hallmark_ids = set(hallmark_data['id'])

    mobile_genes_path = snakemake.params.mobile_genes_path
    mobile_ids = None
    if mobile_genes_path:
        logger.debug(f"Reading MGE file: {mobile_genes_path}")
        mobile_genes_data = pl.read_csv(mobile_genes_path)
        mobile_ids = set(mobile_genes_data['id'])

    processed_data = process_genomes(
        data, circular_contigs, minimum_percentage,
        window_size, max_flank_length, minimum_vscore,
        lgbm_model, thresholds, feature_names,
        use_hallmark, hallmark_ids, mobile_ids, n_cpus
    )

    logger.debug(f"Writing output file: {output_file}")
    processed_data.write_csv(output_file, separator='\t', include_header=True)
    logger.info("Genome context analysis completed.")

if __name__ == "__main__":
    main()
