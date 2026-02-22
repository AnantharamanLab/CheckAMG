#!/usr/bin/env python3

import os
import sys
import resource
from pathlib import Path
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)
os.environ["OMP_NUM_THREADS"] = str(snakemake.threads)
os.environ["OPENBLAS_NUM_THREADS"] = str(snakemake.threads)
os.environ["MKL_NUM_THREADS"] = str(snakemake.threads)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(snakemake.threads)
import polars as pl
import numpy as np
from joblib import load
from numba import njit
import warnings
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

INT_PLACEHOLDER = np.iinfo(np.int32).min

def _as_parquet_path_if_enabled(p, save_to_parquet: bool):
    p = Path(p)
    if save_to_parquet and p.suffix.lower() == ".tsv":
        return p.with_suffix(".parquet")
    return p

def calculate_gene_lengths(data):
    """
    Calculate gene lengths and protein lengths in amino acids.
    """
    data = data.with_columns([
        (pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1).alias('gene_length_bases'),
        ((pl.col('contig_pos_end') - pl.col('contig_pos_start') + 1) / 3).cast(pl.Int32).alias('prot_length_AAs')
    ])
    return data

def calculate_contig_end_distances(data: pl.DataFrame) -> pl.DataFrame:
    """
    For each gene/protein, compute distance to contig ends:
      - contig_left_dist_bases / contig_right_dist_bases in nucleotide bases
      - contig_left_dist_genes / contig_right_dist_genes in genes (0 if first/last gene/protein on contig)

    Notes:
      - Uses contig_pos_start/contig_pos_end for base distances.
      - If 'contig_length' exists, right-end distance uses it.
      - Otherwise, right-end distance uses max(contig_pos_end) per contig (last annotated base)
      - Detects 0-based vs 1-based coordinates using global min(contig_pos_start).
    """
    required = {"contig", "gene_number", "contig_pos_start", "contig_pos_end"}
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns for contig-end distances: {missing}")

    min_start_global = data.select(pl.min("contig_pos_start")).item()
    coord_offset = 0 if min_start_global == 0 else 1

    if "contig_length" in data.columns:
        contig_end_expr = pl.col("contig_length")
    else:
        contig_end_expr = pl.max("contig_pos_end").over("contig")

    idx = pl.int_range(0, pl.len()).over("contig").cast(pl.Int32)

    return (
        data.sort(["contig", "gene_number"])
        .with_columns([
            # base distances
            (pl.col("contig_pos_start") - pl.lit(coord_offset))
            .cast(pl.Int32)
            .alias("contig_left_end_dist"),

            (contig_end_expr - pl.col("contig_pos_end"))
            .cast(pl.Int32)
            .alias("contig_right_end_dist"),

            # gene/protein distances (0 if first/last protein on contig)
            idx.alias("contig_left_end_gene_dist"),
            (pl.len().over("contig").cast(pl.Int32) - 1 - idx).alias("contig_right_end_gene_dist"),
        ])
    )

def calculate_contig_statistics(data, circular_contigs):
    """
    Calculate contig average V-scores/VL-scores and assign a circular_contig flag.
    """
    stats = data.group_by("contig", maintain_order=True).agg([
        pl.col("KEGG_V-score").mean().alias("contig_avg_KEGG_V-score"),
        pl.col("Pfam_V-score").mean().alias("contig_avg_Pfam_V-score"),
        pl.col("PHROG_V-score").mean().alias("contig_avg_PHROG_V-score"),
        # pl.col("VOG_V-score").mean().alias("contig_avg_VOG_V-score"),
        pl.col("KEGG_VL-score").mean().alias("contig_avg_KEGG_VL-score"),
        pl.col("Pfam_VL-score").mean().alias("contig_avg_Pfam_VL-score"),
        pl.col("PHROG_VL-score").mean().alias("contig_avg_PHROG_VL-score"),
        # pl.col("VOG_VL-score").mean().alias("contig_avg_VOG_VL-score"),
    ])
    result = data.join(stats, on="contig")
    result = result.with_columns(pl.col("contig").is_in(circular_contigs).alias("circular_contig"))
    return result


# Prophage identification in disguise as a non-virus identification tool
def identify_viral_nonviral_regions(
    data: pl.DataFrame,
    minimum_vscore: float = 10.0,
    minimum_window_avg_vlscore: float = 3.0,
    max_gap_genes: int = 3,
    min_rough_run: int = 1,
    merge_max_gap_genes: int = 10,
    merge_max_bad_genes: int = 1,
    bridge_vl_margin: float = 0.25,
    max_snap_distance_genes: int = 10,
    merge_min_bridge_genes: int = 2,
    merge_min_bridge_frac: float = 0.5,
    merge_max_consecutive_nonbridge: int = 0,
    min_region_genes: int = 5,
    max_walkback_genes: int = 20,
    snap_max_consecutive_nonbridge: int = 5,
    snap_require_bridge: bool = True,
    use_hallmark: bool = False,
    hallmark_accessions=None,
    n_cpus: int = 1,
) -> pl.DataFrame:
    """
    Add per-gene region calls plus explicit step-tracking booleans.

    Outputs (new/renamed):
    - candidate_region_id (Int32 / null): region id after Step 1 candidate scan (rough, gap-filled).
    - refined_region_id (Int32 / null): region id after Step 2–4 refinement (V-score gated core + boundary ops).
    - viral_region_id (Int32 / null): final region id after Step 5 merging.

    - step1_seed_strict_wvl (bool): gene passes the strict seed rule (any window_avg_*_VL-score >= minimum_window_avg_vlscore).
    - step1_in_candidate_region (bool): gene lies inside any Step 1 candidate region (after allowing max_gap_genes gaps).
    - step2_passes_vscore_gate (bool): gene is V-score compatible (>=1 V-score present AND none of present V-scores < minimum_vscore),
        OR (if use_hallmark=True) gene has any hallmark annotation across KEGG/Pfam/PHROG hmm_id columns.
    - step2_in_refined_core (bool): gene lies in the Step 2 refined core region (before Step 3/4).
    - step3_bridge_ok (bool): gene is bridge-eligible based on window-avg VL (has any window_avg VL and
        max(window_avg VL) >= (minimum_window_avg_vlscore - bridge_vl_margin)).
    - step3_in_walkback_extended (bool): gene lies in the region after Step 3 walkback extension (before Step 4 snap).
    - step4_in_snapped_region (bool): gene lies in the region after Step 4 snap-to-strong-V
        OR (if use_hallmark=True) snap-to-hallmark.
    - step5_in_merged_region (bool): gene lies in the final merged region set (after Step 5).

    Step-by-step logic (per contig):
    Step 1 (candidate scan; strict window-avg VL seeds):
        - step1_seed_strict_wvl is true if any window_avg_{KEGG,Pfam,PHROG}_VL-score >= minimum_window_avg_vlscore.
        - Build candidate regions by connecting seed hits allowing gaps up to max_gap_genes.
        - Keep regions with at least min_rough_run seed hits.
        - candidate_region_id and step1_in_candidate_region are assigned here.

    Step 2 (refine each candidate to a V-score supported core):
        - step2_passes_vscore_gate is true if at least one V-score is present AND none of the present V-scores are < minimum_vscore.
        - If use_hallmark is enabled, step2_passes_vscore_gate is instead true if the gene has any hallmark annotation
          (any of KEGG_hmm_id/Pfam_hmm_id/PHROG_hmm_id in hallmark_accessions).
        - Within each candidate, find the leftmost and rightmost genes satisfying the gate to form a "core".
        - refined_region_id (pre-boundary ops) and step2_in_refined_core are assigned here.

    Step 3 (walkback boundary extension with hard distance limit):
        - step3_bridge_ok is defined as above (bridge-eligible by window-avg VL).
        - Extend each core boundary outward while step3_bridge_ok holds, but never more than max_walkback_genes on a side.
        - Require the walked shoulder to contain at least one anchor gene with any V-score >= minimum_vscore;
          if use_hallmark is enabled, require at least one hallmark anchor instead.
        - step3_in_walkback_extended is assigned here.

    Step 4 (optional snap to nearby strong-V genes with path constraints):
        - Only applies on a side that Step 3 actually extended.
        - Look outward <= max_snap_distance_genes for a "strong" gene with all three V-scores present and each >= minimum_vscore.
          If use_hallmark is enabled, "strong" is any hallmark-annotated gene.
        - Allow passing at most one "bad" gene (has any V-score present and at least one present V-score < minimum_vscore).
          If use_hallmark is enabled, "bad" genes are not used to block snapping.
        - Additional anti-overextension constraints:
            * The path between the current boundary and the strong gene cannot contain a long run of non-bridge genes
            (max run length <= snap_max_consecutive_nonbridge).
            * If snap_require_bridge is True, the path must include at least one bridge gene (step3_bridge_ok).
        - step4_in_snapped_region is assigned here.

    Step 5 (merge adjacent refined regions only if the gap is genuinely bridge-like):
        - Before merging, drop any refined regions shorter than min_region_genes.
        - Consider consecutive surviving refined regions on the contig.
        - If the inter-region gap (in genes) <= merge_max_gap_genes, merge only if ALL of the following hold:
            * The gap contains at least merge_min_bridge_genes bridge genes (step3_bridge_ok).
            * The fraction of bridge genes in the gap is >= merge_min_bridge_frac.
            * The longest consecutive run of non-bridge genes in the gap is <= merge_max_consecutive_nonbridge.
            * The number of "bad" genes in the gap is <= merge_max_bad_genes.
        - viral_region_id and step5_in_merged_region are assigned here.

    Minimum region size:
    - min_region_genes enforces the minimum number of genes a refined region must contain to be eligible for Step 5 merging
        and to be reported as a final viral_region_id.
    """
    required = {"contig", "gene_number", "KEGG_V-score", "Pfam_V-score", "PHROG_V-score"}
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise ValueError(f"Missing required columns for viral region features: {missing}")

    if use_hallmark:
        hmm_required = {"KEGG_hmm_id", "Pfam_hmm_id", "PHROG_hmm_id"}
        hmm_missing = [c for c in hmm_required if c not in data.columns]
        if hmm_missing:
            raise ValueError(
                f"use_hallmark=True requires hmm_id columns for hallmark gating, missing: {hmm_missing}"
            )
        if hallmark_accessions is None:
            raise ValueError("use_hallmark=True requires hallmark_accessions to be provided")

    wvl_cols = [
        "window_avg_KEGG_VL-score",
        "window_avg_Pfam_VL-score",
        "window_avg_PHROG_VL-score",
    ]
    missing_wvl = [c for c in wvl_cols if c not in data.columns]
    if missing_wvl:
        raise ValueError(
            "Missing required window-average VL-score columns for candidate region scan: "
            f"{missing_wvl}"
        )

    if minimum_window_avg_vlscore is None:
        minimum_window_avg_vlscore = minimum_vscore

    if max_gap_genes < 0:
        raise ValueError("max_gap_genes must be >= 0")
    if min_rough_run < 1:
        raise ValueError("min_rough_run must be >= 1")
    if max_snap_distance_genes < 0:
        raise ValueError("max_snap_distance_genes must be >= 0")
    if merge_max_gap_genes < 0:
        raise ValueError("merge_max_gap_genes must be >= 0")
    if merge_max_bad_genes < 0:
        raise ValueError("merge_max_bad_genes must be >= 0")
    if bridge_vl_margin < 0:
        raise ValueError("bridge_vl_margin must be >= 0")
    if merge_min_bridge_genes < 0:
        raise ValueError("merge_min_bridge_genes must be >= 0")
    if not (0.0 <= merge_min_bridge_frac <= 1.0):
        raise ValueError("merge_min_bridge_frac must be between 0 and 1")
    if merge_max_consecutive_nonbridge < 0:
        raise ValueError("merge_max_consecutive_nonbridge must be >= 0")
    if min_region_genes < 1:
        raise ValueError("min_region_genes must be >= 1")
    if max_walkback_genes < 0:
        raise ValueError("max_walkback_genes must be >= 0")
    if snap_max_consecutive_nonbridge < 0:
        raise ValueError("snap_max_consecutive_nonbridge must be >= 0")

    group_cols = ["genome", "contig"] if "genome" in data.columns else ["contig"]

    data = data.with_columns(
        [
            pl.col("KEGG_V-score").cast(pl.Float64, strict=False),
            pl.col("Pfam_V-score").cast(pl.Float64, strict=False),
            pl.col("PHROG_V-score").cast(pl.Float64, strict=False),
            pl.col("window_avg_KEGG_VL-score").cast(pl.Float64, strict=False),
            pl.col("window_avg_Pfam_VL-score").cast(pl.Float64, strict=False),
            pl.col("window_avg_PHROG_VL-score").cast(pl.Float64, strict=False),
        ]
    )

    def _find_candidate_regions(seed: np.ndarray, max_gap: int, min_true_hits: int):
        n = seed.size
        regions = []
        i = 0
        while i < n:
            if not seed[i]:
                i += 1
                continue
            start = i
            true_hits = 0
            last_true = i
            gap = 0
            while i < n:
                if seed[i]:
                    true_hits += 1
                    last_true = i
                    gap = 0
                    i += 1
                    continue
                gap += 1
                if gap <= max_gap:
                    i += 1
                    continue
                break
            end = last_true
            if true_hits >= min_true_hits:
                regions.append((start, end))
        return regions

    def _max_consecutive_true(mask: np.ndarray) -> int:
        if mask.size == 0:
            return 0
        m = mask.astype(np.int8)
        dm = np.diff(np.r_[0, m, 0])
        starts = np.flatnonzero(dm == 1)
        ends = np.flatnonzero(dm == -1)
        if starts.size == 0:
            return 0
        return int((ends - starts).max())

    def _per_contig(df: pl.DataFrame) -> pl.DataFrame:
        df = df.sort("gene_number")
        n = df.height

        if n == 0:
            return df.with_columns(
                [
                    pl.lit(None, dtype=pl.Int32).alias("candidate_region_id"),
                    pl.lit(None, dtype=pl.Int32).alias("refined_region_id"),
                    pl.lit(None, dtype=pl.Int32).alias("viral_region_id"),
                    pl.lit(False).alias("step1_seed_strict_wvl"),
                    pl.lit(False).alias("step1_in_candidate_region"),
                    pl.lit(False).alias("step2_passes_vscore_gate"),
                    pl.lit(False).alias("step2_in_refined_core"),
                    pl.lit(False).alias("step3_bridge_ok"),
                    pl.lit(False).alias("step3_in_walkback_extended"),
                    pl.lit(False).alias("step4_in_snapped_region"),
                    pl.lit(False).alias("step5_in_merged_region"),
                ]
            )

        kegg_v = df["KEGG_V-score"].to_numpy()
        pfam_v = df["Pfam_V-score"].to_numpy()
        phrog_v = df["PHROG_V-score"].to_numpy()

        kegg_v_nan = np.isnan(kegg_v)
        pfam_v_nan = np.isnan(pfam_v)
        phrog_v_nan = np.isnan(phrog_v)

        any_v_present = (~kegg_v_nan) | (~pfam_v_nan) | (~phrog_v_nan)
        any_v_below = (
            ((~kegg_v_nan) & (kegg_v < minimum_vscore))
            | ((~pfam_v_nan) & (pfam_v < minimum_vscore))
            | ((~phrog_v_nan) & (phrog_v < minimum_vscore))
        )

        kegg_wvl = df["window_avg_KEGG_VL-score"].to_numpy()
        pfam_wvl = df["window_avg_Pfam_VL-score"].to_numpy()
        phrog_wvl = df["window_avg_PHROG_VL-score"].to_numpy()

        wvl_stack = np.vstack([kegg_wvl, pfam_wvl, phrog_wvl]).astype(float)
        wvl_finite = np.isfinite(wvl_stack)
        wvl_any = wvl_finite.any(axis=0)
        wvl_tmp = np.where(wvl_finite, wvl_stack, -np.inf)
        wvl_max = wvl_tmp.max(axis=0)
        wvl_max = np.where(wvl_any, wvl_max, -np.inf)

        strict_thr = float(minimum_window_avg_vlscore)
        bridge_thr = float(minimum_window_avg_vlscore) - float(bridge_vl_margin)

        step1_seed_strict_wvl = wvl_any & (wvl_max >= strict_thr)
        step3_bridge_ok = wvl_any & (wvl_max >= bridge_thr)

        if use_hallmark:
            hallmark_set = hallmark_accessions if isinstance(hallmark_accessions, set) else set(hallmark_accessions)
            kegg_hmm = df["KEGG_hmm_id"].to_list()
            pfam_hmm = df["Pfam_hmm_id"].to_list()
            phrog_hmm = df["PHROG_hmm_id"].to_list()

            hallmark_any = np.array(
                [
                    ((k is not None) and (k in hallmark_set))
                    or ((p is not None) and (p in hallmark_set))
                    or ((h is not None) and (h in hallmark_set))
                    for k, p, h in zip(kegg_hmm, pfam_hmm, phrog_hmm)
                ],
                dtype=np.bool_,
            )

            step2_passes_vscore_gate = hallmark_any
            bad_gene = np.zeros(n, dtype=np.bool_)  # do not block snapping based on V-score thresholds in hallmark mode
            v_anchor = hallmark_any
            all3_ge_min = hallmark_any  # "strong" target becomes any hallmark-annotated gene
        else:
            step2_passes_vscore_gate = any_v_present & (~any_v_below)
            bad_gene = any_v_present & any_v_below

            v_stack = np.vstack([kegg_v, pfam_v, phrog_v]).astype(float)
            v_finite = np.isfinite(v_stack)
            v_any = v_finite.any(axis=0)
            v_tmp = np.where(v_finite, v_stack, -np.inf)
            v_max = v_tmp.max(axis=0)
            v_max = np.where(v_any, v_max, -np.inf)
            v_anchor = v_any & (v_max >= float(minimum_vscore))

            all3_present = (~kegg_v_nan) & (~pfam_v_nan) & (~phrog_v_nan)
            all3_ge_min = all3_present & (kegg_v >= minimum_vscore) & (pfam_v >= minimum_vscore) & (phrog_v >= minimum_vscore)

        candidate_regions = _find_candidate_regions(step1_seed_strict_wvl, max_gap_genes, min_rough_run)
        candidate_id = np.full(n, -1, dtype=np.int32)
        step1_in_candidate = np.zeros(n, dtype=bool)
        for rid, (a, b) in enumerate(candidate_regions):
            candidate_id[a:b + 1] = rid
            step1_in_candidate[a:b + 1] = True

        refined_regions_all = []
        refined_id = np.full(n, -1, dtype=np.int32)

        step2_in_refined_core = np.zeros(n, dtype=bool)
        step3_in_walkback = np.zeros(n, dtype=bool)
        step4_in_snapped = np.zeros(n, dtype=bool)

        def _walkback_left(ref_left: int) -> tuple[int, bool]:
            i = ref_left - 1
            if i < 0 or max_walkback_genes == 0:
                return ref_left, False
            L_ext = ref_left
            walked = 0
            while i >= 0 and step3_bridge_ok[i] and walked < max_walkback_genes:
                L_ext = i
                i -= 1
                walked += 1
            cand = np.flatnonzero(v_anchor[L_ext:ref_left + 1])
            if cand.size == 0:
                return ref_left, False
            new_left = L_ext + int(cand[0])
            return new_left, (new_left < ref_left)

        def _walkback_right(ref_right: int) -> tuple[int, bool]:
            i = ref_right + 1
            if i >= n or max_walkback_genes == 0:
                return ref_right, False
            R_ext = ref_right
            walked = 0
            while i < n and step3_bridge_ok[i] and walked < max_walkback_genes:
                R_ext = i
                i += 1
                walked += 1
            cand = np.flatnonzero(v_anchor[ref_right:R_ext + 1])
            if cand.size == 0:
                return ref_right, False
            new_right = ref_right + int(cand[-1])
            return new_right, (new_right > ref_right)

        def _snap_left_if_strong(new_left: int, did_extend: bool) -> int:
            if (not did_extend) or max_snap_distance_genes == 0:
                return new_left
            end = max(0, new_left - max_snap_distance_genes)
            bad_seen = 0
            best = None
            for j in range(new_left - 1, end - 1, -1):
                if bad_gene[j]:
                    bad_seen += 1
                    if bad_seen > 1:
                        break
                if all3_ge_min[j]:
                    path = step3_bridge_ok[j:new_left]
                    if snap_require_bridge and int(path.sum()) == 0:
                        continue
                    if _max_consecutive_true(~path) > int(snap_max_consecutive_nonbridge):
                        continue
                    best = j
            return best if best is not None else new_left

        def _snap_right_if_strong(new_right: int, did_extend: bool) -> int:
            if (not did_extend) or max_snap_distance_genes == 0:
                return new_right
            end = min(n - 1, new_right + max_snap_distance_genes)
            bad_seen = 0
            best = None
            for j in range(new_right + 1, end + 1):
                if bad_gene[j]:
                    bad_seen += 1
                    if bad_seen > 1:
                        break
                if all3_ge_min[j]:
                    path = step3_bridge_ok[new_right + 1:j + 1]
                    if snap_require_bridge and int(path.sum()) == 0:
                        continue
                    if _max_consecutive_true(~path) > int(snap_max_consecutive_nonbridge):
                        continue
                    best = j
            return best if best is not None else new_right

        for rid, (a, b) in enumerate(candidate_regions):
            core_left = None
            for i in range(a, b + 1):
                if step2_passes_vscore_gate[i]:
                    core_left = i
                    break
            core_right = None
            for i in range(b, a - 1, -1):
                if step2_passes_vscore_gate[i]:
                    core_right = i
                    break
            if core_left is None or core_right is None or core_left > core_right:
                continue

            step2_in_refined_core[core_left:core_right + 1] = True

            new_left, did_extend_left = _walkback_left(core_left)
            new_right, did_extend_right = _walkback_right(core_right)
            step3_in_walkback[new_left:new_right + 1] = True

            new_left = _snap_left_if_strong(new_left, did_extend_left)
            new_right = _snap_right_if_strong(new_right, did_extend_right)
            step4_in_snapped[new_left:new_right + 1] = True

            refined_regions_all.append((new_left, new_right))

        for rrid, (L, R) in enumerate(refined_regions_all):
            refined_id[L:R + 1] = rrid

        refined_regions = [
            (L, R) for (L, R) in refined_regions_all
            if (R - L + 1) >= min_region_genes
        ]

        merged_regions = []
        if refined_regions:
            refined_regions_sorted = sorted(refined_regions, key=lambda x: (x[0], x[1]))
            curL, curR = refined_regions_sorted[0]

            for nxtL, nxtR in refined_regions_sorted[1:]:
                gapL = curR + 1
                gapR = nxtL - 1
                gap_len = gapR - gapL + 1

                can_merge = False
                if gap_len <= 0:
                    can_merge = True
                elif gap_len <= merge_max_gap_genes:
                    gap_bridge = step3_bridge_ok[gapL:gapR + 1]
                    gap_bad = bad_gene[gapL:gapR + 1]

                    bridge_count = int(gap_bridge.sum())
                    bridge_frac = bridge_count / float(gap_len)

                    nonbridge = ~gap_bridge
                    max_dead_run = _max_consecutive_true(nonbridge)

                    if (
                        bridge_count >= int(merge_min_bridge_genes)
                        and bridge_frac >= float(merge_min_bridge_frac)
                        and max_dead_run <= int(merge_max_consecutive_nonbridge)
                        and int(gap_bad.sum()) <= int(merge_max_bad_genes)
                    ):
                        can_merge = True

                if can_merge:
                    curR = max(curR, nxtR)
                else:
                    merged_regions.append((curL, curR))
                    curL, curR = nxtL, nxtR

            merged_regions.append((curL, curR))

        final_id = np.full(n, -1, dtype=np.int32)
        step5_in_merged = np.zeros(n, dtype=bool)
        for vid, (L, R) in enumerate(merged_regions):
            final_id[L:R + 1] = vid
            step5_in_merged[L:R + 1] = True

        def _to_nullable_int32(arr: np.ndarray) -> list:
            return [None if int(x) < 0 else int(x) for x in arr.tolist()]

        return df.with_columns(
            [
                pl.Series("candidate_region_id", _to_nullable_int32(candidate_id), dtype=pl.Int32),
                pl.Series("refined_region_id", _to_nullable_int32(refined_id), dtype=pl.Int32),
                pl.Series("viral_region_id", _to_nullable_int32(final_id), dtype=pl.Int32),
                pl.Series("step1_seed_strict_wvl", step1_seed_strict_wvl.tolist()),
                pl.Series("step1_in_candidate_region", step1_in_candidate.tolist()),
                pl.Series("step2_passes_vscore_gate", step2_passes_vscore_gate.tolist()),
                pl.Series("step2_in_refined_core", step2_in_refined_core.tolist()),
                pl.Series("step3_bridge_ok", step3_bridge_ok.tolist()),
                pl.Series("step3_in_walkback_extended", step3_in_walkback.tolist()),
                pl.Series("step4_in_snapped_region", step4_in_snapped.tolist()),
                pl.Series("step5_in_merged_region", step5_in_merged.tolist()),
            ]
        )

    contig_dfs = data.partition_by(group_cols, maintain_order=True)

    if n_cpus is None or n_cpus < 2 or len(contig_dfs) < 2:
        out = [_per_contig(df) for df in contig_dfs]
        return pl.concat(out, how="vertical")

    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [executor.submit(_per_contig, df) for df in contig_dfs]
        out = [f.result() for f in as_completed(futures)]

    return pl.concat(out, how="vertical")

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

def process_window_statistics_for_contig_df(df: pl.DataFrame, window_size, minimum_percentage) -> pl.DataFrame:
    df = df.sort("gene_number")
    lengths = df["gene_length_bases"].to_numpy()
    if len(lengths) == 0:
        return df

    kegg_vl = df["KEGG_VL-score"].to_numpy()
    pfam_vl = df["Pfam_VL-score"].to_numpy()
    phrog_vl = df["PHROG_VL-score"].to_numpy()
    # vog_vl = df["VOG_VL-score"].to_numpy()

    kegg_v = df["KEGG_V-score"].to_numpy()
    pfam_v = df["Pfam_V-score"].to_numpy()
    phrog_v = df["PHROG_V-score"].to_numpy()
    # vog_v = df["VOG_V-score"].to_numpy()

    df = df.with_columns([
        pl.Series("window_avg_KEGG_VL-score", window_avg(kegg_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_Pfam_VL-score", window_avg(pfam_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_PHROG_VL-score", window_avg(phrog_vl, lengths, window_size, minimum_percentage)),
        # pl.Series("window_avg_VOG_VL-score", window_avg(vog_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_KEGG_V-score", window_avg(kegg_v, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_Pfam_V-score", window_avg(pfam_v, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_PHROG_V-score", window_avg(phrog_v, lengths, window_size, minimum_percentage)),
        # pl.Series("window_avg_VOG_V-score", window_avg(vog_v, lengths, window_size, minimum_percentage)),
    ])
    return df

def calculate_window_statistics(data, window_size, minimum_percentage, n_cpus):
    data = data.sort(["contig", "gene_number"])
    contig_dfs = data.partition_by("contig", maintain_order=True)

    if n_cpus is None or n_cpus < 2 or len(contig_dfs) < 2:
        out = [
            process_window_statistics_for_contig_df(df, window_size, minimum_percentage)
            for df in tqdm(contig_dfs, total=len(contig_dfs), desc="Calculating sliding-window averages", unit="contig")
        ]
        return pl.concat(out, how="vertical")

    with ThreadPoolExecutor(max_workers=n_cpus) as ex:
        futs = [ex.submit(process_window_statistics_for_contig_df, df, window_size, minimum_percentage) for df in contig_dfs]
        out = [f.result() for f in tqdm(as_completed(futs), total=len(futs), desc="Calculating sliding-window averages", unit="contig")]
    return pl.concat(out, how="vertical")

def viral_origin_confidence_lgbm(df, lgbm_model, thresholds, feature_names, n_threads=1, batch_size=200_000):
    """
    df: pl.DataFrame (must have columns matching those used in model training)
    lgbm_model: fitted sklearn Light GBM or CalibratedClassifierCV
    thresholds: dict, { 'high': {'threshold': float, ...}, 'medium': {'threshold': float, ...}, 'low': {'threshold': float, ...} }
    feature_names: list of feature columns (ordered)
    Returns polars DataFrame with added columns 'LGBM_viral_prob' and 'Viral_Origin_Confidence'
    """
    # Ensure LightGBM is allowed to use threads
    try:
        lgbm_model.set_params(n_jobs=int(n_threads))
    except Exception:
        pass
    try:
        booster = getattr(lgbm_model, "booster_", None)
        if booster is not None:
            booster.params["num_threads"] = int(n_threads)
    except Exception:
        pass

    # Ensure required features are present, fill missing with NaN
    missing = [c for c in feature_names if c not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(np.nan, dtype=pl.Float32).alias(c) for c in missing])

    feat_df = df.select([pl.col(c).cast(pl.Float32, strict=False) for c in feature_names])
    X_all = feat_df.to_numpy()
    n = X_all.shape[0]

    y_proba_all = np.empty(n, dtype=np.float32)

    # Silence the sklearn "feature names" warning for this call site
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"X does not have valid feature names, but LGBMClassifier was fitted with feature names",
            category=UserWarning,
            module=r"sklearn\.utils\.validation",
        )

        if batch_size is None or batch_size <= 0 or n <= batch_size:
            y_proba_all[:] = lgbm_model.predict_proba(X_all)[:, 1].astype(np.float32)
        else:
            n_batches = (n + batch_size - 1) // batch_size
            for b in tqdm(range(n_batches), desc="Fitting LGBM", unit="batch"):
                start = b * batch_size
                end = min(n, start + batch_size)
                y_proba_all[start:end] = lgbm_model.predict_proba(X_all[start:end])[:, 1].astype(np.float32)

    conf = np.full(n, 'low', dtype=object)
    conf[y_proba_all >= thresholds['medium']['threshold']] = 'medium'
    conf[y_proba_all >= thresholds['high']['threshold']] = 'high'

    df = df.with_columns([
        pl.Series('LGBM_viral_prob', y_proba_all),
        pl.Series('Viral_Origin_Confidence', conf)
    ])

    return df

@njit
def flank_distance_vscores(lengths, scores, minimum_vscore):
    n = len(lengths)
    positions = np.zeros(n, dtype=np.float64)
    for i in range(1, n):
        positions[i] = positions[i-1] + lengths[i-1]

    left_dist = np.full(n, INT_PLACEHOLDER, dtype=np.float64)
    right_dist = np.full(n, INT_PLACEHOLDER, dtype=np.float64)

    last_idx = -1
    for i in range(n):
        if last_idx != -1:
            left_dist[i] = positions[i] - positions[last_idx]
        if not np.isnan(scores[i]) and scores[i] >= minimum_vscore:
            last_idx = i

    next_idx = -1
    for i in range(n - 1, -1, -1):
        if next_idx != -1:
            right_dist[i] = positions[next_idx] - positions[i]
        if not np.isnan(scores[i]) and scores[i] >= minimum_vscore:
            next_idx = i

    return left_dist, right_dist

def verify_flanking_vscores(contig_data, minimum_vscore):
    """
    For each gene, calculate nucleotide distance to nearest left and right gene meeting V-score threshold.
    Adds columns: KEGG_viral_left_dist, KEGG_viral_right_dist, etc.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_scores = contig_data['KEGG_V-score'].to_numpy()
    pfam_scores = contig_data['Pfam_V-score'].to_numpy()
    phrog_scores = contig_data['PHROG_V-score'].to_numpy()
    # vog_scores = contig_data['VOG_V-score'].to_numpy()

    kegg_left, kegg_right = flank_distance_vscores(lengths, kegg_scores, minimum_vscore)
    pfam_left, pfam_right = flank_distance_vscores(lengths, pfam_scores, minimum_vscore)
    phrog_left, phrog_right = flank_distance_vscores(lengths, phrog_scores, minimum_vscore)
    # vog_left, vog_right = flank_distance_vscores(lengths, vog_scores, minimum_vscore)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str]*len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_viral_left_dist": kegg_left.astype(np.int32),
        "KEGG_viral_right_dist": kegg_right.astype(np.int32),
        "Pfam_viral_left_dist": pfam_left.astype(np.int32),
        "Pfam_viral_right_dist": pfam_right.astype(np.int32),
        "PHROG_viral_left_dist": phrog_left.astype(np.int32),
        "PHROG_viral_right_dist": phrog_right.astype(np.int32),
        # "VOG_viral_left_dist": vog_left.astype(np.int32),
        # "VOG_viral_right_dist": vog_right.astype(np.int32),
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

@njit
def flank_distance_in_set(lengths, in_set):
    n = len(lengths)
    positions = np.zeros(n, dtype=np.float64)
    for i in range(1, n):
        positions[i] = positions[i-1] + lengths[i-1]

    left_dist = np.full(n, INT_PLACEHOLDER, dtype=np.float64)
    right_dist = np.full(n, INT_PLACEHOLDER, dtype=np.float64)

    last_idx = -1
    for i in range(n):
        if last_idx != -1:
            left_dist[i] = positions[i] - positions[last_idx]
        if in_set[i]:
            last_idx = i

    next_idx = -1
    for i in range(n - 1, -1, -1):
        if next_idx != -1:
            right_dist[i] = positions[next_idx] - positions[i]
        if in_set[i]:
            next_idx = i

    return left_dist, right_dist

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

def verify_flanking_hallmark(contig_data, hallmark_accessions):
    """
    For each gene, calculate nucleotide distance to nearest left/right hallmark gene.
    Adds columns: KEGG_hallmark_left_dist, KEGG_hallmark_right_dist, etc.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()
    # vog_hmm = contig_data['VOG_hmm_id'].to_list()

    kegg_arr = create_in_set_array(kegg_hmm, hallmark_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, hallmark_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, hallmark_accessions)
    # vog_arr = create_in_set_array(vog_hmm, hallmark_accessions)

    kegg_left, kegg_right = flank_distance_in_set(lengths, kegg_arr)
    pfam_left, pfam_right = flank_distance_in_set(lengths, pfam_arr)
    phrog_left, phrog_right = flank_distance_in_set(lengths, phrog_arr)
    # vog_left, vog_right = flank_distance_in_set(lengths, vog_arr)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str]*len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_viral_left_dist": kegg_left.astype(np.int32),
        "KEGG_viral_right_dist": kegg_right.astype(np.int32),
        "Pfam_viral_left_dist": pfam_left.astype(np.int32),
        "Pfam_viral_right_dist": pfam_right.astype(np.int32),
        "PHROG_viral_left_dist": phrog_left.astype(np.int32),
        "PHROG_viral_right_dist": phrog_right.astype(np.int32),
        # "VOG_viral_left_dist": vog_left.astype(np.int32),
        # "VOG_viral_right_dist": vog_right.astype(np.int32),
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

@njit
def flank_nearest_mge_scores(scores, mge_mask):
    n = len(scores)
    left = np.full(n, np.nan, dtype=np.float32)
    right = np.full(n, np.nan, dtype=np.float32)

    last_val = np.nan
    for i in range(n):
        left[i] = last_val
        if mge_mask[i] and not np.isnan(scores[i]):
            last_val = np.float32(scores[i])

    next_val = np.nan
    for i in range(n - 1, -1, -1):
        right[i] = next_val
        if mge_mask[i] and not np.isnan(scores[i]):
            next_val = np.float32(scores[i])

    return left, right

def report_flanking_mge_vscores(contig_data, mobile_accessions):
    """
    For each gene, report the KEGG/Pfam/PHROG V/VL-score of the nearest left/right MGE gene.
    Adds columns: <DB>_V-score_left_MGE, <DB>_V-score_right_MGE, <DB>_VL-score_left_MGE, etc.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    add_cols = {}
    for db in ["KEGG", "Pfam", "PHROG"]:
    # for db in ["KEGG", "Pfam", "PHROG", "VOG"]:
        hmm_col = f"{db}_hmm_id"
        mge_mask = np.array([x in mobile_accessions if x is not None else False for x in contig_data[hmm_col].to_list()], dtype=np.bool_)
        for sfx in ["V-score", "VL-score"]:
            scores = contig_data[f"{db}_{sfx}"].to_numpy()
            left, right = flank_nearest_mge_scores(scores, mge_mask)
            add_cols[f"{db}_{sfx}_left_MGE"] = left
            add_cols[f"{db}_{sfx}_right_MGE"] = right
    add_df = pl.DataFrame({
        "contig": contig_data["contig"],
        "gene_number": contig_data["gene_number"],
        **add_cols
    })
    return contig_data.join(add_df, on=["contig", "gene_number"])

def check_flanking_mge(contig_data, mobile_accessions):
    """
    For each gene, calculate nucleotide distance to nearest left/right MGE gene.
    Adds columns: KEGG_MGE_left_dist, KEGG_MGE_right_dist, etc.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()
    # vog_hmm = contig_data['VOG_hmm_id'].to_list()

    kegg_arr = create_in_set_array(kegg_hmm, mobile_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, mobile_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, mobile_accessions)
    # vog_arr = create_in_set_array(vog_hmm, mobile_accessions)

    kegg_left, kegg_right = flank_distance_in_set(lengths, kegg_arr)
    pfam_left, pfam_right = flank_distance_in_set(lengths, pfam_arr)
    phrog_left, phrog_right = flank_distance_in_set(lengths, phrog_arr)
    # vog_left, vog_right = flank_distance_in_set(lengths, vog_arr)

    contig_str = contig_data["contig"].to_list()[0]
    flanks_df = pl.DataFrame({
        "contig": [contig_str] * len(lengths),
        "gene_number": contig_data["gene_number"].to_list(),
        "KEGG_MGE_left_dist": kegg_left.astype(np.int32),
        "KEGG_MGE_right_dist": kegg_right.astype(np.int32),
        "Pfam_MGE_left_dist": pfam_left.astype(np.int32),
        "Pfam_MGE_right_dist": pfam_right.astype(np.int32),
        "PHROG_MGE_left_dist": phrog_left.astype(np.int32),
        "PHROG_MGE_right_dist": phrog_right.astype(np.int32),
        # "VOG_MGE_left_dist": vog_left.astype(np.int32),
        # "VOG_MGE_right_dist": vog_right.astype(np.int32),
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])

def process_genomes(data,
                    circular_contigs, minimum_percentage,
                    window_size,
                    lgbm_model, thresholds, feature_names,
                    minimum_vscore=10.0,
                    minimum_window_avg_vlscore=3.0,
                    max_gap_genes=3, min_rough_run=1,
                    max_snap_distance_genes=10,
                    use_hallmark=False,
                    hallmark_accessions=None, mobile_accessions=None,
                    n_cpus=1, mem_limit=10):
    logger.debug(f"Calculating lengths for {data.shape[0]:,} genes.")
    logger.debug(f"Data before calculating gene lengths: {data.head()}")
    data = calculate_gene_lengths(data)
    data_orig = data

    logger.info("Calculating distance to contig ends.")
    logger.debug(f"Data before calculating contig-end distances: {data.head()}")
    data = calculate_contig_end_distances(data)

    logger.info("Calculating contig statistics.")
    logger.debug(f"Data before calculating contig statistics: {data.head()}")
    data = calculate_contig_statistics(data, circular_contigs)

    logger.info("Calculating window statistics.")
    logger.debug(f"Data before calculating window statistics: {data.head()}")
    logger.debug(f"Column dtypes before conversion: {data.schema}")
    score_columns = [
        "KEGG_V-score","KEGG_VL-score","Pfam_V-score","Pfam_VL-score","PHROG_V-score","PHROG_VL-score",#"VOG_V-score","VOG_VL-score",
        "contig_avg_KEGG_V-score","contig_avg_Pfam_V-score","contig_avg_PHROG_V-score",#"contig_avg_VOG_V-score",
        "contig_avg_KEGG_VL-score","contig_avg_Pfam_VL-score", "contig_avg_PHROG_VL-score",#"contig_avg_VOG_VL-score"
    ]
    for col in score_columns:
        if col in data.columns:
            data = data.with_columns(pl.col(col).cast(pl.Float64, strict=False))
    logger.debug(f"Column dtypes after conversion: {data.schema}")

    # Parallel window statistics calculated per contig.
    data = calculate_window_statistics(data, window_size, minimum_percentage, n_cpus)
    data = data.unique()

    logger.info("Identifying ambiguous viral/nonviral contig regions.")
    logger.debug(f"Data before calculating nonviral region features: {data.head()}")
    data =  identify_viral_nonviral_regions(
        data=data,
        minimum_vscore=minimum_vscore,
        minimum_window_avg_vlscore=minimum_window_avg_vlscore,
        max_gap_genes=max_gap_genes,
        min_rough_run=min_rough_run,
        max_snap_distance_genes=max_snap_distance_genes,
        use_hallmark=use_hallmark,
        hallmark_accessions=hallmark_accessions,
        n_cpus=n_cpus
    )

    # Parallel verification of flanking regions by partitioning by contig.
    if use_hallmark and hallmark_accessions is not None:
        logger.info("Calculating distance of nearest flanking hallmark genes.")
        logger.debug(f"Data before verifying flanking hallmark genes: {data.head()}")
        contig_dfs = data.partition_by("contig")
        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            futures = [
                executor.submit(verify_flanking_hallmark, df, hallmark_accessions)
                for df in contig_dfs
            ]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanks for viral hallmarks", unit="contig")]
        data = pl.concat(results, how="vertical")
    else:
        logger.info(f"Calculating distance of nearest flanking V-score={minimum_vscore} genes.")
        logger.debug(f"Data before verifying flanking V-scores: {data.head()}")
        contig_dfs = data.partition_by("contig")
        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            futures = [
                executor.submit(verify_flanking_vscores, df, minimum_vscore)
                for df in contig_dfs
            ]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc=f"Checking flanks for V-score={minimum_vscore}", unit="contig")]
        data = pl.concat(results, how="vertical")

    logger.info("Calculating distance of nearest mobile genetic element genes.")
    logger.debug(f"Data before checking for mobile genetic element genes: {data.head()}")

    contig_dfs = data_orig.partition_by("contig")
    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(check_flanking_mge, df, mobile_accessions)
            for df in contig_dfs
        ]
        results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc=f"Checking flanks for mobile genes", unit="contig")]

    # Identify distance columns in the concatenated results
    dist_df = pl.concat(results, how="vertical")
    dist_cols = [c for c in dist_df.columns if c.endswith("_dist")]

    if not dist_cols:
        raise ValueError("No flanking MGE distance columns found, cannot merge.")

    # Sub-frame with only the protein identifier and distance columns
    dist_df = dist_df.select(["protein", *dist_cols])

    # Verify all distance columns have the same length
    heights = {dist_df[c].len() for c in dist_cols}
    if len(heights) > 1:
        raise ValueError(
            f"Flanking MGE distance columns have different heights: {sorted(heights)}, cannot merge."
        )

    dist_height = heights.pop()
    if dist_height != data.height:
        raise ValueError(
            f"Flanking MGE distance columns have height {dist_height}, "
            f"but original data has {data.height}, cannot merge."
        )

    # Append distance columns to the original DataFrame
    data = data.hstack([dist_df[c] for c in dist_cols])

    dist_cols = [col for col in data.columns if col.endswith('_dist')]
    logger.debug(f"Columns that will be recast to float32: {dist_cols}")
    logger.debug(f"Data before recasting distance columns: {data.head()}")
    for col in dist_cols:
        data = data.with_columns(
            pl.when(pl.col(col) == INT_PLACEHOLDER)
            .then(np.nan)
            .otherwise(pl.col(col))
            .alias(col).cast(pl.Float32)
        )
    logger.debug(f"Data after recasting distance columns: {data.head()}")

    logger.info("Calculating flanking MGE V/VL-scores.")
    logger.debug(f"Data before reporting flanking scores: {data.head()}")

    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [executor.submit(report_flanking_mge_vscores, df, mobile_accessions) for df in contig_dfs]
        mge_score_results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanking MGE V/VL-scores", unit="contig")]

    mge_score_df = pl.concat(mge_score_results, how="vertical")

    # Now join these to the `data` frame:
    data = data.join(mge_score_df.select([c for c in mge_score_df.columns if c not in data.columns or c in ["contig", "gene_number"]]), on=["contig", "gene_number"], how="left")

    logger.info("Assigning viral origin confidence using LightGBM with genome context features.")
    logger.debug(f"Data before assigning viral origin confidence: {data.head()}")
    data = viral_origin_confidence_lgbm(data, lgbm_model, thresholds, feature_names, n_threads=n_cpus)

    data = data.unique().sort(["genome", "contig", "gene_number"])
    return data

def main():
    minimum_percentage = snakemake.params.annotation_percent_threshold
    window_size = snakemake.params.window_size
    minimum_vscore = snakemake.params.minimum_flank_vscore
    minimum_window_avg_vlscore = snakemake.params.minimum_window_avg_vlscore
    max_gap_genes = snakemake.params.max_gap_genes
    min_rough_run = snakemake.params.min_rough_run
    max_snap_distance_genes = snakemake.params.max_snap_distance_genes
    lgbm_model = load(snakemake.params.lgbm_model)
    feature_names = list(load(snakemake.params.feature_names))
    thresholds = load(snakemake.params.thresholds)
    outparent = snakemake.params.outparent
    save_to_parquet = snakemake.params.save_to_parquet
    n_cpus = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)

    logger.info("Starting genome context analysis...")
    os.makedirs(outparent, exist_ok=True)

    input_file = _as_parquet_path_if_enabled(snakemake.params.gene_index_annotated, save_to_parquet)
    output_file = _as_parquet_path_if_enabled(snakemake.params.context_table, save_to_parquet)
    circular_contigs_file = _as_parquet_path_if_enabled(snakemake.params.circular_contigs, save_to_parquet)

    if not os.path.exists(circular_contigs_file) or os.path.getsize(circular_contigs_file) == 0:
        circular_contigs = set()
        logger.warning("No results found for checking circular contigs. All values for 'circular_contig' will be False.")
        logger.debug(f"Reading input file: {input_file}")
        if save_to_parquet:
            data = pl.read_parquet(input_file)
        else:
            data = pl.read_csv(input_file, separator='\t')
    else:
        logger.debug(f"Reading input files: {input_file} and {circular_contigs_file}")
        if save_to_parquet:
            circular_contigs = set(pl.read_parquet(circular_contigs_file)['contig'].to_list())
            data = pl.read_parquet(input_file)
        else:
            circular_contigs = set(pl.read_csv(circular_contigs_file, separator='\t')['contig'].to_list())
            data = pl.read_csv(input_file, separator='\t')

    logger.debug(f"Loaded data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")
    data = data.sort(["contig", "gene_number"]).unique()
    logger.debug(f"Unique data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")
    logger.debug(f"Data before processing: {data.head()}")

    # Reference hallmark and mobile genes are always TSV, not parquet
    use_hallmark = snakemake.params.use_hallmark
    hallmark_path = snakemake.params.hallmark_path
    hallmark_ids = None
    if use_hallmark:
        logger.debug(f"Reading hallmark file: {hallmark_path}")
        hallmark_data = pl.read_csv(hallmark_path, separator='\t')
        hallmark_ids = set(hallmark_data['id'])

    mobile_genes_path = snakemake.params.mobile_genes_path
    mobile_ids = None
    if mobile_genes_path:
        logger.debug(f"Reading MGE file: {mobile_genes_path}")
        mobile_genes_data = pl.read_csv(mobile_genes_path, separator='\t')
        mobile_ids = set(mobile_genes_data['id'])

    processed_data = process_genomes(
        data,
        circular_contigs, minimum_percentage,
        window_size,
        lgbm_model, thresholds, feature_names,
        minimum_vscore,
        minimum_window_avg_vlscore,
        max_gap_genes, min_rough_run,
        max_snap_distance_genes,
        use_hallmark, hallmark_ids, mobile_ids,
        n_cpus, mem_limit
    )

    logger.debug(f"Writing output file: {output_file}")
    if save_to_parquet:
        processed_data.write_parquet(output_file)
    else:
        processed_data.write_csv(output_file, separator='\t', include_header=True)
    logger.info("Genome context analysis completed.")

if __name__ == "__main__":
    main()