#!/usr/bin/env python3

import os
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

warnings.filterwarnings("ignore", category=RuntimeWarning, message="Mean of empty slice")

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit
from CheckAMG.scripts.annotate.utils import as_parquet_path_if_enabled

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Analyze the genomic context of annotations", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

INT_PLACEHOLDER = np.iinfo(np.int32).min


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
    For each gene/protein, compute distance to contig ends.
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
            (pl.col("contig_pos_start") - pl.lit(coord_offset))
            .cast(pl.Int32)
            .alias("contig_left_end_dist"),

            (contig_end_expr - pl.col("contig_pos_end"))
            .cast(pl.Int32)
            .alias("contig_right_end_dist"),

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
        pl.col("KEGG_VL-score").mean().alias("contig_avg_KEGG_VL-score"),
        pl.col("Pfam_VL-score").mean().alias("contig_avg_Pfam_VL-score"),
        pl.col("PHROG_VL-score").mean().alias("contig_avg_PHROG_VL-score"),
    ])
    result = data.join(stats, on="contig")
    result = result.with_columns(pl.col("contig").is_in(circular_contigs).alias("circular_contig"))
    return result


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

    data = data.with_columns([
        pl.col("KEGG_V-score").cast(pl.Float64, strict=False),
        pl.col("Pfam_V-score").cast(pl.Float64, strict=False),
        pl.col("PHROG_V-score").cast(pl.Float64, strict=False),
        pl.col("window_avg_KEGG_VL-score").cast(pl.Float64, strict=False),
        pl.col("window_avg_Pfam_VL-score").cast(pl.Float64, strict=False),
        pl.col("window_avg_PHROG_VL-score").cast(pl.Float64, strict=False),
    ])

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
            return df.with_columns([
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
            ])

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
            bad_gene = np.zeros(n, dtype=np.bool_)
            v_anchor = hallmark_any
            all3_ge_min = hallmark_any
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

        return df.with_columns([
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
        ])

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

    kegg_v = df["KEGG_V-score"].to_numpy()
    pfam_v = df["Pfam_V-score"].to_numpy()
    phrog_v = df["PHROG_V-score"].to_numpy()

    df = df.with_columns([
        pl.Series("window_avg_KEGG_VL-score", window_avg(kegg_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_Pfam_VL-score", window_avg(pfam_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_PHROG_VL-score", window_avg(phrog_vl, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_KEGG_V-score", window_avg(kegg_v, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_Pfam_V-score", window_avg(pfam_v, lengths, window_size, minimum_percentage)),
        pl.Series("window_avg_PHROG_V-score", window_avg(phrog_v, lengths, window_size, minimum_percentage)),
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

    missing = [c for c in feature_names if c not in df.columns]
    if missing:
        df = df.with_columns([pl.lit(np.nan, dtype=pl.Float32).alias(c) for c in missing])

    feat_df = df.select([pl.col(c).cast(pl.Float32, strict=False) for c in feature_names])
    X_all = feat_df.to_numpy()
    n = X_all.shape[0]

    y_proba_all = np.empty(n, dtype=np.float32)

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

    left_dist = np.full(n, np.iinfo(np.int32).min, dtype=np.float64)
    right_dist = np.full(n, np.iinfo(np.int32).min, dtype=np.float64)

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
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_scores = contig_data['KEGG_V-score'].to_numpy()
    pfam_scores = contig_data['Pfam_V-score'].to_numpy()
    phrog_scores = contig_data['PHROG_V-score'].to_numpy()

    kegg_left, kegg_right = flank_distance_vscores(lengths, kegg_scores, minimum_vscore)
    pfam_left, pfam_right = flank_distance_vscores(lengths, pfam_scores, minimum_vscore)
    phrog_left, phrog_right = flank_distance_vscores(lengths, phrog_scores, minimum_vscore)

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
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])


@njit
def flank_distance_in_set(lengths, in_set):
    n = len(lengths)
    positions = np.zeros(n, dtype=np.float64)
    for i in range(1, n):
        positions[i] = positions[i-1] + lengths[i-1]

    left_dist = np.full(n, np.iinfo(np.int32).min, dtype=np.float64)
    right_dist = np.full(n, np.iinfo(np.int32).min, dtype=np.float64)

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
    arr = np.zeros(len(hmm_ids), dtype=np.int64)
    for i, val in enumerate(hmm_ids):
        if val is not None and val in valid_set:
            arr[i] = 1
    return arr


def verify_flanking_hallmark(contig_data, hallmark_accessions):
    """
    For each gene, calculate nucleotide distance to nearest left/right hallmark gene.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()

    kegg_arr = create_in_set_array(kegg_hmm, hallmark_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, hallmark_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, hallmark_accessions)

    kegg_left, kegg_right = flank_distance_in_set(lengths, kegg_arr)
    pfam_left, pfam_right = flank_distance_in_set(lengths, pfam_arr)
    phrog_left, phrog_right = flank_distance_in_set(lengths, phrog_arr)

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
    })
    return contig_data.join(flanks_df, on=["contig", "gene_number"])


@njit
def flank_nearest_mge_scores(scores, distances, mge_mask):
    """
    For each gene, report the score of the nearest left/right MGE gene.
    Only written when a valid MGE distance exists (distances[i] != INT_PLACEHOLDER sentinel).
    Scores remain NaN when no MGE neighbor exists on that side.
    """
    n = len(scores)
    left = np.full(n, np.nan, dtype=np.float32)
    right = np.full(n, np.nan, dtype=np.float32)

    # Only carry a score forward when the corresponding distance is valid,
    # i.e. there is an actual MGE to the left. Track whether last_idx was set
    # to know if left[i] should be populated at all.
    last_idx = -1
    for i in range(n):
        if last_idx != -1:
            left[i] = np.float32(scores[last_idx])
        if mge_mask[i] and not np.isnan(scores[i]):
            last_idx = i

    next_idx = -1
    for i in range(n - 1, -1, -1):
        if next_idx != -1:
            right[i] = np.float32(scores[next_idx])
        if mge_mask[i] and not np.isnan(scores[i]):
            next_idx = i

    return left, right


def report_flanking_mge_vscores(contig_data, mobile_accessions):
    """
    For each gene, report the KEGG/Pfam/PHROG V/VL-score of the nearest left/right MGE gene.
    Score is null when no MGE exists on that side within the contig.
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    add_cols = {}
    for db in ["KEGG", "Pfam", "PHROG"]:
        hmm_col = f"{db}_hmm_id"
        mge_mask = np.asarray(
            [x in mobile_accessions if x is not None else False for x in contig_data[hmm_col].to_list()],
            dtype=np.bool_
        )
        # Compute distances once per db to reuse as the validity mask for scores
        dist_left, dist_right = flank_distance_in_set(lengths, mge_mask.astype(np.int64))

        for sfx in ["V-score", "VL-score"]:
            scores = contig_data[f"{db}_{sfx}"].to_numpy().astype(np.float32)
            left, right = flank_nearest_mge_scores(scores, dist_left, mge_mask)

            # null out any score where the corresponding distance is still the
            # INT_PLACEHOLDER sentinel (meaning no MGE neighbor on that side).
            placeholder = float(np.iinfo(np.int32).min)
            left = np.where(dist_left == placeholder, np.nan, left).astype(np.float32)
            right = np.where(dist_right == placeholder, np.nan, right).astype(np.float32)

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
    """
    contig_data = contig_data.sort(["contig", "gene_number"])
    lengths = contig_data['gene_length_bases'].to_numpy()
    kegg_hmm = contig_data['KEGG_hmm_id'].to_list()
    pfam_hmm = contig_data['Pfam_hmm_id'].to_list()
    phrog_hmm = contig_data['PHROG_hmm_id'].to_list()

    kegg_arr = create_in_set_array(kegg_hmm, mobile_accessions)
    pfam_arr = create_in_set_array(pfam_hmm, mobile_accessions)
    phrog_arr = create_in_set_array(phrog_hmm, mobile_accessions)

    kegg_left, kegg_right = flank_distance_in_set(lengths, kegg_arr)
    pfam_left, pfam_right = flank_distance_in_set(lengths, pfam_arr)
    phrog_left, phrog_right = flank_distance_in_set(lengths, phrog_arr)

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
    })
    # Join back so the result is keyed on (contig, gene_number), preserving order
    return contig_data.join(flanks_df, on=["contig", "gene_number"])


def _replace_placeholder_with_null(data: pl.DataFrame, placeholder: int = INT_PLACEHOLDER) -> pl.DataFrame:
    """
    Convert INT_PLACEHOLDER sentinel values to null in all *_dist columns,
    then cast to Float32. Operates on all columns ending in '_dist'.
    """
    dist_cols = [col for col in data.columns if col.endswith('_dist')]
    logger.debug(f"Recasting {len(dist_cols)} distance columns to Float32 with null for placeholder.")
    for col in dist_cols:
        # Guard: if a column was read back from CSV as Utf8 (all-null edge case),
        # cast directly to Float32 rather than attempting a numeric comparison.
        if data[col].dtype == pl.Utf8:
            data = data.with_columns(
                pl.col(col).cast(pl.Float32, strict=False).alias(col)
            )
        else:
            data = data.with_columns(
                pl.when(pl.col(col) == placeholder)
                .then(pl.lit(None, dtype=pl.Float32))
                .otherwise(pl.col(col).cast(pl.Float32))
                .alias(col)
            )
    return data


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
                    n_cpus=1):

    logger.debug(f"Calculating lengths for {data.shape[0]:,} genes.")
    data = calculate_gene_lengths(data)
    # Snapshot taken after gene-length calculation: V/VL scores unchanged from
    # input, gene_length_bases present. Used for MGE distance/score functions
    # which only need input-level columns.
    data_orig = data

    logger.info("Calculating distance to contig ends.")
    data = calculate_contig_end_distances(data)

    logger.info("Calculating contig statistics.")
    data = calculate_contig_statistics(data, circular_contigs)

    logger.info("Calculating window statistics.")
    score_columns = [
        "KEGG_V-score", "KEGG_VL-score", "Pfam_V-score", "Pfam_VL-score",
        "PHROG_V-score", "PHROG_VL-score",
        "contig_avg_KEGG_V-score", "contig_avg_Pfam_V-score", "contig_avg_PHROG_V-score",
        "contig_avg_KEGG_VL-score", "contig_avg_Pfam_VL-score", "contig_avg_PHROG_VL-score",
    ]
    for col in score_columns:
        if col in data.columns:
            data = data.with_columns(pl.col(col).cast(pl.Float64, strict=False))

    data = calculate_window_statistics(data, window_size, minimum_percentage, n_cpus)
    data = data.unique()

    logger.info("Identifying ambiguous viral/nonviral contig regions.")
    data = identify_viral_nonviral_regions(
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

    # Flanking viral/hallmark distances (joined by key)
    if use_hallmark and hallmark_accessions is not None:
        logger.info("Calculating distance of nearest flanking hallmark genes.")
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
        contig_dfs = data.partition_by("contig")
        with ThreadPoolExecutor(max_workers=n_cpus) as executor:
            futures = [
                executor.submit(verify_flanking_vscores, df, minimum_vscore)
                for df in contig_dfs
            ]
            results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc=f"Checking flanks for V-score={minimum_vscore}", unit="contig")]
        data = pl.concat(results, how="vertical")

    # Partition data_orig (stable, input-level columns) for MGE work
    logger.info("Calculating distance of nearest mobile genetic element genes.")
    contig_dfs_orig = data_orig.partition_by("contig")

    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [
            executor.submit(check_flanking_mge, df, mobile_accessions)
            for df in contig_dfs_orig
        ]
        mge_dist_results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanks for mobile genes", unit="contig")]

    mge_dist_df = pl.concat(mge_dist_results, how="vertical")
    mge_dist_cols = [c for c in mge_dist_df.columns if c.endswith("_dist") and "MGE" in c]

    if not mge_dist_cols:
        raise ValueError("No flanking MGE distance columns found.")

    # Join on 'protein' (keyed) instead of positional hstack
    mge_dist_df = mge_dist_df.select(["protein", *mge_dist_cols])
    data = data.join(mge_dist_df, on="protein", how="left")

    # Replace INT_PLACEHOLDER with null across ALL _dist columns
    # Covers contig-end, viral/hallmark, and MGE distance columns in one pass.
    data = _replace_placeholder_with_null(data, placeholder=INT_PLACEHOLDER)

    # Reuse same data_orig partitions so MGE V/VL-score join keys match exactly.
    # report_flanking_mge_vscores only reads V/VL scores + hmm_ids, both present
    # in data_orig.
    logger.info("Calculating flanking MGE V/VL-scores.")
    with ThreadPoolExecutor(max_workers=n_cpus) as executor:
        futures = [executor.submit(report_flanking_mge_vscores, df, mobile_accessions) for df in contig_dfs_orig]
        mge_score_results = [f.result() for f in tqdm(as_completed(futures), total=len(futures), desc="Checking flanking MGE V/VL-scores", unit="contig")]

    mge_score_df = pl.concat(mge_score_results, how="vertical")
    score_cols_to_add = [c for c in mge_score_df.columns if c not in data.columns or c in ["contig", "gene_number"]]
    data = data.join(mge_score_df.select(score_cols_to_add), on=["contig", "gene_number"], how="left")

    # Null-mask MGE score columns where the corresponding distance is null.
    # flank_nearest_mge_scores can carry a score from beyond the contig boundary
    # after INT_PLACEHOLDER replacement; this ensures score == null iff dist == null.
    for db in ["KEGG", "Pfam", "PHROG"]:
        for side in ["left", "right"]:
            dist_col = f"{db}_MGE_{side}_dist"
            for sfx in ["V-score", "VL-score"]:
                score_col = f"{db}_{sfx}_{side}_MGE"
                if dist_col in data.columns and score_col in data.columns:
                    data = data.with_columns(
                        pl.when(pl.col(dist_col).is_null())
                        .then(pl.lit(None, dtype=pl.Float32))
                        .otherwise(pl.col(score_col))
                        .alias(score_col)
                    )

    logger.info("Assigning viral origin confidence using LightGBM with genome context features.")
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

    logger.info("Starting genome context analysis...")
    os.makedirs(outparent, exist_ok=True)

    input_file = as_parquet_path_if_enabled(snakemake.params.gene_index_annotated, save_to_parquet)
    output_file = as_parquet_path_if_enabled(snakemake.params.context_table, save_to_parquet)
    circular_contigs_file = as_parquet_path_if_enabled(snakemake.params.circular_contigs, save_to_parquet)

    if not os.path.exists(circular_contigs_file) or os.path.getsize(circular_contigs_file) == 0:
        circular_contigs = set()
        logger.warning("No results found for checking circular contigs. All values for 'circular_contig' will be False.")
        if save_to_parquet:
            data = pl.read_parquet(input_file)
        else:
            data = pl.read_csv(input_file, separator='\t')
    else:
        if save_to_parquet:
            circular_contigs = set(pl.read_parquet(circular_contigs_file)['contig'].to_list())
            data = pl.read_parquet(input_file)
        else:
            circular_contigs = set(pl.read_csv(circular_contigs_file, separator='\t')['contig'].to_list())
            data = pl.read_csv(input_file, separator='\t')

    logger.debug(f"Loaded data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")
    data = data.sort(["contig", "gene_number"]).unique()
    logger.debug(f"Unique data with {data.shape[0]:,} rows and {data.shape[1]:,} columns.")

    use_hallmark = snakemake.params.use_hallmark
    hallmark_path = snakemake.params.hallmark_path
    hallmark_ids = None
    if use_hallmark:
        hallmark_data = pl.read_csv(hallmark_path, separator='\t')
        hallmark_ids = set(hallmark_data['id'])

    mobile_genes_path = snakemake.params.mobile_genes_path
    mobile_ids = None
    if mobile_genes_path:
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
        n_cpus
    )

    if save_to_parquet:
        processed_data.write_parquet(output_file)
    else:
        processed_data.write_csv(output_file, separator='\t', include_header=True)
    logger.info("Genome context analysis completed.")


if __name__ == "__main__":
    main()