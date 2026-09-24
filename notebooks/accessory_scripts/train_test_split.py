#!/usr/bin/env python3

from __future__ import annotations

import argparse
import hashlib
import logging
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import polars as pl
from concurrent.futures import ProcessPoolExecutor, as_completed


def _stable_hash_u64(obj: Union[str, int, float, tuple, bytes]) -> int:
    if not isinstance(obj, (bytes, bytearray)):
        obj = repr(obj).encode("utf-8")
    return int.from_bytes(hashlib.sha256(obj).digest()[:8], "little", signed=False)


def _seed_sequence(master_seed: int, tag: Union[str, int, tuple]) -> np.random.SeedSequence:
    mix = int(master_seed) ^ (_stable_hash_u64(tag) & 0x7FFF_FFFF_FFFF_FFFF)
    return np.random.SeedSequence(mix)


def read_cluster_file_to_contig_block_map(cluster_path: str, sep: str = "\t") -> dict[str, int]:
    contig_block_map: dict[str, int] = {}
    with open(cluster_path, "r") as fh:
        for cid, raw in enumerate(fh):
            line = raw.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(sep) if p.strip()]
            if not parts:
                continue
            is_pair_format = (len(parts) % 2 == 0 and parts[0] in {"Virus", "Host", "MGE"})
            if is_pair_format:
                it = iter(parts)
                for _, ctg in zip(it, it):
                    if ctg not in contig_block_map:
                        contig_block_map[ctg] = cid
            else:
                for ctg in parts:
                    if ctg not in contig_block_map:
                        contig_block_map[ctg] = cid
    return contig_block_map


def attach_block_id_fast(
    df: pl.LazyFrame,
    contig_block_map: dict[str, int],
    contig_col: str = "Contig",
) -> pl.LazyFrame:
    # Sort by key so the DataFrame row order is fully deterministic regardless
    # of dict insertion order.
    sorted_pairs = sorted(contig_block_map.items())
    known_ctgs = [k for k, _ in sorted_pairs]
    known_bids = [int(v) for _, v in sorted_pairs]

    if known_ctgs:
        base = pl.DataFrame(
            {contig_col: known_ctgs, "_block_id_known": known_bids},
            schema={contig_col: pl.Utf8, "_block_id_known": pl.UInt64},
        ).lazy()
        dfj = df.join(base, on=contig_col, how="left")
    else:
        dfj = df.with_columns(pl.lit(None, dtype=pl.UInt64).alias("_block_id_known"))

    topbit = 1 << 63
    mask63 = (1 << 63) - 1

    h = pl.col(contig_col).hash(seed=0).cast(pl.UInt64)
    unseen_expr = (h & pl.lit(mask63).cast(pl.UInt64)) | pl.lit(topbit).cast(pl.UInt64)

    return (
        dfj.with_columns(
            pl.when(pl.col("_block_id_known").is_not_null())
            .then(pl.col("_block_id_known"))
            .otherwise(unseen_expr)
            .alias("_block_id")
        )
        .drop("_block_id_known")
    )


def contig_train_test_split_block_atomic(
    df: pl.LazyFrame,
    contig_block_map: dict[str, int],
    train_frac: float,
    seed: int,
    contig_col: str = "Contig",
) -> Tuple[pl.DataFrame, pl.DataFrame]:
    ss = _seed_sequence(seed, ("contig_train_test_split_block_atomic", int(seed)))
    rng = np.random.default_rng(ss)

    df_wb = attach_block_id_fast(df, contig_block_map=contig_block_map, contig_col=contig_col)

    block_sizes = (
        df_wb.group_by("_block_id")
        .agg(pl.len().alias("n"))
        .sort("_block_id")
        .collect()
    )

    bids = block_sizes["_block_id"].to_numpy()
    sizes = block_sizes["n"].to_numpy()
    total = int(sizes.sum())
    target = float(train_frac) * float(total)

    order = np.arange(bids.size, dtype=np.int64)
    rng.shuffle(order)

    selected = np.zeros(bids.size, dtype=bool)
    cum = 0.0
    for i in order:
        s = float(sizes[i])
        if cum + s <= target:
            selected[i] = True
            cum += s

    remaining = np.where(~selected)[0]
    if remaining.size:
        idx = remaining[np.argmin(np.abs((cum + sizes[remaining]) - target))]
        if abs(cum + sizes[idx] - target) < abs(cum - target):
            selected[idx] = True

    train_bids = bids[selected]
    train_bid_list = train_bids.tolist()

    df_e = df_wb.collect().sort([contig_col, "contig_pos_start", "contig_pos_end"])
    train_df = df_e.filter(pl.col("_block_id").is_in(train_bid_list)).drop("_block_id")
    test_df = df_e.filter(~pl.col("_block_id").is_in(train_bid_list)).drop("_block_id")
    return train_df, test_df


def _select_contigs_for_target(contig_ids: np.ndarray, n_rows: np.ndarray, target_count: int, seed: int) -> np.ndarray:
    if target_count <= 0 or contig_ids.size == 0:
        return np.array([], dtype=object)

    ss = _seed_sequence(seed, ("select_contigs_for_target", int(seed)))
    rng = np.random.default_rng(ss)

    perm = np.arange(contig_ids.size, dtype=np.int64)
    rng.shuffle(perm)

    cids = contig_ids[perm]
    rows = n_rows[perm]

    cumsum = np.cumsum(rows)
    k = int(np.searchsorted(cumsum, target_count, side="right"))
    k = min(k, cids.size - 1)

    cand = [k]
    if k > 0:
        cand.append(k - 1)
    k = min(cand, key=lambda i: abs(int(cumsum[i]) - int(target_count)))
    return cids[: k + 1]


def enforce_balanced_training_set(
    train_df: pl.DataFrame,
    seed: int,
    contig_col: str = "Contig",
    label_col: str = "True Positive",
) -> pl.DataFrame:
    if train_df.is_empty():
        return train_df

    for c in (contig_col, label_col):
        if c not in train_df.columns:
            raise ValueError(f"Missing required column: {c}")

    contig_labels = (
        train_df.group_by(contig_col)
        .agg(
            pl.col(label_col).any().alias("has_pos"),
            pl.len().alias("n_rows"),
        )
        .with_columns(pl.when(pl.col("has_pos")).then(pl.lit("pos")).otherwise(pl.lit("neg")).alias("contig_class"))
        .select([contig_col, "contig_class", "n_rows"])
        .sort(contig_col)
    )

    totals = contig_labels.group_by("contig_class").agg(pl.sum("n_rows").alias("n")).sort("contig_class")
    tot_map = {row[0]: int(row[1]) for row in totals.iter_rows()}
    n_pos = tot_map.get("pos", 0)
    n_neg = tot_map.get("neg", 0)
    if n_pos == 0 or n_neg == 0 or n_pos == n_neg:
        return train_df.sort(["Source", contig_col, "contig_pos_start", "contig_pos_end"])

    if n_pos < n_neg:
        majority = "neg"
        target = n_pos
    else:
        majority = "pos"
        target = n_neg

    maj = contig_labels.filter(pl.col("contig_class") == majority).select([contig_col, "n_rows"])
    maj_ids = np.array(maj[contig_col].to_list(), dtype=object)
    maj_rows = np.array(maj["n_rows"].to_list(), dtype=np.int64)

    rs = int(seed) + int(_stable_hash_u64(("train_balance_contig", majority)) % 1_000_003)
    keep_contigs = _select_contigs_for_target(maj_ids, maj_rows, int(target), int(rs))
    # Use a sorted list, not a set, so the sequence passed to is_in() is
    # deterministic across runs (set iteration order varies with PYTHONHASHSEED).
    keep_list = sorted(str(c) for c in keep_contigs.tolist())

    df_majority = train_df.filter(pl.col(contig_col).cast(pl.Utf8).is_in(keep_list))

    minority_contigs = (
        contig_labels.filter(pl.col("contig_class") != majority)
        .select(contig_col)
        .to_series()
        .to_list()
    )
    df_minority = train_df.filter(pl.col(contig_col).is_in(minority_contigs))

    out = pl.concat([df_minority, df_majority]).sort(["Source", contig_col, "contig_pos_start", "contig_pos_end"])
    return out


def maximize_test_pool_with_source_bounds(
    test_pool: pl.DataFrame,
    contig_block_map: dict[str, int],
    seed: int,
    source: str,
    min_frac: float,
    max_frac: float,
) -> pl.DataFrame:
    if test_pool.is_empty():
        return test_pool

    if source not in {"Virus", "Host", "MGE"}:
        raise ValueError(f"Invalid source: {source}")

    df = attach_block_id_fast(test_pool.lazy(), contig_block_map=contig_block_map, contig_col="Contig").collect()

    block_stats = (
        df.group_by(["_block_id", "Source"])
        .agg(pl.len().alias("n"))
        .pivot(values="n", index="_block_id", on="Source")
        .fill_null(0)
    )

    for col in ("Virus", "Host", "MGE"):
        if col not in block_stats.columns:
            block_stats = block_stats.with_columns(pl.lit(0).cast(pl.Int64).alias(col))
        else:
            block_stats = block_stats.with_columns(pl.col(col).cast(pl.Int64))

    total_expr = pl.col("Virus") + pl.col("Host") + pl.col("MGE")
    src_expr = pl.col(source)

    block_stats = (
        block_stats
        .with_columns(total_expr.alias("n_total"))
        .with_columns(
            pl.when(total_expr > 0)
            .then(src_expr / total_expr)
            .otherwise(0.0)
            .alias("sfrac")
        )
        .select(["_block_id", "Virus", "Host", "MGE", "n_total", "sfrac"])
        .sort("_block_id")
    )

    tot_s = int(block_stats[source].sum())
    tot_t = int(block_stats["n_total"].sum())
    if tot_t <= 0:
        return df.drop("_block_id").sort(["Source", "Contig", "contig_pos_start", "contig_pos_end"])

    frac0 = tot_s / tot_t
    if min_frac <= frac0 <= max_frac:
        return df.drop("_block_id").sort(["Source", "Contig", "contig_pos_start", "contig_pos_end"])

    b_ids = block_stats["_block_id"].to_list()
    tie = [int(_stable_hash_u64(("pool_bound", source, int(seed), int(b)))) for b in b_ids]
    block_stats = block_stats.with_columns(pl.Series("_tie", tie, dtype=pl.UInt64))

    if frac0 < min_frac:
        ordered = block_stats.sort(["sfrac", "n_total", "_tie"], descending=[False, False, False])
        target = float(min_frac)
    else:
        ordered = block_stats.sort(["sfrac", "n_total", "_tie"], descending=[True, False, False])
        target = float(max_frac)

    remove_blocks: List[int] = []
    cur_s = tot_s
    cur_t = tot_t

    best_remove_blocks: List[int] = []
    best_dist = abs(frac0 - target)
    best_t = cur_t

    for row in ordered.iter_rows(named=True):
        if cur_t <= 0:
            break
        cur_frac = cur_s / cur_t
        if min_frac <= cur_frac <= max_frac:
            break

        b = int(row["_block_id"])
        bs = int(row[source])
        bt = int(row["n_total"])
        if bt <= 0:
            continue

        new_s = cur_s - bs
        new_t = cur_t - bt
        if new_t <= 0:
            continue
        new_frac = new_s / new_t

        remove_blocks.append(b)
        cur_s = new_s
        cur_t = new_t

        dist = abs(new_frac - target)
        improved = (dist < best_dist) or (dist == best_dist and new_t > best_t)
        if improved:
            best_dist = dist
            best_t = new_t
            best_remove_blocks = list(remove_blocks)

    keep_df = df
    if best_remove_blocks:
        keep_df = keep_df.filter(~pl.col("_block_id").is_in(best_remove_blocks))

    return keep_df.drop("_block_id").sort(["Source", "Contig", "contig_pos_start", "contig_pos_end"])


def build_test_provirus(
    test_pool_raw: pl.DataFrame,
    contig_type_col: str = "region_contig_type",
    contig_col: str = "Contig",
    source_col: str = "Source",
) -> pl.DataFrame:
    if test_pool_raw.is_empty():
        return test_pool_raw

    for c in (contig_type_col, contig_col, source_col):
        if c not in test_pool_raw.columns:
            raise ValueError(f"Missing required column for test_provirus: {c}")

    df_mixed = test_pool_raw.filter(pl.col(contig_type_col) == "chromosome_mixed")
    if df_mixed.is_empty():
        return df_mixed

    keep = (
        df_mixed.group_by(contig_col)
        .agg(
            (pl.col(source_col) == "Host").any().alias("has_host"),
            (pl.col(source_col) == "Virus").any().alias("has_virus"),
        )
        .sort(contig_col)  # group_by output is unordered; sort before filtering/joining
        .filter(pl.col("has_host") & pl.col("has_virus"))
        .select(contig_col)
    )

    out = (
        df_mixed.join(keep, on=contig_col, how="inner")
        .sort([source_col, contig_col, "contig_pos_start", "contig_pos_end"])
    )
    return out


def _pool_init(chosen: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = chosen.size
    in_list = np.flatnonzero(chosen).astype(np.int64, copy=False)
    out_list = np.flatnonzero(~chosen).astype(np.int64, copy=False)
    in_pos = np.full(n, -1, dtype=np.int64)
    out_pos = np.full(n, -1, dtype=np.int64)
    in_pos[in_list] = np.arange(in_list.size, dtype=np.int64)
    out_pos[out_list] = np.arange(out_list.size, dtype=np.int64)
    return in_list, out_list, in_pos, out_pos


def _pool_remove(idx: int, lst: np.ndarray, pos: np.ndarray) -> np.ndarray:
    p = pos[idx]
    if p < 0:
        return lst
    last = lst[-1]
    lst[p] = last
    pos[last] = p
    pos[idx] = -1
    return lst[:-1]


def _pool_add(idx: int, lst: np.ndarray, pos: np.ndarray) -> np.ndarray:
    pos[idx] = lst.size
    return np.append(lst, idx)


def solve_contig_subset_for_distribution_fast(
    contrib: np.ndarray,
    totals: np.ndarray,
    req: np.ndarray,
    desired_total: int,
    seed: int,
    min_total: int,
    max_total: int,
    n_restarts: int,
    max_iters: int,
    stagnation_limit: int,
    move_probs: Tuple[float, float, float],
) -> Tuple[Optional[np.ndarray], float, int]:
    n = int(totals.size)
    if n == 0:
        return None, float("inf"), 0

    totals_i64 = totals.astype(np.int64, copy=False)
    contrib_i64 = contrib.astype(np.int64, copy=False)
    req_f64 = req.astype(np.float64, copy=False)

    frac = np.zeros((n, 3), dtype=np.float64)
    nz = totals_i64 > 0
    frac[nz, :] = contrib_i64[nz, :].astype(np.float64) / totals_i64[nz].astype(np.float64)[:, None]
    dist0 = np.max(np.abs(frac - req_f64[None, :]), axis=1)
    order_by_dist = np.argsort(dist0, kind="mergesort")

    def _maxerr(cur: np.ndarray, tot: int) -> float:
        if tot <= 0:
            return float("inf")
        return float(np.max(np.abs((cur.astype(np.float64) / float(tot)) - req_f64)))

    def _score(cur: np.ndarray, tot: int) -> Tuple[float, float, int]:
        maxerr = _maxerr(cur, tot)
        size_pen = (float(desired_total) - float(tot)) / float(max(1, desired_total))
        size_pen = max(0.0, size_pen)
        return (maxerr, size_pen, -int(tot))

    p1, p2, p3 = move_probs
    ps = p1 + p2 + p3
    if ps <= 0:
        p1, p2, p3 = 0.70, 0.15, 0.15
        ps = 1.0
    p1, p2, p3 = p1 / ps, p2 / ps, p3 / ps

    ss = _seed_sequence(seed, ("restart_stream", int(seed)))
    restart_seeds = ss.spawn(int(n_restarts))

    best_chosen: Optional[np.ndarray] = None
    best_tuple: Tuple[float, float, int] = (float("inf"), float("inf"), 0)

    for r in range(int(n_restarts)):
        rng = np.random.default_rng(restart_seeds[r])

        chosen = np.zeros(n, dtype=bool)
        cur = np.zeros(3, dtype=np.int64)
        cur_total = 0

        perm = np.arange(n, dtype=np.int64)
        rng.shuffle(perm)

        for i in perm:
            t = int(totals_i64[i])
            if t <= 0:
                continue
            if cur_total + t > max_total:
                continue
            chosen[i] = True
            cur += contrib_i64[i]
            cur_total += t
            if cur_total >= min_total:
                break

        if cur_total < min_total:
            continue

        for i in order_by_dist:
            if chosen[i]:
                continue
            t = int(totals_i64[i])
            if t <= 0:
                continue
            if cur_total + t > max_total:
                continue
            if cur_total >= desired_total:
                break
            chosen[i] = True
            cur += contrib_i64[i]
            cur_total += t

        cur_tuple = _score(cur, cur_total)

        in_list, out_list, in_pos, out_pos = _pool_init(chosen)
        if in_list.size == 0 or out_list.size == 0:
            continue

        no_improve = 0
        for _ in range(int(max_iters)):
            if no_improve >= int(stagnation_limit):
                break
            if in_list.size == 0 or out_list.size == 0:
                break

            u = float(rng.random())

            if u < p1:
                i = int(in_list[rng.integers(0, in_list.size)])
                j = int(out_list[rng.integers(0, out_list.size)])
                t2 = cur_total - int(totals_i64[i]) + int(totals_i64[j])
                if t2 < min_total or t2 > max_total:
                    no_improve += 1
                    continue
                c2 = cur - contrib_i64[i] + contrib_i64[j]
                tup2 = _score(c2, t2)
                if tup2 < cur_tuple:
                    chosen[i] = False
                    chosen[j] = True

                    in_list = _pool_remove(i, in_list, in_pos)
                    out_list = _pool_add(i, out_list, out_pos)

                    out_list = _pool_remove(j, out_list, out_pos)
                    in_list = _pool_add(j, in_list, in_pos)

                    cur = c2
                    cur_total = t2
                    cur_tuple = tup2
                    no_improve = 0
                else:
                    no_improve += 1

            elif u < p1 + p2:
                if out_list.size < 2:
                    no_improve += 1
                    continue
                i = int(in_list[rng.integers(0, in_list.size)])
                j1 = int(out_list[rng.integers(0, out_list.size)])
                j2 = int(out_list[rng.integers(0, out_list.size)])
                if j1 == j2:
                    no_improve += 1
                    continue
                t2 = cur_total - int(totals_i64[i]) + int(totals_i64[j1]) + int(totals_i64[j2])
                if t2 < min_total or t2 > max_total:
                    no_improve += 1
                    continue
                c2 = cur - contrib_i64[i] + contrib_i64[j1] + contrib_i64[j2]
                tup2 = _score(c2, t2)
                if tup2 < cur_tuple:
                    chosen[i] = False
                    chosen[j1] = True
                    chosen[j2] = True

                    in_list = _pool_remove(i, in_list, in_pos)
                    out_list = _pool_add(i, out_list, out_pos)

                    out_list = _pool_remove(j1, out_list, out_pos)
                    in_list = _pool_add(j1, in_list, in_pos)

                    out_list = _pool_remove(j2, out_list, out_pos)
                    in_list = _pool_add(j2, in_list, in_pos)

                    cur = c2
                    cur_total = t2
                    cur_tuple = tup2
                    no_improve = 0
                else:
                    no_improve += 1

            else:
                if in_list.size < 2:
                    no_improve += 1
                    continue
                i1 = int(in_list[rng.integers(0, in_list.size)])
                i2 = int(in_list[rng.integers(0, in_list.size)])
                if i1 == i2:
                    no_improve += 1
                    continue
                j = int(out_list[rng.integers(0, out_list.size)])
                t2 = cur_total - int(totals_i64[i1]) - int(totals_i64[i2]) + int(totals_i64[j])
                if t2 < min_total or t2 > max_total:
                    no_improve += 1
                    continue
                c2 = cur - contrib_i64[i1] - contrib_i64[i2] + contrib_i64[j]
                tup2 = _score(c2, t2)
                if tup2 < cur_tuple:
                    chosen[i1] = False
                    chosen[i2] = False
                    chosen[j] = True

                    in_list = _pool_remove(i1, in_list, in_pos)
                    out_list = _pool_add(i1, out_list, out_pos)

                    in_list = _pool_remove(i2, in_list, in_pos)
                    out_list = _pool_add(i2, out_list, out_pos)

                    out_list = _pool_remove(j, out_list, out_pos)
                    in_list = _pool_add(j, in_list, in_pos)

                    cur = c2
                    cur_total = t2
                    cur_tuple = tup2
                    no_improve = 0
                else:
                    no_improve += 1

        if cur_tuple < best_tuple:
            best_tuple = cur_tuple
            best_chosen = chosen.copy()

        if best_tuple[0] <= 0.005:
            break

    if best_chosen is None:
        return None, float("inf"), 0
    return best_chosen, float(best_tuple[0]), int(-best_tuple[2])


@dataclass
class MemMapBundle:
    contrib_path: str
    totals_path: str
    n: int

    def load(self) -> tuple[np.ndarray, np.ndarray]:
        contrib = np.memmap(self.contrib_path, dtype=np.int64, mode="r", shape=(self.n, 3))
        totals = np.memmap(self.totals_path, dtype=np.int64, mode="r", shape=(self.n,))
        return contrib, totals


def _write_memmaps(contrib: np.ndarray, totals: np.ndarray, tmpdir: Path, tag: str) -> MemMapBundle:
    n = int(totals.size)
    contrib_path = str(tmpdir / f"contrib_i64_{tag}.mmap")
    totals_path = str(tmpdir / f"totals_i64_{tag}.mmap")

    cm = np.memmap(contrib_path, dtype=np.int64, mode="w+", shape=(n, 3))
    tm = np.memmap(totals_path, dtype=np.int64, mode="w+", shape=(n,))

    cm[:] = contrib.astype(np.int64, copy=False)
    tm[:] = totals.astype(np.int64, copy=False)

    cm.flush()
    tm.flush()

    return MemMapBundle(contrib_path=contrib_path, totals_path=totals_path, n=n)


def _dataset_worker(task: dict) -> dict:
    name = str(task["name"])
    seed = int(task["seed"])
    bundle: MemMapBundle = task["bundle"]
    req = np.array(task["req"], dtype=np.float64)

    target_size = int(task["target_size"])
    allow_downsize = bool(task["allow_downsize"])
    downsize_grid = list(task["downsize_grid"])
    total_tol_frac = float(task["total_tol_frac"])
    max_err_threshold = float(task["max_err_threshold"])
    n_restarts = int(task["n_restarts"])
    max_iters = int(task["max_iters"])
    stagnation_limit = int(task["stagnation_limit"])
    move_probs = tuple(task["move_probs"])
    min_target_frac = float(task["min_target_frac"])
    min_pool_frac = float(task["min_pool_frac"])
    min_abs_keep = int(task["min_abs_keep"])

    contrib, totals = bundle.load()

    max_possible = int(totals.sum())
    if max_possible <= 0:
        return {"name": name, "ok": False, "err": float("inf"), "chosen_idx": []}

    base_target = int(min(int(target_size), int(max_possible)))
    if base_target <= 0:
        return {"name": name, "ok": False, "err": float("inf"), "chosen_idx": []}

    if float(req.sum()) <= 0:
        return {"name": name, "ok": False, "err": float("inf"), "chosen_idx": []}
    req = req / float(req.sum())

    pool_floor = int(np.ceil(float(min_pool_frac) * float(max_possible)))
    target_floor = int(np.ceil(float(min_target_frac) * float(base_target)))
    hard_floor = int(max(1, min(max_possible, max(pool_floor, target_floor, min_abs_keep))))

    min_scale = float(hard_floor) / float(base_target) if base_target > 0 else 1.0
    min_scale = min(1.0, max(0.0, min_scale))

    if allow_downsize:
        grid = [float(x) for x in downsize_grid if float(x) >= min_scale]
        grid.append(min_scale)
        grid = sorted(set(grid), reverse=True)
    else:
        grid = [1.0]

    ss_grid = _seed_sequence(seed, ("downsize_grid", name))
    grid_seeds = ss_grid.spawn(len(grid))

    best_chosen = None
    best_err = float("inf")

    for gi, scale in enumerate(grid):
        desired = int(max(hard_floor, round(float(base_target) * float(scale))))
        desired = int(min(max_possible, max(1, desired)))

        min_total = int(max(hard_floor, round(float(desired) * (1.0 - float(total_tol_frac)))))
        max_total = int(min(max_possible, round(float(desired) * (1.0 + float(total_tol_frac)))))
        if min_total > max_total:
            min_total = max(1, min(max_total, min_total))

        gs = int(np.random.default_rng(grid_seeds[gi]).integers(0, 2**63 - 1))

        chosen, err, _ = solve_contig_subset_for_distribution_fast(
            contrib=contrib,
            totals=totals,
            req=req,
            desired_total=int(desired),
            seed=int(gs),
            min_total=int(min_total),
            max_total=int(max_total),
            n_restarts=int(n_restarts),
            max_iters=int(max_iters),
            stagnation_limit=int(stagnation_limit),
            move_probs=move_probs,
        )

        if chosen is None or not np.isfinite(err):
            continue

        if float(err) < best_err:
            best_chosen = chosen
            best_err = float(err)

        if best_err <= max_err_threshold:
            break

    if best_chosen is None:
        return {"name": name, "ok": False, "err": float("inf"), "chosen_idx": []}

    chosen_idx = np.flatnonzero(best_chosen).astype(np.int64).tolist()
    return {"name": name, "ok": True, "err": float(best_err), "chosen_idx": chosen_idx}


def _config_primary(fracs: dict[str, float]) -> str:
    best = max(
        (("Virus", float(fracs.get("Virus", 0.0))), ("Host", float(fracs.get("Host", 0.0))), ("MGE", float(fracs.get("MGE", 0.0)))),
        key=lambda kv: kv[1],
    )
    return str(best[0])


def _is_near_all(name: str) -> bool:
    return "near_all" in name


def _is_enriched(name: str, fracs: dict[str, float]) -> bool:
    primary = _config_primary(fracs)
    return (not _is_near_all(name)) and (float(fracs.get(primary, 0.0)) >= 0.70)


def build_test_pool_stats(
    test_pool: pl.DataFrame,
    prefilter: Optional[dict] = None,
) -> tuple[pl.DataFrame, np.ndarray, np.ndarray]:
    if test_pool.is_empty():
        empty = pl.DataFrame({"Contig": [], "contig_id": []})
        return empty, np.zeros((0, 3), dtype=np.int64), np.zeros((0,), dtype=np.int64)

    stats = (
        test_pool.group_by("Contig")
        .agg(
            (pl.col("Source") == "Virus").sum().cast(pl.Int64).alias("Virus"),
            (pl.col("Source") == "Host").sum().cast(pl.Int64).alias("Host"),
            (pl.col("Source") == "MGE").sum().cast(pl.Int64).alias("MGE"),
        )
        .with_columns((pl.col("Virus") + pl.col("Host") + pl.col("MGE")).alias("n_total"))
        .sort("Contig")  # group_by output is unordered; sort for deterministic row order
    )

    if prefilter is not None:
        primary = str(prefilter.get("primary", "Virus"))
        min_primary = float(prefilter.get("min_primary", 0.0))
        top_k = int(prefilter.get("top_k", 0))
        if primary not in {"Virus", "Host", "MGE"}:
            primary = "Virus"
        stats = stats.with_columns(
            pl.when(pl.col("n_total") > 0)
            .then(pl.col(primary) / pl.col("n_total"))
            .otherwise(0.0)
            .alias("_pfrac")
        ).filter(pl.col("_pfrac") >= min_primary)
        if top_k > 0:
            stats = stats.sort(["_pfrac", "n_total", "Contig"], descending=[True, True, False]).head(top_k)
        stats = stats.drop("_pfrac")

    keep_contigs = stats.select("Contig")
    df = test_pool.join(keep_contigs, on="Contig", how="inner")

    unique_contigs = df.select("Contig").unique().sort("Contig").with_row_index("contig_id")
    df_wid = df.join(unique_contigs, on="Contig", how="left")

    stats2 = (
        df_wid.group_by("contig_id")
        .agg(
            (pl.col("Source") == "Virus").sum().cast(pl.Int64).alias("Virus"),
            (pl.col("Source") == "Host").sum().cast(pl.Int64).alias("Host"),
            (pl.col("Source") == "MGE").sum().cast(pl.Int64).alias("MGE"),
        )
        .with_columns((pl.col("Virus") + pl.col("Host") + pl.col("MGE")).alias("n_total"))
        .sort("contig_id")
    )

    contrib = stats2.select(["Virus", "Host", "MGE"]).to_numpy().astype(np.int64, copy=False)
    totals = stats2["n_total"].to_numpy().astype(np.int64, copy=False)
    return df_wid, contrib, totals


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Fast deterministic block-atomic train/test + contig-atomic distribution test sets.")
    p.add_argument("--input", required=True, help="Input table (parquet or csv/tsv).")
    p.add_argument("--input-format", default="parquet", choices=["parquet", "csv", "tsv"])
    p.add_argument("--cluster-path", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--seed", type=int, default=20250128)
    p.add_argument("--train-frac", type=float, default=0.90)
    p.add_argument("--train-balance", action="store_true")
    p.add_argument("--n-jobs", type=int, default=10)
    p.add_argument("--target-test-size", type=int, default=1500000)
    p.add_argument("--target-test-pool-frac", type=float, default=0.60)
    p.add_argument("--min-viral-frac", type=float, default=0.40)
    p.add_argument("--max-viral-frac", type=float, default=0.60)

    p.add_argument("--allow-downsize", action="store_true")
    p.add_argument("--downsize-grid", default="1.00,0.95,0.90,0.85,0.80,0.75,0.70,0.60,0.50,0.40,0.30,0.20,0.12,0.10")
    p.add_argument("--total-tol-frac", type=float, default=0.05)
    p.add_argument("--max-err-threshold", type=float, default=0.05)
    p.add_argument("--n-restarts", type=int, default=20)
    p.add_argument("--max-iters", type=int, default=50000)
    p.add_argument("--stagnation-limit", type=int, default=8000)
    p.add_argument("--move-probs", default="0.70,0.15,0.15")
    p.add_argument("--min-target-frac", type=float, default=0.10)
    p.add_argument("--min-pool-frac", type=float, default=0.10)
    return p.parse_args()


def load_input(path: str, fmt: str) -> pl.LazyFrame:
    if fmt == "parquet":
        return pl.scan_parquet(path)
    if fmt == "csv":
        return pl.scan_csv(path, separator=",", infer_schema_length=1000)
    if fmt == "tsv":
        return pl.scan_csv(path, separator="\t", infer_schema_length=1000)
    raise ValueError(f"Unknown input format: {fmt}")


def _config_pool_and_prefilter(
    name: str,
    fracs: dict[str, float],
    test_pool_raw: pl.DataFrame,
    test_pool_vbounded: pl.DataFrame,
    contig_block_map: dict[str, int],
    seed: int,
) -> tuple[pl.DataFrame, Optional[dict]]:
    if name == "test_near_all_mge":
        pool = maximize_test_pool_with_source_bounds(
            test_pool=test_pool_raw,
            contig_block_map=contig_block_map,
            seed=seed + 17,
            source="MGE",
            min_frac=0.70,
            max_frac=0.85,
        )
        prefilter = {"primary": "MGE", "min_primary": 0.50, "top_k": 0}
        return pool, prefilter

    if name == "test_mge_enriched":
        pool = maximize_test_pool_with_source_bounds(
            test_pool=test_pool_raw,
            contig_block_map=contig_block_map,
            seed=seed + 29,
            source="MGE",
            min_frac=0.50,
            max_frac=0.75,
        )
        return pool, None

    if name == "test_near_all_virus":
        pool = maximize_test_pool_with_source_bounds(
            test_pool=test_pool_raw,
            contig_block_map=contig_block_map,
            seed=seed + 19,
            source="Virus",
            min_frac=0.80,
            max_frac=0.90,
        )
        return pool, None

    if name == "test_near_all_host":
        pool = maximize_test_pool_with_source_bounds(
            test_pool=test_pool_raw,
            contig_block_map=contig_block_map,
            seed=seed + 23,
            source="Host",
            min_frac=0.80,
            max_frac=0.90,
        )
        return pool, None

    return test_pool_vbounded, None


def _config_solver_overrides(
    name: str,
    fracs: dict[str, float],
    args: argparse.Namespace,
) -> dict:
    primary = _config_primary(fracs)
    base = {
        "allow_downsize": bool(args.allow_downsize),
        "downsize_grid": [float(x) for x in str(args.downsize_grid).split(",") if str(x).strip()],
        "total_tol_frac": float(args.total_tol_frac),
        "max_err_threshold": float(args.max_err_threshold),
        "n_restarts": int(args.n_restarts),
        "max_iters": int(args.max_iters),
        "stagnation_limit": int(args.stagnation_limit),
        "min_target_frac": float(args.min_target_frac),
        "min_pool_frac": float(args.min_pool_frac),
        "min_abs_keep": 150000,
    }

    if _is_near_all(name):
        base.update(
            {
                "allow_downsize": True,
                "downsize_grid": [1.00, 0.80, 0.60, 0.45, 0.35, 0.25, 0.18, 0.12, 0.08, 0.05, 0.03, 0.02],
                "total_tol_frac": max(float(args.total_tol_frac), 0.07),
                "max_err_threshold": min(float(args.max_err_threshold), 0.05),
                "n_restarts": max(int(args.n_restarts), 14),
                "max_iters": max(int(args.max_iters), 55000),
                "stagnation_limit": max(int(args.stagnation_limit), 15000),
                "min_target_frac": min(float(args.min_target_frac), 0.05),
                "min_pool_frac": min(float(args.min_pool_frac), 0.05),
                "min_abs_keep": 50000,
            }
        )
        return base

    if name == "test_mge_enriched":
        base.update(
            {
                "allow_downsize": True,
                "downsize_grid": [1.00, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.33, 0.25, 0.20, 0.15, 0.12, 0.10],
                "total_tol_frac": max(float(args.total_tol_frac), 0.06),
                "max_err_threshold": min(float(args.max_err_threshold), 0.06),
                "n_restarts": max(int(args.n_restarts), 14),
                "max_iters": max(int(args.max_iters), 50000),
                "stagnation_limit": max(int(args.stagnation_limit), 14000),
                "min_target_frac": max(float(args.min_target_frac), 0.12),
                "min_pool_frac": max(float(args.min_pool_frac), 0.08),
                "min_abs_keep": 75000,
            }
        )
        return base

    if _is_enriched(name, fracs):
        base.update(
            {
                "allow_downsize": True,
                "downsize_grid": [1.00, 0.85, 0.70, 0.60, 0.50, 0.40, 0.33, 0.25, 0.20, 0.15, 0.12, 0.10],
                "total_tol_frac": max(float(args.total_tol_frac), 0.06),
                "max_err_threshold": min(float(args.max_err_threshold), 0.05),
                "n_restarts": max(int(args.n_restarts), 12),
                "max_iters": max(int(args.max_iters), 45000),
                "stagnation_limit": max(int(args.stagnation_limit), 12000),
                "min_target_frac": max(float(args.min_target_frac), 0.15),
                "min_pool_frac": max(float(args.min_pool_frac), 0.10),
                "min_abs_keep": 200000,
            }
        )
        return base

    if primary == "MGE":
        base.update({"min_abs_keep": 100000})

    return base


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    log = logging.getLogger("train_test_split_fast")

    df0 = load_input(args.input, args.input_format)
    schema = df0.collect_schema()
    required_cols = {"Source", "Contig", "contig_pos_start", "contig_pos_end", "region_contig_type"}
    missing = [c for c in required_cols if c not in schema]
    if missing:
        log.error("Missing required columns: %s", missing)
        return 2

    contig_block_map = read_cluster_file_to_contig_block_map(args.cluster_path)
    log.info("Loaded %d contigs from cluster map", len(contig_block_map))

    log.info("Splitting train/test (block-atomic)...")
    train_df, test_pool = contig_train_test_split_block_atomic(
        df=df0,
        contig_block_map=contig_block_map,
        train_frac=float(args.train_frac),
        seed=int(args.seed),
        contig_col="Contig",
    )
    log.info("Initial train proteins: %d", train_df.height)
    log.info("Initial test_pool proteins: %d", test_pool.height)

    train_blocks_df = (
        attach_block_id_fast(train_df.lazy(), contig_block_map=contig_block_map, contig_col="Contig")
        .select("_block_id")
        .unique()
        .sort("_block_id")  # unique() output is unordered; sort so to_list() is deterministic
        .collect()
    )
    train_block_list = train_blocks_df["_block_id"].to_list()

    test_pool_raw = (
        attach_block_id_fast(test_pool.lazy(), contig_block_map=contig_block_map, contig_col="Contig")
        .filter(~pl.col("_block_id").is_in(train_block_list))
        .drop("_block_id")
        .collect()
        .sort(["Source", "Contig", "contig_pos_start", "contig_pos_end"])
    )
    log.info("Leakage-pruned test_pool proteins: %d", test_pool_raw.height)

    log.info("Building provirus test set (chromosome_mixed contigs with both Host and Virus)...")
    test_provirus = build_test_provirus(test_pool_raw)
    n_contigs_provirus = int(test_provirus.select("Contig").n_unique()) if not test_provirus.is_empty() else 0
    counts_p = test_provirus.group_by("Source").agg(pl.len().alias("n")).sort("Source") if not test_provirus.is_empty() else pl.DataFrame({"Source": [], "n": []})
    total_p = int(test_provirus.height) or 1
    c_map_p = {str(src): int(n) for src, n in counts_p.iter_rows()} if not test_provirus.is_empty() else {}
    vp = c_map_p.get("Virus", 0) / total_p
    hp = c_map_p.get("Host", 0) / total_p
    mp = c_map_p.get("MGE", 0) / total_p
    log.info("test_provirus: proteins=%d contigs=%d V=%.3f H=%.3f M=%.3f", int(test_provirus.height), int(n_contigs_provirus), float(vp), float(hp), float(mp))
    outp_provirus = outdir / "test_provirus.parquet"
    test_provirus.write_parquet(outp_provirus)
    log.info("Wrote %s", str(outp_provirus))

    log.info("Bounding base test_pool viral fraction (block-atomic) for non-extreme configs...")
    test_pool_vbounded = maximize_test_pool_with_source_bounds(
        test_pool=test_pool_raw,
        contig_block_map=contig_block_map,
        seed=int(args.seed) + int(_stable_hash_u64("test_pool_bounds_virus") % 1_000_003),
        source="Virus",
        min_frac=float(args.min_viral_frac),
        max_frac=float(args.max_viral_frac),
    )
    log.info("Base v-bounded test_pool proteins: %d", test_pool_vbounded.height)

    if args.train_balance:
        if "True Positive" not in train_df.columns:
            log.error("--train-balance requested but input lacks 'True Positive' column")
            return 2
        log.info("Balancing training set (contig-atomic)...")
        train_df = enforce_balanced_training_set(train_df, seed=int(args.seed), contig_col="Contig", label_col="True Positive")
        log.info("Balanced train proteins: %d", train_df.height)

    train_df = train_df.sort(["Source", "Contig", "contig_pos_start", "contig_pos_end"])
    train_out = outdir / "train.parquet"
    train_df.write_parquet(train_out)
    log.info("Wrote %s", str(train_out))

    src_counts = train_df.group_by("Source").len().rename({"len": "n"}).sort("Source")
    counts_map = {"Virus": 0, "Host": 0, "MGE": 0}
    for src, n in src_counts.iter_rows():
        counts_map[str(src)] = int(n)
    total_train = int(sum(counts_map.values())) or 1
    train_fracs = {k: counts_map[k] / total_train for k in ("Virus", "Host", "MGE")}

    configs: List[Tuple[str, Dict[str, float]]] = [
        ("test_input", train_fracs),
        ("test_equal_pos", {"MGE": 0.25, "Virus": 0.50, "Host": 0.25}),
        ("test_equal_source", {"MGE": 0.333, "Virus": 0.333, "Host": 0.334}),
        ("test_half_virus_host", {"MGE": 0.0, "Virus": 0.50, "Host": 0.50}),
        ("test_virus_enriched", {"MGE": 0.125, "Virus": 0.75, "Host": 0.125}),
        ("test_mge_enriched", {"MGE": 0.75, "Virus": 0.125, "Host": 0.125}),
        ("test_host_enriched", {"MGE": 0.125, "Virus": 0.125, "Host": 0.75}),
        ("test_near_all_virus", {"MGE": 0.05, "Virus": 0.90, "Host": 0.05}),
        # ("test_near_all_mge", {"MGE": 0.90, "Virus": 0.05, "Host": 0.05}),
        ("test_near_all_host", {"MGE": 0.05, "Virus": 0.05, "Host": 0.90}),
        ("test_provirus", {}), # special: written above; not solved to a target composition
    ]

    move_probs = [float(x) for x in str(args.move_probs).split(",") if str(x).strip()]
    if len(move_probs) != 3:
        log.error("--move-probs must have 3 comma-separated floats")
        return 2
    move_probs_t = tuple(move_probs)

    available_total = int(test_pool_raw.height)
    target_test_size = int(min(int(args.target_test_size), int(available_total * float(args.target_test_pool_frac))))
    log.info("Target test size (pre-downsize): %d", target_test_size)

    tmpdir_obj = tempfile.TemporaryDirectory(prefix="train_test_split_mmap_", dir=str(outdir))
    tmpdir = Path(tmpdir_obj.name)

    tasks = []
    config_payloads: dict[str, dict] = {}

    log.info("Preparing per-config pools, contig stats, and memmaps...")
    for idx, (name, fracs) in enumerate(configs):
        if name == "test_provirus":
            continue  # special: already produced; not solved to a target composition

        req = np.array([float(fracs.get(s, 0.0)) for s in ("Virus", "Host", "MGE")], dtype=np.float64)
        if float(req.sum()) <= 0:
            continue
        conf_seed = int(args.seed) + int(_stable_hash_u64(("config", name, idx)) % 1_000_003)

        pool_for_config, prefilter = _config_pool_and_prefilter(
            name=name,
            fracs=fracs,
            test_pool_raw=test_pool_raw,
            test_pool_vbounded=test_pool_vbounded,
            contig_block_map=contig_block_map,
            seed=conf_seed,
        )

        df_wid, contrib, totals = build_test_pool_stats(pool_for_config, prefilter=prefilter)
        log.info("%s: pool proteins=%d contigs=%d", name, int(df_wid.height), int(totals.size))

        if totals.size == 0 or int(totals.sum()) <= 0:
            log.warning("%s: empty stats after pool/prefilter, skipping", name)
            continue

        bundle = _write_memmaps(contrib=contrib, totals=totals, tmpdir=tmpdir, tag=f"{idx}_{name}")
        config_payloads[name] = {"df_wid": df_wid, "bundle": bundle}

        overrides = _config_solver_overrides(name=name, fracs=fracs, args=args)
        tasks.append(
            {
                "name": name,
                "seed": conf_seed,
                "bundle": bundle,
                "req": (req / float(req.sum())).tolist(),
                "target_size": int(target_test_size),
                "allow_downsize": bool(overrides["allow_downsize"]),
                "downsize_grid": list(overrides["downsize_grid"]),
                "total_tol_frac": float(overrides["total_tol_frac"]),
                "max_err_threshold": float(overrides["max_err_threshold"]),
                "n_restarts": int(overrides["n_restarts"]),
                "max_iters": int(overrides["max_iters"]),
                "stagnation_limit": int(overrides["stagnation_limit"]),
                "move_probs": move_probs_t,
                "min_target_frac": float(overrides["min_target_frac"]),
                "min_pool_frac": float(overrides["min_pool_frac"]),
                "min_abs_keep": int(overrides["min_abs_keep"]),
            }
        )

    if not tasks:
        log.error("No runnable configs after preprocessing")
        tmpdir_obj.cleanup()
        return 2

    log.info("Solving %d configs with %d processes...", len(tasks), int(args.n_jobs))
    results = []
    with ProcessPoolExecutor(max_workers=int(args.n_jobs)) as ex:
        futs = [ex.submit(_dataset_worker, t) for t in tasks]
        for fut in as_completed(futs):
            results.append(fut.result())

    order = {name: i for i, (name, _) in enumerate(configs)}
    results = sorted(results, key=lambda r: order.get(r["name"], 10**9))

    for r in results:
        name = r["name"]
        if not r["ok"]:
            log.warning("No solution for %s", name)
            continue

        chosen_idx = np.array(r["chosen_idx"], dtype=np.int64)
        if chosen_idx.size == 0:
            log.warning("Empty selection for %s", name)
            continue

        df_wid = config_payloads[name]["df_wid"]
        df_sel = (
            df_wid.filter(pl.col("contig_id").is_in(chosen_idx.tolist()))
            .drop("contig_id")
            .sort(["Dataset", "Contig", "contig_pos_start", "contig_pos_end"])
        )

        counts = df_sel.group_by("Source").agg(pl.len().alias("n")).sort("Source")
        total = int(df_sel.height) or 1
        c_map = {str(src): int(n) for src, n in counts.iter_rows()}
        vf = c_map.get("Virus", 0) / total
        hf = c_map.get("Host", 0) / total
        mf = c_map.get("MGE", 0) / total

        log.info("%s: proteins=%d err=%.3f V=%.3f H=%.3f M=%.3f", name, int(df_sel.height), float(r["err"]), vf, hf, mf)

        outp = outdir / f"{name}.parquet"
        df_sel.write_parquet(outp)
        log.info("Wrote %s", str(outp))

    tmpdir_obj.cleanup()
    log.info("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())