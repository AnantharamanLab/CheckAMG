#!/usr/bin/env python3

import argparse
import json
import pickle
import pandas as pd
import numpy as np
import glob
import os
import polars as pl
import lightgbm as lgb
from lightgbm import LGBMClassifier
from skopt import BayesSearchCV
import math
import csv
import logging
import warnings

from joblib import Parallel, delayed

from sklearn.model_selection import GroupKFold, GroupShuffleSplit, cross_val_score
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    confusion_matrix, precision_score, recall_score, f1_score, accuracy_score,
    balanced_accuracy_score, matthews_corrcoef, precision_recall_curve,
    auc, roc_curve, average_precision_score, roc_auc_score
)

warnings.filterwarnings("ignore")

log_level = logging.INFO
logging.basicConfig(
    level=log_level,
    format="%(asctime)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger()

try:
    from sklearn.model_selection import StratifiedGroupKFold
    _HAVE_SGKF = True
except Exception:
    StratifiedGroupKFold = None
    _HAVE_SGKF = False

FEATURE_COLS = [
    "Pfam_V-score", "Pfam_VL-score", "KEGG_V-score", "KEGG_VL-score",
    "PHROG_V-score", "PHROG_VL-score",
    "contig_avg_KEGG_V-score", "contig_avg_KEGG_VL-score",
    "contig_avg_Pfam_V-score", "contig_avg_Pfam_VL-score",
    "contig_avg_PHROG_V-score", "contig_avg_PHROG_VL-score",
    "window_avg_KEGG_VL-score", "window_avg_Pfam_VL-score", "window_avg_PHROG_VL-score",
    "KEGG_viral_left_dist", "KEGG_viral_right_dist",
    "Pfam_viral_left_dist", "Pfam_viral_right_dist",
    "PHROG_viral_left_dist", "PHROG_viral_right_dist",
    "KEGG_MGE_left_dist", "KEGG_MGE_right_dist",
    "Pfam_MGE_left_dist", "Pfam_MGE_right_dist",
    "PHROG_MGE_left_dist", "PHROG_MGE_right_dist",
    "KEGG_V-score_left_MGE", "KEGG_V-score_right_MGE",
    "KEGG_VL-score_left_MGE", "KEGG_VL-score_right_MGE",
    "Pfam_V-score_left_MGE", "Pfam_V-score_right_MGE",
    "Pfam_VL-score_left_MGE", "Pfam_VL-score_right_MGE",
    "PHROG_V-score_left_MGE", "PHROG_V-score_right_MGE",
    "PHROG_VL-score_left_MGE", "PHROG_VL-score_right_MGE",
    "circular_contig",
]

ALLOWED_META_COLS = [
    "Protein", "Contig", "Genome", "Protein Classification",
    "gene_number", "contig_pos_start", "contig_pos_end", "length", "frame",
    "True Positive", "Source",
]


def _set_global_thread_env(max_threads: int):
    max_threads = max(1, int(max_threads))
    for k, v in {
        "POLARS_MAX_THREADS": str(max_threads),
        "OMP_NUM_THREADS": str(max_threads),
        "OPENBLAS_NUM_THREADS": str(max_threads),
        "MKL_NUM_THREADS": str(max_threads),
        "VECLIB_MAXIMUM_THREADS": str(max_threads),
        "NUMEXPR_MAX_THREADS": str(max_threads),
    }.items():
        os.environ[k] = v


def sample_by_contig(df, frac=0.08, group_col="Contig", random_state=20260201):
    if frac <= 0 or frac >= 1.0:
        return df.copy()
    rng = np.random.RandomState(random_state)
    contig_sizes = df.groupby(group_col).size().reset_index(name="n")
    contig_sizes = contig_sizes.sample(frac=1.0, random_state=random_state).reset_index(drop=True)
    total = len(df)
    target = int(total * frac)
    picked = []
    acc = 0
    for _, r in contig_sizes.iterrows():
        picked.append(r[group_col])
        acc += int(r["n"])
        if acc >= target:
            break
    out = df[df[group_col].isin(picked)].reset_index(drop=True)
    logger.info(f"Sampled {len(out)} rows (~{len(out)/total:.3%}) by contig (target frac={frac})")
    return out


def _parse_precision_floors(spec_list):
    """
    Parse a list of 'dataset_name:min_precision' strings into a dict.
    E.g. ['test_provirus:0.90', 'test_near_all_host:0.80'] ->
         {'test_provirus': 0.90, 'test_near_all_host': 0.80}
    """
    if not spec_list:
        return {}
    out = {}
    for s in spec_list:
        s = s.strip()
        if not s:
            continue
        parts = s.rsplit(":", 1)
        if len(parts) != 2:
            raise ValueError(
                f"Bad --*-precision-floor spec '{s}': expected 'dataset_name:min_precision' "
                f"(e.g. 'test_provirus:0.90')"
            )
        name = parts[0].strip()
        try:
            val = float(parts[1].strip())
        except ValueError:
            raise ValueError(f"Bad precision value in '{s}': '{parts[1]}' is not a float")
        if not (0.0 <= val <= 1.0):
            raise ValueError(f"Precision floor must be in [0,1], got {val} for dataset '{name}'")
        out[name] = val
    return out


class ViralConfidencePredictor:
    def __init__(
        self,
        max_threads=50,
        reliability_n_bins=15,
        reliability_strategy="uniform",
        bayes_candidates_in_parallel=1,
    ):
        self.model = None
        self.feature_names = None
        self.confidence_thresholds = None
        self.training_results = None
        self.validation_results = None
        self.validation_set = None
        self.training_dataset = None
        self.cleaned_datasets = None

        self.max_threads = max(1, int(max_threads))
        self.bayes_candidates_in_parallel = max(1, int(bayes_candidates_in_parallel))

        self.reliability_n_bins = int(reliability_n_bins)
        self.reliability_strategy = str(reliability_strategy).lower().strip()
        if self.reliability_strategy not in ("uniform", "quantile"):
            logger.warning(f"Unknown reliability_strategy='{reliability_strategy}', defaulting to 'uniform'")
            self.reliability_strategy = "uniform"

        self.cv_jobs = 1
        self.lgb_threads = 1

        _set_global_thread_env(self.max_threads)
        logger.info(
            f"Thread budget: max_threads={self.max_threads}; "
            f"candidates_in_parallel={self.bayes_candidates_in_parallel}"
        )

    def _configure_parallelism_for_cv(self, n_splits: int, candidates_in_parallel: int = 1):
        n_splits = max(1, int(n_splits))
        candidates_in_parallel = max(1, int(candidates_in_parallel))
        desired_jobs = n_splits * candidates_in_parallel
        self.cv_jobs = min(self.max_threads, desired_jobs)
        self.lgb_threads = max(1, self.max_threads // self.cv_jobs)
        while self.cv_jobs * self.lgb_threads > self.max_threads and self.lgb_threads > 1:
            self.lgb_threads -= 1
        logger.info(
            f"Parallelism for CV: n_splits={n_splits}, candidates_in_parallel={candidates_in_parallel} "
            f"-> cv_jobs={self.cv_jobs}, lgb_threads={self.lgb_threads} "
            f"(product={self.cv_jobs * self.lgb_threads} <= {self.max_threads})"
        )

    def prepare_features(self, df):
        return df[FEATURE_COLS]

    def calculate_confusion_metrics(self, y_true, y_pred):
        accuracy = accuracy_score(y_true, y_pred)
        balanced_accuracy = balanced_accuracy_score(y_true, y_pred)
        precision = precision_score(y_true, y_pred, zero_division=0)
        recall = recall_score(y_true, y_pred, zero_division=0)
        f1 = f1_score(y_true, y_pred, zero_division=0)
        mcc = matthews_corrcoef(y_true, y_pred)
        cm = confusion_matrix(y_true, y_pred)
        tn, fp, fn, tp = cm.ravel()
        fdr = fp / (fp + tp) if (fp + tp) > 0 else 0
        specificity = tn / (tn + fp) if (tn + fp) > 0 else 0
        npv = tn / (tn + fn) if (tn + fn) > 0 else 0
        return {
            "accuracy": accuracy, "balanced_accuracy": balanced_accuracy,
            "tpr": recall, "tnr": specificity,
            "precision": precision, "recall": recall, "f1_score": f1,
            "fdr": fdr, "mcc": mcc, "specificity": specificity, "npv": npv,
            "confusion_matrix": cm, "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        }

    def _compute_reliability(self, y_true, y_prob, n_bins=15, strategy="uniform"):
        y_prob = np.asarray(y_prob, dtype=float)
        y_true = np.asarray(y_true, dtype=int)
        n = y_prob.size
        if n == 0:
            return {
                "prob_pred": [], "prob_true": [], "bin_edges": [],
                "bin_counts": [], "ece": np.nan, "n_bins_used": 0, "strategy": strategy,
            }
        if strategy == "quantile":
            q = np.linspace(0.0, 1.0, n_bins + 1)
            edges = np.quantile(y_prob, q)
            edges = np.unique(edges)
            if edges.size < 3:
                edges = np.linspace(0.0, 1.0, n_bins + 1)
                strategy_used = "uniform_fallback"
            else:
                strategy_used = "quantile"
        else:
            edges = np.linspace(0.0, 1.0, n_bins + 1)
            strategy_used = "uniform"
        idx = np.digitize(y_prob, edges, right=True) - 1
        idx = np.clip(idx, 0, edges.size - 2)
        prob_pred, prob_true, counts = [], [], []
        ece_num = 0.0
        for b in range(edges.size - 1):
            mask = idx == b
            cnt = int(np.sum(mask))
            if cnt == 0:
                continue
            conf = float(y_prob[mask].mean())
            acc = float(y_true[mask].mean())
            prob_pred.append(conf)
            prob_true.append(acc)
            counts.append(cnt)
            ece_num += abs(acc - conf) * cnt
        ece = ece_num / float(n) if n > 0 else np.nan
        return {
            "prob_pred": np.asarray(prob_pred, dtype=float),
            "prob_true": np.asarray(prob_true, dtype=float),
            "bin_edges": np.asarray(edges, dtype=float),
            "bin_counts": np.asarray(counts, dtype=int),
            "ece": float(ece),
            "n_bins_used": int(len(prob_pred)),
            "strategy": strategy_used,
        }

    def make_sample_weights(self, df, w_neg_mge=1.0, w_neg_host=1.0, w_pos_virus=1.0):
        logger.info(
            f"Using sample weights: w_neg_mge={w_neg_mge}, w_neg_host={w_neg_host}, w_pos_virus={w_pos_virus}"
        )
        y = df["True Positive"].astype(int).values
        src = df["Source"].values
        w = np.ones(len(df), dtype=float)
        w[y == 1] = w_pos_virus
        w[(y == 0) & (src == "Host")] = w_neg_host
        w[(y == 0) & (src == "MGE")] = w_neg_mge
        return w

    def _fit_base_lgbm(
        self, X, y,
        best_params=None,
        random_state=20260201,
        sample_weight=None,
        groups=None,
        X_val=None, y_val=None, w_val=None,
        use_early_stopping=True,
        early_stopping_rounds=200,
    ):
        if _HAVE_SGKF:
            cv_splitter = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=random_state)
        else:
            cv_splitter = GroupKFold(n_splits=5)
            logger.info("Note: StratifiedGroupKFold not available; using GroupKFold.")

        try:
            n_splits = cv_splitter.get_n_splits(groups=groups)
        except Exception:
            n_splits = getattr(cv_splitter, "n_splits", 5)

        self._configure_parallelism_for_cv(n_splits, candidates_in_parallel=self.bayes_candidates_in_parallel)

        if best_params:
            cfg = None
            best_score = None
            if isinstance(best_params, dict):
                if "best_config" in best_params and isinstance(best_params["best_config"], dict):
                    cfg = best_params["best_config"]
                    best_score = best_params.get("best_score") or best_params.get("cv_score")
                elif "model_config" in best_params and isinstance(best_params["model_config"], dict):
                    cfg = best_params["model_config"]
                    best_score = best_params.get("best_score") or best_params.get("cv_score")
                elif "best_params" in best_params and isinstance(best_params["best_params"], dict):
                    cfg = best_params["best_params"]
                    best_score = best_params.get("best_score") or best_params.get("cv_score")
                else:
                    cfg = {k: v for k, v in best_params.items()
                           if k not in ("best_score", "cv_score", "best_config", "model_config", "best_params")}
                    best_score = best_params.get("best_score") or best_params.get("cv_score")
            if cfg is None:
                raise ValueError(
                    "best_params provided but could not extract a param dict. "
                    "Pass {'best_config': {...}, 'best_score': ...} or a raw param dict."
                )
            cfg = {str(k): (v.item() if hasattr(v, "item") else v) for k, v in cfg.items()}
            logger.info(f"Using provided best params (n_keys={len(cfg)}); best_score={best_score}")
            try:
                base = LGBMClassifier(
                    **cfg,
                    objective="binary",
                    metric="binary_logloss",
                    n_jobs=self.max_threads,
                    random_state=random_state,
                    verbosity=-1,
                    use_missing=True,
                )
                fit_kwargs = {}
                if use_early_stopping and (X_val is not None) and (y_val is not None):
                    fit_kwargs.update({
                        "eval_set": [(X_val, y_val)],
                        "eval_metric": "aucpr",
                        "eval_sample_weight": [w_val] if w_val is not None else None,
                        "callbacks": [lgb.early_stopping(early_stopping_rounds, verbose=False)],
                    })
                base.fit(X, y, sample_weight=sample_weight, **fit_kwargs)
            except TypeError as e:
                logger.error("LightGBM TypeError. Likely an invalid hyperparameter name.")
                logger.error(f"Error: {e}")
                logger.error(f"Keys passed to LGBMClassifier: {list(cfg.keys())}")
                raise
            return base, cfg, best_score

        lgb_est = LGBMClassifier(
            objective="binary", metric="binary_logloss", boosting_type="gbdt",
            n_jobs=self.lgb_threads, random_state=random_state, verbosity=-1, use_missing=True,
        )
        search_space = {
            "num_leaves": (16, 256),
            "max_depth": (-1, 12),
            "learning_rate": (0.01, 0.2, "log-uniform"),
            "n_estimators": (300, 1200),
            "subsample": (0.5, 1.0, "uniform"),
            "colsample_bytree": (0.4, 1.0, "uniform"),
            "subsample_freq": (0, 10),
            "reg_lambda": (1e-3, 5.0, "log-uniform"),
            "reg_alpha": (1e-3, 2.0, "log-uniform"),
            "min_gain_to_split": (0.0, 1.0, "uniform"),
            "max_bin": (63, 255),
        }
        pre_dispatch = min(2 * self.cv_jobs, 64)
        bayes_cv = BayesSearchCV(
            estimator=lgb_est,
            search_spaces=search_space,
            n_iter=40,
            cv=cv_splitter,
            scoring="average_precision",
            n_jobs=self.lgb_threads,
            pre_dispatch=pre_dispatch,
            verbose=0,
            random_state=random_state,
        )
        fit_kwargs = {}
        if sample_weight is not None:
            fit_kwargs["sample_weight"] = sample_weight
        bayes_cv.fit(X, y, groups=groups, **fit_kwargs)
        logger.info(f"Best CV-AP: {bayes_cv.best_score_:.3f}")
        logger.info(f"Best params: {bayes_cv.best_params_}")
        best_est = bayes_cv.best_estimator_
        best_est.set_params(n_jobs=self.max_threads)
        fit_kwargs2 = {"sample_weight": sample_weight} if sample_weight is not None else {}
        if use_early_stopping and (X_val is not None) and (y_val is not None):
            fit_kwargs2.update({
                "eval_set": [(X_val, y_val)],
                "eval_metric": "aucpr",
                "eval_sample_weight": [w_val] if w_val is not None else None,
                "callbacks": [lgb.early_stopping(early_stopping_rounds, verbose=False)],
            })
        best_est.fit(X, y, **fit_kwargs2)
        return best_est, bayes_cv.best_params_, bayes_cv.best_score_

    def _wilson_lower_bound(self, k, n, z=1.96):
        if n == 0:
            return np.nan
        phat = k / n
        denom = 1 + (z ** 2) / n
        centre = phat + (z ** 2) / (2 * n)
        margin = z * np.sqrt((phat * (1 - phat) + (z ** 2) / (4 * n)) / n)
        return (centre - margin) / denom

    def _dataset_metrics_at_t(self, y_true, y_proba, t):
        pred = (y_proba >= t)
        k = int(pred.sum())
        tp = int(((y_true == 1) & pred).sum())
        P = int((y_true == 1).sum())
        precision = tp / k if k > 0 else 0.0
        recall = tp / P if P > 0 else 0.0
        prec_lb = self._wilson_lower_bound(tp, k, z=1.96) if k > 0 else 0.0
        return {"k": k, "tp": tp, "P": P, "precision": precision, "recall": recall, "precision_lb": prec_lb}

    # ------------------------------------------------------------------
    # THRESHOLD SELECTION
    #
    # Hard precision floors per dataset (precision_floors) act as a
    # precondition: a threshold t is only *eligible* if every named
    # dataset meets its floor.  Among eligible thresholds the existing
    # n_pass / mean_recall objective runs unchanged.
    #
    # This is the right abstraction for provirus / near_all_host because
    # those datasets live in a different operating regime (low positive
    # prevalence, local-context signal) and should never be silently
    # ignored in favour of majority-dataset pass counts.
    # ------------------------------------------------------------------

    def _precision_floors_met(self, t, evaluation_results, precision_floors):
        """Return True iff threshold t satisfies every precision floor constraint."""
        if not precision_floors:
            return True
        for ds_name, min_prec in precision_floors.items():
            if ds_name not in evaluation_results:
                logger.warning(
                    f"precision_floors: dataset '{ds_name}' not found in evaluation_results; "
                    f"floor ignored for this dataset."
                )
                continue
            y = np.asarray(evaluation_results[ds_name]["true_labels"])
            p = np.asarray(evaluation_results[ds_name]["probabilities"])
            pred = p >= t
            k = int(pred.sum())
            tp = int(((y == 1) & pred).sum())
            prec = tp / k if k > 0 else 0.0
            if prec < min_prec:
                return False
        return True

    def _score_threshold_over_datasets(
        self,
        evaluation_results,
        t,
        r_target,
        p_target,
        min_pos_calls=100,
        use_wilson=True,
    ):
        per_ds = {}
        passes = []
        for name, res in evaluation_results.items():
            y = np.asarray(res["true_labels"])
            p = np.asarray(res["probabilities"])
            m = self._dataset_metrics_at_t(y, p, t)
            prec_ok = (m["precision_lb"] >= p_target) if use_wilson else (m["precision"] >= p_target)
            rec_ok = (m["recall"] >= r_target)
            k_ok = (m["k"] >= min_pos_calls) if min_pos_calls is not None else True
            ok = prec_ok and rec_ok and k_ok
            per_ds[name] = {**m, "pass": ok}
            passes.append(ok)
        n_pass = int(np.sum(passes))
        mean_recall_on_pass = float(
            np.mean([per_ds[n]["recall"] for n in per_ds if per_ds[n]["pass"]])
        ) if n_pass > 0 else 0.0
        return n_pass, mean_recall_on_pass, per_ds

    def auto_pick_cutoffs_across_datasets(
        self,
        evaluation_results,
        high_targets=(0.50, 0.95),
        med_targets=(0.75, 0.90),
        grid_step=0.001,
        min_pos_calls=100,
        use_wilson=True,
        # ------------------------------------------------------------------
        # Hard precision floors: {dataset_name: min_precision_required}
        # A threshold is *ineligible* if it violates any floor.
        # Separate floors for high- and medium-confidence tiers so you can
        # apply tighter constraints to the high tier.
        # ------------------------------------------------------------------
        high_precision_floors=None,   # e.g. {"test_provirus": 0.90}
        med_precision_floors=None,    # e.g. {"test_provirus": 0.85}
    ):
        """
        Pick high- and medium-confidence thresholds.

        For each tier the search space is first narrowed to thresholds that
        satisfy all precision_floors constraints for that tier (hard gates).
        Among eligible thresholds the original n_pass / mean_recall objective
        selects the best candidate.  If no threshold satisfies the floors,
        the method falls back to the unconstrained search and logs a warning.
        """
        thresholds = np.arange(0.0, 1.0 + 1e-12, grid_step)
        high_precision_floors = high_precision_floors or {}
        med_precision_floors = med_precision_floors or {}

        # Log active floor constraints
        if high_precision_floors:
            logger.info(
                f"High-tier precision floors: "
                + ", ".join(f"{k}>={v:.3f}" for k, v in high_precision_floors.items())
            )
        if med_precision_floors:
            logger.info(
                f"Medium-tier precision floors: "
                + ", ".join(f"{k}>={v:.3f}" for k, v in med_precision_floors.items())
            )

        def _choose(thr_list, r_target, p_target, precision_floors):
            best = {"t": None, "n_pass": -1, "mean_recall_on_pass": -1.0, "per_dataset": None}

            eligible = [
                t for t in thr_list
                if self._precision_floors_met(t, evaluation_results, precision_floors)
            ]

            if len(eligible) == 0 and precision_floors:
                logger.warning(
                    "No threshold satisfies the precision floors "
                    f"({precision_floors}). Falling back to unconstrained search."
                )
                eligible = list(thr_list)
            elif precision_floors:
                logger.info(
                    f"Precision floors reduce eligible thresholds: "
                    f"{len(eligible)}/{len(thr_list)} candidates remain "
                    f"(range [{eligible[0]:.3f}, {eligible[-1]:.3f}])"
                )

            for t in eligible:
                n_pass, mean_rec, per_ds = self._score_threshold_over_datasets(
                    evaluation_results, t, r_target, p_target,
                    min_pos_calls=min_pos_calls, use_wilson=use_wilson,
                )
                if (
                    n_pass > best["n_pass"]
                    or (n_pass == best["n_pass"] and mean_rec > best["mean_recall_on_pass"])
                    or (
                        n_pass == best["n_pass"]
                        and np.isclose(mean_rec, best["mean_recall_on_pass"])
                        and (best["t"] is None or t < best["t"])
                    )
                ):
                    best = {
                        "t": float(t),
                        "n_pass": n_pass,
                        "mean_recall_on_pass": float(mean_rec),
                        "per_dataset": per_ds,
                    }
            return best

        high = _choose(thresholds, *high_targets, high_precision_floors)
        med = _choose(thresholds, *med_targets, med_precision_floors)

        # Tighten: if all datasets pass for high, walk to the minimum feasible t
        total_ds = len(evaluation_results)
        if high["n_pass"] == total_ds:
            feas = [
                t for t in thresholds
                if self._precision_floors_met(t, evaluation_results, high_precision_floors)
                and self._score_threshold_over_datasets(
                    evaluation_results, t, *high_targets,
                    min_pos_calls=min_pos_calls, use_wilson=use_wilson,
                )[0] == total_ds
            ]
            if feas:
                high_t = min(feas)
                n_pass, mean_rec, per_ds = self._score_threshold_over_datasets(
                    evaluation_results, high_t, *high_targets,
                    min_pos_calls=min_pos_calls, use_wilson=use_wilson,
                )
                high = {
                    "t": float(high_t),
                    "n_pass": n_pass,
                    "mean_recall_on_pass": float(mean_rec),
                    "per_dataset": per_ds,
                }

        # Log per-dataset precision at chosen thresholds so floors are auditable
        if high_precision_floors and high["t"] is not None:
            logger.info(f"High threshold {high['t']:.3f} — precision on floor-constrained datasets:")
            for ds_name, min_prec in high_precision_floors.items():
                if ds_name in evaluation_results:
                    y = np.asarray(evaluation_results[ds_name]["true_labels"])
                    p = np.asarray(evaluation_results[ds_name]["probabilities"])
                    m = self._dataset_metrics_at_t(y, p, high["t"])
                    status = "OK" if m["precision"] >= min_prec else "FLOOR_VIOLATED"
                    logger.info(
                        f"  {ds_name}: precision={m['precision']:.3f} (floor>={min_prec:.3f}) "
                        f"recall={m['recall']:.3f}  [{status}]"
                    )

        if med_precision_floors and med["t"] is not None:
            logger.info(f"Medium threshold {med['t']:.3f} — precision on floor-constrained datasets:")
            for ds_name, min_prec in med_precision_floors.items():
                if ds_name in evaluation_results:
                    y = np.asarray(evaluation_results[ds_name]["true_labels"])
                    p = np.asarray(evaluation_results[ds_name]["probabilities"])
                    m = self._dataset_metrics_at_t(y, p, med["t"])
                    status = "OK" if m["precision"] >= min_prec else "FLOOR_VIOLATED"
                    logger.info(
                        f"  {ds_name}: precision={m['precision']:.3f} (floor>={min_prec:.3f}) "
                        f"recall={m['recall']:.3f}  [{status}]"
                    )

        return {"high": high, "medium": med}

    def train_model(
        self,
        datasets_dict,
        validation_frac=0.20,
        best_params=None,
        w_neg_mge=1.0, w_neg_host=1.0, w_pos_virus=1.0,
        group_col="Contig",
        use_early_stopping=True,
        early_stopping_rounds=200,
    ):
        logger.info("TRAINING VIRAL CONFIDENCE LGBM MODEL")
        if "train" not in datasets_dict:
            raise ValueError("Training dataset with key 'train' not found in datasets_dict")
        training_df = datasets_dict["train"]
        self.training_dataset = "train"

        if group_col not in training_df.columns:
            raise ValueError(f"Grouping column '{group_col}' not in training dataframe.")

        gss = GroupShuffleSplit(n_splits=1, test_size=validation_frac, random_state=20260201)
        groups_all = training_df[group_col].values
        idx = np.arange(len(training_df))
        fit_pos, val_pos = next(gss.split(idx, groups=groups_all))
        df_fit = training_df.iloc[fit_pos].reset_index(drop=True)
        df_val = training_df.iloc[val_pos].reset_index(drop=True)

        X_fit = self.prepare_features(df_fit)
        y_fit = df_fit["True Positive"]
        X_val = self.prepare_features(df_val)
        y_val = df_val["True Positive"]
        self.feature_names = X_fit.columns.tolist()

        logger.info(f"Fit: {len(X_fit)} | Val: {len(X_val)}")
        logger.info(
            f"Class (fit): {y_fit.value_counts().to_dict()} | (val): {y_val.value_counts().to_dict()}"
        )

        w_fit = self.make_sample_weights(
            df_fit, w_neg_mge=w_neg_mge, w_neg_host=w_neg_host, w_pos_virus=w_pos_virus
        )
        w_val = self.make_sample_weights(
            df_val, w_neg_mge=w_neg_mge, w_neg_host=w_neg_host, w_pos_virus=w_pos_virus
        )

        base_model, best_config, best_score = self._fit_base_lgbm(
            X_fit, y_fit,
            best_params=best_params,
            random_state=20260201,
            sample_weight=w_fit,
            groups=df_fit[group_col].values,
            X_val=X_val, y_val=y_val, w_val=w_val,
            use_early_stopping=use_early_stopping,
            early_stopping_rounds=early_stopping_rounds,
        )

        importance_df = pd.DataFrame({
            "feature": self.feature_names,
            "importance": base_model.feature_importances_,
        }).sort_values("importance", ascending=False)

        self.model = base_model
        self.training_results = {
            "model_config": best_config,
            "cv_score": best_score,
            "feature_importance": importance_df,
            "training_composition": {
                "n_total": len(df_fit),
                "n_viral": int(y_fit.sum()),
                "viral_fraction": float(y_fit.mean()),
            },
            "validation_composition": {
                "n_total": len(X_val),
                "n_viral": int(y_val.sum()),
                "viral_fraction": float(y_val.mean()),
            },
        }

        logger.info("LGBM MODEL SUMMARY")
        logger.info(f"Model configuration: {best_config}")
        logger.info(
            f"Best CV AUPRC score (base): {best_score:.3f}"
            if best_score is not None else "Best CV AUPRC score (base): N/A"
        )
        logger.info(f"Number of features: {len(self.feature_names)}")
        logger.info("Top 15 Most Important Features:")
        for i, (_, row) in enumerate(importance_df.head(15).iterrows()):
            logger.info(f"  {i+1:2d}. {row['feature']:30s}: {row['importance']:.3f}")
        return self.training_results

    def calculate_confidence_performance(self, y_true, y_pred, level_name, threshold_used=None):
        m = self.calculate_confusion_metrics(y_true, y_pred)
        m["level"] = level_name
        m["n_predictions"] = int(y_pred.sum())
        m["threshold"] = float(threshold_used) if threshold_used is not None else np.nan
        return m

    def evaluate_on_all_datasets(self, datasets_dict):
        logger.info("EVALUATING MODEL ON ALL DATASETS")
        if self.model is None:
            raise ValueError("No model available. Train the model first.")

        evaluation_results = {}
        skip_keys = {"train", "validation"}

        for dataset_name, df in datasets_dict.items():
            if dataset_name in skip_keys:
                continue
            logger.info(f"Evaluating on {dataset_name}")

            X_test = self.prepare_features(df)
            y_test = df["True Positive"].astype(int).values
            p = self.model.predict_proba(X_test)[:, 1]

            results = {
                "dataset_name": dataset_name,
                "n_total": int(len(y_test)),
                "n_viral": int(y_test.sum()),
                "viral_fraction": float(np.mean(y_test)),
                "probabilities": p,
                "true_labels": y_test,
            }

            auprc = float(average_precision_score(y_test, p))
            auroc = float(roc_auc_score(y_test, p))
            prec, rec, _ = precision_recall_curve(y_test, p)
            fpr, tpr, _ = roc_curve(y_test, p)

            results["pr_curve"] = {"precision": prec, "recall": rec, "auprc": float(auprc)}
            results["roc_curve"] = {"fpr": fpr, "tpr": tpr, "auroc": float(auroc)}

            rel = self._compute_reliability(
                y_true=y_test, y_prob=p,
                n_bins=self.reliability_n_bins, strategy=self.reliability_strategy,
            )
            results["reliability"] = rel
            logger.info(
                f"Calibration ECE (strategy={rel['strategy']}, bins_used={rel['n_bins_used']}): {rel['ece']:.4f}"
            )

            if self.confidence_thresholds is not None:
                t_high = self.confidence_thresholds["high"]["threshold"]
                t_med = self.confidence_thresholds["medium"]["threshold"]
                t_low = self.confidence_thresholds["low"]["threshold"]

                res_high = self.calculate_confidence_performance(
                    y_test, (p >= t_high), "high", threshold_used=t_high
                )
                res_med = self.calculate_confidence_performance(
                    y_test, (p >= t_med), "medium", threshold_used=t_med
                )
                res_low = self.calculate_confidence_performance(
                    y_test, (p >= t_low), "low", threshold_used=t_low
                )

                results["high_confidence"] = res_high
                results["medium_confidence"] = res_med
                results["low_confidence"] = res_low

                for conf_level in ["high_confidence", "medium_confidence", "low_confidence"]:
                    m = results[conf_level]
                    logger.info(f"{conf_level.replace('_', ' ').title()}:")
                    logger.info(
                        f"  Threshold: {m['threshold']:.3f} | Predictions: {m['n_predictions']}"
                    )
                    logger.info(
                        f"  Accuracy: {m['accuracy']:.3f} | Balanced Accuracy: {m['balanced_accuracy']:.3f}"
                    )
                    logger.info(
                        f"  Precision: {m['precision']:.3f} | Recall: {m['recall']:.3f} | F1: {m['f1_score']:.3f}"
                    )
                    logger.info(f"  FDR: {m['fdr']:.3f} | MCC: {m['mcc']:.3f}")

            evaluation_results[dataset_name] = results

        return evaluation_results

    def check_fdr_target(self, conf_level, fdr_value):
        targets = {"high_confidence": 0.05, "medium_confidence": 0.1, "low_confidence": 0.1}
        return fdr_value <= targets.get(conf_level, 1.0)

    def check_recall_target(self, conf_level, recall_value):
        targets = {"high_confidence": 0.50, "medium_confidence": 0.75, "low_confidence": 1.0}
        return recall_value >= targets.get(conf_level, 0.0)

    def generate_performance_summary_table(self, evaluation_results):
        logger.info("PERFORMANCE SUMMARY TABLE")
        if self.confidence_thresholds is None:
            logger.info("Thresholds not set; summary table will only be meaningful after auto-pick.")
        summary_rows = []

        for dataset_name, results in evaluation_results.items():
            df = self.cleaned_datasets.get(dataset_name)
            if df is not None and "Source" in df.columns:
                n_total = len(df)
                viral_frac = (df["Source"] == "Virus").sum() / n_total if n_total > 0 else np.nan
                host_frac = (df["Source"] == "Host").sum() / n_total if n_total > 0 else np.nan
                mge_frac = (df["Source"] == "MGE").sum() / n_total if n_total > 0 else np.nan
            else:
                viral_frac, host_frac, mge_frac = np.nan, np.nan, np.nan

            for k in ["high_confidence", "medium_confidence", "low_confidence"]:
                if k not in results:
                    continue
                metrics = results[k]
                summary_rows.append({
                    "Dataset": dataset_name,
                    "Viral_Fraction_num": viral_frac,
                    "MGE_Fraction_num": mge_frac,
                    "Host_Fraction_num": host_frac,
                    "Confidence_Level": k.replace("_", " ").title(),
                    "pseudo_conf_sort": {"high_confidence": 0, "medium_confidence": 1, "low_confidence": 2}[k],
                    "Threshold": f"{metrics['threshold']:.3f}",
                    "Predictions_Made": metrics["n_predictions"],
                    "Accuracy": f"{metrics['accuracy']:.3f}",
                    "Balanced_Accuracy": f"{metrics['balanced_accuracy']:.3f}",
                    "TPR": f"{metrics['tpr']:.3f}",
                    "TNR": f"{metrics['tnr']:.3f}",
                    "Precision": f"{metrics['precision']:.3f}",
                    "Recall": f"{metrics['recall']:.3f}",
                    "F1_Score": f"{metrics['f1_score']:.3f}",
                    "FDR": f"{metrics['fdr']:.3f}",
                    "MCC": f"{metrics['mcc']:.3f}",
                    "Meets_FDR_Target": self.check_fdr_target(k, metrics["fdr"]),
                    "Meets_Recall_Target": self.check_recall_target(k, metrics["recall"]),
                })

        if not summary_rows:
            logger.info("No confidence metrics available to summarize yet.")
            return pd.DataFrame()

        summary_df = pd.DataFrame(summary_rows)
        summary_df = summary_df.sort_values(
            by=["Viral_Fraction_num", "MGE_Fraction_num", "Host_Fraction_num", "pseudo_conf_sort"],
            ascending=[True, True, True, True],
        ).reset_index(drop=True)

        summary_df["Viral_Fraction"] = summary_df["Viral_Fraction_num"].apply(
            lambda x: f"{x:.1%}" if pd.notnull(x) else "NA"
        )
        summary_df["MGE_Fraction"] = summary_df["MGE_Fraction_num"].apply(
            lambda x: f"{x:.1%}" if pd.notnull(x) else "NA"
        )
        summary_df["Host_Fraction"] = summary_df["Host_Fraction_num"].apply(
            lambda x: f"{x:.1%}" if pd.notnull(x) else "NA"
        )

        summary_df = summary_df.drop(
            columns=["Viral_Fraction_num", "MGE_Fraction_num", "Host_Fraction_num", "pseudo_conf_sort"]
        )
        summary_df = summary_df[[
            "Dataset", "Viral_Fraction", "MGE_Fraction", "Host_Fraction",
            "Confidence_Level", "Threshold", "Predictions_Made",
            "Accuracy", "Balanced_Accuracy", "TPR", "TNR",
            "Precision", "Recall", "F1_Score", "FDR", "MCC",
            "Meets_FDR_Target", "Meets_Recall_Target",
        ]]

        logger.info(f"\n{summary_df.to_string(index=False)}")
        return summary_df

    def get_model_summary(self):
        if self.model is None:
            raise ValueError("Must train model first!")
        return {
            "training_dataset": self.training_dataset,
            "feature_names": self.feature_names,
            "confidence_thresholds": (
                None if self.confidence_thresholds is None else {
                    "high": self.confidence_thresholds["high"]["threshold"],
                    "medium": self.confidence_thresholds["medium"]["threshold"],
                    "low": self.confidence_thresholds["low"]["threshold"],
                }
            ),
            "model_config": self.training_results["model_config"],
            "cv_score": self.training_results["cv_score"],
            "feature_importance": self.training_results["feature_importance"],
        }

    def diagnostic_report(
        self, eval_df, train_df=None, top_k=5, n_shuffle=5,
        sample=50000, random_state=20260201,
    ):
        rng = np.random.RandomState(random_state)
        X_eval = self.prepare_features(eval_df)
        y_eval = eval_df["True Positive"].astype(int).values
        model = self.model
        if model is None:
            raise ValueError("No model available. Train the model first.")

        p = model.predict_proba(X_eval)[:, 1]
        prec, rec, _ = precision_recall_curve(y_eval, p)
        auprc = auc(rec, prec)
        logger.info(f"PR-curve AUPRC (raw probs): {auprc:.3f}")

        available_threads = int(getattr(self, "max_threads", 50) or 50)
        model_threads = None
        try:
            model_threads = int(getattr(model, "n_jobs", None)) if getattr(model, "n_jobs", None) is not None else None
        except Exception:
            model_threads = None

        if model_threads is None:
            try:
                booster = getattr(model, "booster_", None)
                if booster is not None:
                    params = getattr(booster, "params", {})
                    if isinstance(params, dict) and "num_threads" in params:
                        model_threads = int(params.get("num_threads") or 1)
            except Exception:
                model_threads = None

        if model_threads is None:
            model_threads = available_threads
        model_threads = max(1, int(model_threads))

        max_perm_candidates = 8
        perm_jobs = max(1, min(max_perm_candidates, available_threads // model_threads))
        if perm_jobs * model_threads > available_threads:
            perm_jobs = 1

        perm = permutation_importance(
            model, X_eval, y_eval, n_repeats=5, random_state=random_state, n_jobs=perm_jobs
        )
        idx = perm.importances_mean.argsort()[::-1][:top_k]
        logger.info(f"Top-{top_k} permutation importance:")
        for i in idx:
            logger.info(f"{self.feature_names[i]:30s}: {perm.importances_mean[i]:.4f}")

        if train_df is None:
            logger.info("Ablation skipped (train_df not provided)")
            return

        n_total = len(train_df)
        if sample is None or sample >= n_total:
            df_sample = train_df
        else:
            frac = sample / n_total
            df_sample = (
                train_df
                .groupby("True Positive", group_keys=False)
                .apply(lambda g: g.sample(frac=frac, random_state=random_state))
            )

        X_full = self.prepare_features(df_sample)
        y_full = df_sample["True Positive"].astype(int).values
        groups = df_sample["Contig"].values

        cv_folds = 3
        cv_obj = (
            StratifiedGroupKFold(n_splits=cv_folds, shuffle=True, random_state=random_state)
            if _HAVE_SGKF else GroupKFold(n_splits=cv_folds)
        )
        cv_jobs = max(1, min(available_threads, cv_folds))
        threads_per_fit = max(1, available_threads // cv_jobs)
        if model_threads > threads_per_fit:
            cv_jobs = max(1, available_threads // model_threads)
            threads_per_fit = max(1, available_threads // cv_jobs)

        base = LGBMClassifier(
            objective="binary", n_estimators=200, learning_rate=0.1, num_leaves=64,
            n_jobs=threads_per_fit, random_state=random_state,
            scale_pos_weight=(len(y_full) - y_full.sum()) / max(y_full.sum(), 1), verbosity=-1,
        )
        try:
            f1_base = cross_val_score(
                base, X_full, y_full, cv=cv_obj, groups=groups, scoring="f1", n_jobs=cv_jobs
            ).mean()
        except Exception as e:
            logger.warning(f"Cross-val failed with n_jobs={cv_jobs}: {e}; retrying n_jobs=1")
            f1_base = cross_val_score(
                base, X_full, y_full, cv=cv_obj, groups=groups, scoring="f1", n_jobs=1
            ).mean()

        X_drop = X_full.drop(columns=[self.feature_names[i] for i in idx])
        try:
            f1_drop = cross_val_score(
                base, X_drop, y_full, cv=cv_obj, groups=groups, scoring="f1", n_jobs=cv_jobs
            ).mean()
        except Exception as e:
            logger.warning(f"Cross-val (drop) failed with n_jobs={cv_jobs}: {e}; retrying n_jobs=1")
            f1_drop = cross_val_score(
                base, X_drop, y_full, cv=cv_obj, groups=groups, scoring="f1", n_jobs=1
            ).mean()

        logger.info(
            f"Ablation test (drop top {top_k} feats): "
            f"F1 full {f1_base:.3f} -> F1 ablated {f1_drop:.3f} Delta {f1_drop - f1_base:+.3f}"
        )

    def run_complete_analysis(
        self,
        train_test_datasets_dict,
        best_params=None,
        w_neg_mge=1.0, w_neg_host=1.0, w_pos_virus=1.0,
        group_col="Contig",
        use_early_stopping=True,
        early_stopping_rounds=200,
        high_recall_target=0.50,
        high_precision_target=0.95,
        med_recall_target=0.75,
        med_precision_target=0.90,
        high_on_mge_datasets=False,
        # --- new: hard precision floors per dataset, per tier ---
        high_precision_floors=None,   # dict[str, float]
        med_precision_floors=None,    # dict[str, float]
    ):
        logger.info("RUNNING ANALYSIS")
        self.cleaned_datasets = train_test_datasets_dict

        training_results = self.train_model(
            train_test_datasets_dict,
            best_params=best_params,
            w_neg_mge=w_neg_mge, w_neg_host=w_neg_host, w_pos_virus=w_pos_virus,
            group_col=group_col,
            use_early_stopping=use_early_stopping,
            early_stopping_rounds=early_stopping_rounds,
        )

        evaluation_results = self.evaluate_on_all_datasets(train_test_datasets_dict)

        # Select threshold datasets for high tier
        if high_on_mge_datasets:
            mge_keys = [k for k in evaluation_results.keys() if "mge" in k.lower()]
            if mge_keys:
                eval_for_high = {k: evaluation_results[k] for k in mge_keys}
                logger.info(f"High threshold derived from MGE datasets only: {mge_keys}")
            else:
                eval_for_high = evaluation_results
                logger.info("No datasets matched 'mge'; using all datasets for high threshold.")
        else:
            eval_for_high = evaluation_results

        cutoffs_all = self.auto_pick_cutoffs_across_datasets(
            evaluation_results,
            high_targets=(high_recall_target, high_precision_target),
            med_targets=(med_recall_target, med_precision_target),
            grid_step=0.001,
            min_pos_calls=100,
            use_wilson=True,
            high_precision_floors=high_precision_floors or {},
            med_precision_floors=med_precision_floors or {},
        )

        high_t = cutoffs_all["high"]["t"]
        high_metrics = cutoffs_all["high"]

        if high_on_mge_datasets:
            mge_keys = [k for k in evaluation_results.keys() if "mge" in k.lower()]
            if mge_keys:
                sub = {k: evaluation_results[k] for k in mge_keys}
                cutoffs_mge = self.auto_pick_cutoffs_across_datasets(
                    sub,
                    high_targets=(high_recall_target, high_precision_target),
                    med_targets=(med_recall_target, med_precision_target),
                    grid_step=0.001,
                    min_pos_calls=100,
                    use_wilson=True,
                    high_precision_floors=high_precision_floors or {},
                    med_precision_floors=med_precision_floors or {},
                )
                high_t = cutoffs_mge["high"]["t"]
                high_metrics = cutoffs_mge["high"]
                logger.info(
                    f"High threshold derived from MGE datasets only: {high_t:.3f} "
                    f"(was {cutoffs_all['high']['t']:.3f})"
                )

        self.confidence_thresholds = {
            "high": {"threshold": high_t, "metrics": high_metrics},
            "medium": {"threshold": cutoffs_all["medium"]["t"], "metrics": cutoffs_all["medium"]},
            "low": {"threshold": 0.0, "metrics": {}},
        }
        logger.info("Chosen thresholds:")
        logger.info(f"  HIGH={self.confidence_thresholds['high']['threshold']:.3f}")
        logger.info(f"  MEDIUM={self.confidence_thresholds['medium']['threshold']:.3f}")
        logger.info("  LOW=0.0")

        evaluation_results = self.evaluate_on_all_datasets(train_test_datasets_dict)
        summary_df = self.generate_performance_summary_table(evaluation_results)
        model_summary = self.get_model_summary()

        complete_results = {
            "training_results": training_results,
            "evaluation_results": evaluation_results,
            "summary_df": summary_df,
            "model_summary": model_summary,
            "trained_model": self,
            "datasets": train_test_datasets_dict,
        }

        logger.info("TRAINING COMPLETE")
        logger.info(f"Model trained on: {self.training_dataset}")
        logger.info(f"Evaluated on {len(train_test_datasets_dict)} datasets")
        logger.info(
            f"Confidence thresholds: "
            f"High>={self.confidence_thresholds['high']['threshold']:.3f}, "
            f"Medium>={self.confidence_thresholds['medium']['threshold']:.3f}, "
            f"Low>=0.000"
        )

        return complete_results


# ---------------------------------------------------------------------------
# Weight tuning helpers (unchanged except threading through precision_floors)
# ---------------------------------------------------------------------------

def compute_objective_mean_auprc(
    eval_res,
    ignore_substr=None,
    favor_dataset=None,
    favor_weight=2.0,
):
    ignore_substr_l = str(ignore_substr).lower().strip() if ignore_substr else None
    favor_dataset_l = str(favor_dataset).lower().strip() if favor_dataset else None

    per_ds = {}
    num = 0.0
    den = 0.0

    for k, v in eval_res.items():
        try:
            auprc = float(v["pr_curve"]["auprc"])
        except Exception:
            continue
        k_l = str(k).lower()
        if ignore_substr_l and ignore_substr_l in k_l:
            continue
        w = 1.0
        if favor_dataset_l and (k_l == favor_dataset_l):
            w = float(favor_weight)
        per_ds[k] = auprc
        num += w * auprc
        den += w

    if den <= 0:
        return -math.inf, math.nan, per_ds

    mean_w = float(num / den)
    return mean_w, mean_w, per_ds


def _parse_range(spec, name):
    if spec is None:
        return None
    s = str(spec).strip()
    if s == "":
        return None
    if (":" not in s) and ("," not in s):
        try:
            v = float(s)
        except Exception as e:
            raise ValueError(f"Invalid {name} range '{spec}': {e}")
        return [v]
    if ":" in s:
        parts = [p.strip() for p in s.split(":") if p.strip() != ""]
    else:
        parts = [p.strip() for p in s.split(",") if p.strip() != ""]
    if len(parts) not in (2, 3):
        raise ValueError(f"Invalid {name} range '{spec}'. Expected min:max[:step].")
    lo = float(parts[0])
    hi = float(parts[1])
    if hi < lo:
        raise ValueError(f"Invalid {name} range '{spec}': max < min")
    step = float(parts[2]) if len(parts) == 3 else 1.0
    if step <= 0:
        raise ValueError(f"Invalid {name} range '{spec}': step must be > 0")
    vals = []
    x = lo
    n_guard = 0
    while x <= hi + 1e-12:
        vals.append(float(x))
        x += step
        n_guard += 1
        if n_guard > 100000:
            raise ValueError(f"{name} range '{spec}' exploded (too many values).")
    seen = set()
    out = []
    for v in vals:
        key = round(v, 12)
        if key in seen:
            continue
        seen.add(key)
        out.append(v)
    return out


def _eval_one_weight_combo(
    datasets_dict, best_params, predictor_ctor_args,
    wm, wh, wp, group_col, early_stop, early_stopping_rounds,
    high_recall_target, high_precision_target,
    med_recall_target, med_precision_target,
    high_on_mge_datasets,
    ignore_substr_for_score, favor_dataset, favor_weight,
    high_precision_floors=None, med_precision_floors=None,
):
    predictor = ViralConfidencePredictor(**predictor_ctor_args)
    res = predictor.run_complete_analysis(
        train_test_datasets_dict=datasets_dict,
        best_params=best_params,
        w_neg_mge=float(wm), w_neg_host=float(wh), w_pos_virus=float(wp),
        group_col=group_col,
        use_early_stopping=early_stop,
        early_stopping_rounds=early_stopping_rounds,
        high_recall_target=high_recall_target,
        high_precision_target=high_precision_target,
        med_recall_target=med_recall_target,
        med_precision_target=med_precision_target,
        high_on_mge_datasets=high_on_mge_datasets,
        high_precision_floors=high_precision_floors or {},
        med_precision_floors=med_precision_floors or {},
    )
    eval_res = res["evaluation_results"]
    score, mean_all, per_ds = compute_objective_mean_auprc(
        eval_res,
        ignore_substr=ignore_substr_for_score,
        favor_dataset=favor_dataset,
        favor_weight=favor_weight,
    )
    return float(score), float(mean_all), per_ds


def tune_weights_iterative(
    datasets_dict, predictor_ctor_args,
    best_params=None, init_weights=(1.0, 1.0, 1.0),
    range_w_neg_mge=None, range_w_neg_host=None, range_w_pos_virus=None,
    max_iters=5, out_csv="weight_iterative_results.csv",
    group_col="Contig", early_stop=True, early_stopping_rounds=100,
    high_recall_target=0.50, high_precision_target=0.95,
    med_recall_target=0.75, med_precision_target=0.90,
    high_on_mge_datasets=False,
    ignore_substr_for_score=None, favor_dataset=None, favor_weight=2.0,
    parallel_jobs=1, total_threads=50,
    high_precision_floors=None, med_precision_floors=None,
):
    w_neg_mge, w_neg_host, w_pos_virus = map(float, init_weights)

    grids = {
        "w_neg_mge": [w_neg_mge] if range_w_neg_mge is None else list(range_w_neg_mge),
        "w_neg_host": [w_neg_host] if range_w_neg_host is None else list(range_w_neg_host),
        "w_pos_virus": [w_pos_virus] if range_w_pos_virus is None else list(range_w_pos_virus),
    }

    parallel_jobs = max(1, int(parallel_jobs))
    total_threads = max(1, int(total_threads))

    def _worker_args_for_njobs(n_jobs):
        n_jobs = max(1, int(n_jobs))
        per_worker_threads = max(1, total_threads // n_jobs)
        args2 = dict(predictor_ctor_args)
        args2["max_threads"] = per_worker_threads
        return args2, per_worker_threads

    def _eval_weights_once(wm, wh, wp, worker_ctor_args):
        return _eval_one_weight_combo(
            datasets_dict=datasets_dict, best_params=best_params,
            predictor_ctor_args=worker_ctor_args,
            wm=wm, wh=wh, wp=wp, group_col=group_col,
            early_stop=early_stop, early_stopping_rounds=early_stopping_rounds,
            high_recall_target=high_recall_target, high_precision_target=high_precision_target,
            med_recall_target=med_recall_target, med_precision_target=med_precision_target,
            high_on_mge_datasets=high_on_mge_datasets,
            ignore_substr_for_score=ignore_substr_for_score,
            favor_dataset=favor_dataset, favor_weight=favor_weight,
            high_precision_floors=high_precision_floors,
            med_precision_floors=med_precision_floors,
        )

    history_rows = []
    base_ctor_args, _ = _worker_args_for_njobs(1)
    best_score, best_mean_all, _ = _eval_weights_once(w_neg_mge, w_neg_host, w_pos_virus, base_ctor_args)
    logger.info(
        f"Iterative tune init: w_neg_mge={w_neg_mge} w_neg_host={w_neg_host} w_pos_virus={w_pos_virus} "
        f"-> score={best_score:.6f} mean_all={best_mean_all:.6f}"
    )

    for it in range(1, int(max_iters) + 1):
        improved_any = False
        for param_name in ("w_neg_mge", "w_neg_host", "w_pos_virus"):
            cand = grids[param_name]
            if cand is None or len(cand) == 0:
                continue
            local_best = {
                "score": best_score, "mean_all": best_mean_all,
                "val": {"w_neg_mge": w_neg_mge, "w_neg_host": w_neg_host, "w_pos_virus": w_pos_virus}[param_name],
            }
            tasks = []
            for v in cand:
                wm, wh, wp = w_neg_mge, w_neg_host, w_pos_virus
                if param_name == "w_neg_mge":
                    wm = float(v)
                elif param_name == "w_neg_host":
                    wh = float(v)
                else:
                    wp = float(v)
                tasks.append((float(v), float(wm), float(wh), float(wp)))

            n_jobs = min(parallel_jobs, len(tasks))
            worker_ctor_args, per_worker_threads = _worker_args_for_njobs(n_jobs)
            logger.info(
                f"[iter {it}] scanning {param_name} over {len(tasks)} candidates "
                f"with parallel_jobs={n_jobs}, threads_per_worker={per_worker_threads}"
            )

            if n_jobs == 1:
                results = []
                for v, wm, wh, wp in tasks:
                    score, mean_all, _ = _eval_weights_once(wm, wh, wp, worker_ctor_args)
                    results.append((v, wm, wh, wp, score, mean_all))
            else:
                results = Parallel(n_jobs=n_jobs, backend="loky")(
                    delayed(
                        lambda _v, _wm, _wh, _wp: (_v, _wm, _wh, _wp)
                        + _eval_weights_once(_wm, _wh, _wp, worker_ctor_args)[:2]
                    )(v, wm, wh, wp)
                    for (v, wm, wh, wp) in tasks
                )

            for (v, wm, wh, wp, score, mean_all) in results:
                history_rows.append({
                    "iter": it, "param": param_name, "value": float(v),
                    "w_neg_mge": float(wm), "w_neg_host": float(wh), "w_pos_virus": float(wp),
                    "score": float(score), "mean_all": float(mean_all),
                    "ignore_substr_for_score": str(ignore_substr_for_score) if ignore_substr_for_score else "",
                    "favor_dataset": str(favor_dataset) if favor_dataset else "",
                    "favor_weight": float(favor_weight) if favor_dataset else "",
                })
                logger.info(
                    f"[iter {it}] try {param_name}={float(v)} with "
                    f"(w_neg_mge={wm}, w_neg_host={wh}, w_pos_virus={wp}) "
                    f"-> score={score:.6f} mean_all={mean_all:.6f}"
                )
                if score > local_best["score"]:
                    local_best.update({"score": float(score), "mean_all": float(mean_all), "val": float(v)})

            if local_best["score"] > best_score:
                improved_any = True
                best_score = float(local_best["score"])
                best_mean_all = float(local_best["mean_all"])
                if param_name == "w_neg_mge":
                    w_neg_mge = float(local_best["val"])
                elif param_name == "w_neg_host":
                    w_neg_host = float(local_best["val"])
                else:
                    w_pos_virus = float(local_best["val"])
                logger.info(
                    f"[iter {it}] ACCEPT {param_name}={local_best['val']} -> "
                    f"best score={best_score:.6f} mean_all={best_mean_all:.6f} "
                    f"weights=(mge={w_neg_mge}, host={w_neg_host}, pos={w_pos_virus})"
                )

        if not improved_any:
            logger.info(f"Iterative tuning converged at iter {it}: no improvements found.")
            break

    keys = [
        "iter", "param", "value", "w_neg_mge", "w_neg_host", "w_pos_virus",
        "score", "mean_all", "ignore_substr_for_score", "favor_dataset", "favor_weight",
    ]
    with open(out_csv, "w", newline="") as cf:
        writer = csv.DictWriter(cf, fieldnames=keys)
        writer.writeheader()
        for r in history_rows:
            writer.writerow({k: r.get(k, "") for k in keys})

    best_obj = {
        "w_neg_mge": float(w_neg_mge), "w_neg_host": float(w_neg_host), "w_pos_virus": float(w_pos_virus),
        "best_score": float(best_score), "mean_all": float(best_mean_all),
        "ignore_substr_for_score": str(ignore_substr_for_score) if ignore_substr_for_score else None,
        "favor_dataset": str(favor_dataset) if favor_dataset else None,
        "favor_weight": float(favor_weight) if favor_dataset else None,
        "max_iters": int(max_iters), "parallel_jobs": int(parallel_jobs), "total_threads": int(total_threads),
    }
    logger.info(
        f"Iterative tuning done. Best weights: "
        f"w_neg_mge={best_obj['w_neg_mge']} w_neg_host={best_obj['w_neg_host']} w_pos_virus={best_obj['w_pos_virus']} "
        f"score={best_obj['best_score']:.6f} mean_all={best_obj['mean_all']:.6f}. "
        f"Saved history to {out_csv}"
    )
    return best_obj, history_rows


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_dfs(dataset_parent_path, input_format="auto"):
    logger.info("LOADING TRAIN AND TEST TABLES")
    input_format = str(input_format).lower().strip()
    if input_format not in ("auto", "tsv", "parquet"):
        raise ValueError(f"--input-format must be one of: auto, tsv, parquet (got '{input_format}')")

    patterns = []
    if input_format in ("auto", "tsv"):
        patterns.append(os.path.join(dataset_parent_path, "*.tsv"))
    if input_format in ("auto", "parquet"):
        patterns.append(os.path.join(dataset_parent_path, "*.parquet"))

    dataset_paths = []
    for pat in patterns:
        dataset_paths.extend(glob.glob(pat))
    dataset_paths = sorted(set(dataset_paths))

    if not dataset_paths:
        raise FileNotFoundError(
            f"No datasets found in '{dataset_parent_path}'. "
            f"Expected files matching: *.tsv and/or *.parquet (input_format={input_format})."
        )

    train_test_dfs = {}
    expected = set(FEATURE_COLS)

    for dataset_path in dataset_paths:
        base = os.path.basename(dataset_path)
        ext = os.path.splitext(base)[1].lower()
        dataset_name = base.replace("lgbm.", "").replace(ext, "")
        logger.info(f"Loading {dataset_name} from {dataset_path}")

        if ext == ".parquet":
            df = pl.read_parquet(dataset_path).to_pandas()
        elif ext == ".tsv":
            with open(dataset_path, "r", encoding="utf-8") as fh:
                header_line = fh.readline().rstrip("\n")
            file_cols = [c.strip() for c in header_line.split("\t")]
            logger.info(f"{dataset_name}: file has {len(file_cols)} columns")
            df = pl.read_csv(dataset_path, separator="\t", infer_schema_length=10000).to_pandas()
        else:
            raise ValueError(f"Unsupported file extension '{ext}' for {dataset_path}")

        missing = [c for c in expected if c not in df.columns and c not in ALLOWED_META_COLS]
        extra = [c for c in df.columns if c not in expected and c not in ALLOWED_META_COLS]
        if missing:
            logger.warning(f"{dataset_name}: missing expected feature columns: {missing}")
        if extra:
            logger.info(f"{dataset_name}: extra columns not expected: {extra}")

        for col in FEATURE_COLS:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors="coerce")

        train_test_dfs[dataset_name] = df

    return train_test_dfs


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(args):
    _set_global_thread_env(args.threads)

    train_test_dfs = load_dfs(args.dataset_path, input_format=args.input_format)

    best_params = None
    if args.best_params is not None and os.path.exists(args.best_params):
        with open(args.best_params, "r") as f:
            j = json.load(f)
            if isinstance(j, dict):
                if "best_config" in j:
                    best_params = j
                elif "model_config" in j:
                    best_params = {"best_config": j["model_config"], "best_score": j.get("best_score", j.get("cv_score"))}
                elif "best_params" in j:
                    best_params = {"best_config": j["best_params"], "best_score": j.get("best_score", j.get("cv_score"))}
                else:
                    best_params = j
        logger.info("Loaded best model params from JSON")

    # Parse precision floor specs
    high_precision_floors = _parse_precision_floors(args.high_precision_floor or [])
    med_precision_floors = _parse_precision_floors(args.med_precision_floor or [])
    if high_precision_floors:
        logger.info(f"Parsed high precision floors: {high_precision_floors}")
    if med_precision_floors:
        logger.info(f"Parsed medium precision floors: {med_precision_floors}")

    predictor_args = {
        "max_threads": args.threads,
        "reliability_n_bins": args.reliability_bins,
        "reliability_strategy": args.reliability_strategy,
        "bayes_candidates_in_parallel": args.bayes_candidates_in_parallel,
    }

    if args.tune_weights:
        range_mge = _parse_range(args.tune_w_neg_mge_range, "w_neg_mge")
        range_host = _parse_range(args.tune_w_neg_host_range, "w_neg_host")
        range_pos = _parse_range(args.tune_w_pos_virus_range, "w_pos_virus")

        try:
            init_w = [float(x) for x in args.weights.split(",")]
            if len(init_w) != 3:
                raise ValueError("Need exactly 3 values")
            init_w = tuple(init_w)
        except Exception:
            init_w = (1.0, 1.0, 1.0)

        out_csv = (
            os.path.join(args.tune_outdir, os.path.basename(args.output) + ".weights_iterative.csv")
            if args.tune_outdir else os.path.splitext(args.output)[0] + ".weights_iterative.csv"
        )

        best_obj, history_rows = tune_weights_iterative(
            datasets_dict=train_test_dfs,
            predictor_ctor_args=predictor_args,
            best_params=best_params,
            init_weights=init_w,
            range_w_neg_mge=range_mge,
            range_w_neg_host=range_host,
            range_w_pos_virus=range_pos,
            max_iters=int(args.tune_max_iters),
            out_csv=out_csv,
            group_col=args.group_col,
            early_stop=args.early_stop,
            early_stopping_rounds=int(args.early_stopping_rounds),
            high_recall_target=float(args.high_recall_target),
            high_precision_target=float(args.high_precision_target),
            med_recall_target=float(args.med_recall_target),
            med_precision_target=float(args.med_precision_target),
            high_on_mge_datasets=bool(args.high_on_mge_datasets),
            ignore_substr_for_score=args.tune_ignore_substr_for_score,
            favor_dataset=args.tune_favor_dataset,
            favor_weight=float(args.tune_favor_weight),
            parallel_jobs=int(args.tune_parallel_jobs),
            total_threads=int(args.threads),
            high_precision_floors=high_precision_floors,
            med_precision_floors=med_precision_floors,
        )

        if args.tune_outdir:
            os.makedirs(args.tune_outdir, exist_ok=True)
            hist_pickle = os.path.join(args.tune_outdir, "weight_iter_history.pkl")
            with open(hist_pickle, "wb") as fh:
                pickle.dump(history_rows, fh)
            best_json = os.path.join(args.tune_outdir, "best_weights.json")
            with open(best_json, "w") as fh:
                json.dump(best_obj, fh, indent=2)
            logger.info(f"Saved iterative tune history to {hist_pickle}")
            logger.info(f"Saved best weights to {best_json}")

        chosen = (best_obj["w_neg_mge"], best_obj["w_neg_host"], best_obj["w_pos_virus"])
        logger.info(f"Retraining final model with tuned weights: mge={chosen[0]} host={chosen[1]} pos={chosen[2]}")

        predictor = ViralConfidencePredictor(**predictor_args)
        results = predictor.run_complete_analysis(
            train_test_datasets_dict=train_test_dfs,
            best_params=best_params,
            w_neg_mge=float(chosen[0]), w_neg_host=float(chosen[1]), w_pos_virus=float(chosen[2]),
            group_col=args.group_col,
            use_early_stopping=args.early_stop,
            early_stopping_rounds=args.early_stopping_rounds,
            high_recall_target=args.high_recall_target,
            high_precision_target=args.high_precision_target,
            med_recall_target=args.med_recall_target,
            med_precision_target=args.med_precision_target,
            high_on_mge_datasets=args.high_on_mge_datasets,
            high_precision_floors=high_precision_floors,
            med_precision_floors=med_precision_floors,
        )

    else:
        try:
            w_vals = [float(x) for x in args.weights.split(",")]
            w_neg_mge, w_neg_host, w_pos_virus = w_vals[0], w_vals[1], w_vals[2]
        except Exception:
            logger.warning("Could not parse --weights; falling back to 1.0,1.0,1.0")
            w_neg_mge, w_neg_host, w_pos_virus = 1.0, 1.0, 1.0

        predictor = ViralConfidencePredictor(**predictor_args)
        results = predictor.run_complete_analysis(
            train_test_datasets_dict=train_test_dfs,
            best_params=best_params,
            w_neg_mge=w_neg_mge, w_neg_host=w_neg_host, w_pos_virus=w_pos_virus,
            group_col=args.group_col,
            use_early_stopping=args.early_stop,
            early_stopping_rounds=args.early_stopping_rounds,
            high_recall_target=args.high_recall_target,
            high_precision_target=args.high_precision_target,
            med_recall_target=args.med_recall_target,
            med_precision_target=args.med_precision_target,
            high_on_mge_datasets=args.high_on_mge_datasets,
            high_precision_floors=high_precision_floors,
            med_precision_floors=med_precision_floors,
        )

    logger.info(f"Writing results to {args.output}")
    with open(args.output, "wb") as f:
        pickle.dump(results, f)

    best_params_out = os.path.splitext(args.output)[0] + ".best_params.json"
    model_cfg = results["training_results"]["model_config"]
    cv_score = results["training_results"].get("cv_score", results["training_results"].get("best_score"))
    out_obj = {"best_config": model_cfg, "best_score": cv_score}
    with open(best_params_out, "w") as f:
        json.dump(out_obj, f, indent=2)
    logger.info(f"Wrote best params to {best_params_out}")
    logger.info("ALL DONE")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train and evaluate Viral Confidence LGBM model.")
    parser.add_argument("-d", "--dataset-path", required=True,
        help="Path to directory containing lgbm.*.tsv / lgbm.*.parquet datasets.")
    parser.add_argument("--input-format", choices=["auto", "tsv", "parquet"], default="auto",
        help="Input dataset format in --dataset-path. Expects lgbm.*.tsv and/or lgbm.*.parquet (default: auto).")
    parser.add_argument("-t", "--threads", type=int, default=50,
        help="Total thread budget across all workers (LightGBM + joblib).")
    parser.add_argument("-b", "--best-params", default=None,
        help="Path to JSON file with best_params (optional).")
    parser.add_argument("-w", "--weights", default="1.0,1.0,1.0",
        help="Comma-separated sample weights: w_neg_mge,w_neg_host,w_pos_virus (default: 1.0,1.0,1.0).")
    parser.add_argument("--group-col", choices=["Contig", "Genome"], default="Contig",
        help="Column to use for grouping in GroupShuffleSplit and CV (default: Contig).")
    parser.add_argument("--early-stopping-rounds", type=int, default=200,
        help="Rounds for LightGBM early stopping when eval_set is used (default: 200).")
    parser.add_argument("--early-stop", dest="early_stop", action="store_true",
        help="Enable weighted early stopping with eval set (default: enabled).")
    parser.add_argument("--no-early-stop", dest="early_stop", action="store_false",
        help="Disable early stopping.")
    parser.set_defaults(early_stop=True)

    # Global precision / recall targets
    parser.add_argument("--high-recall-target", type=float, default=0.50,
        help="Recall target for HIGH threshold selection (default: 0.50).")
    parser.add_argument("--high-precision-target", type=float, default=0.95,
        help="Precision target for HIGH threshold selection (default: 0.95).")
    parser.add_argument("--med-recall-target", type=float, default=0.75,
        help="Recall target for MEDIUM threshold selection (default: 0.75).")
    parser.add_argument("--med-precision-target", type=float, default=0.90,
        help="Precision target for MEDIUM threshold selection (default: 0.90).")
    parser.add_argument("--high-on-mge-datasets", action="store_true",
        help="If set, derive HIGH threshold using only datasets whose key contains 'mge'.")

    # Per-dataset precision floor constraints
    parser.add_argument(
        "--high-precision-floor",
        action="append",
        metavar="DATASET:MIN_PREC",
        dest="high_precision_floor",
        default=None,
        help=(
            "Hard precision floor for high-confidence threshold selection. "
            "Format: 'dataset_name:min_precision' (e.g. 'test_provirus:0.90'). "
            "Repeatable. A threshold is ineligible if any floor is violated."
        ),
    )
    parser.add_argument(
        "--med-precision-floor",
        action="append",
        metavar="DATASET:MIN_PREC",
        dest="med_precision_floor",
        default=None,
        help=(
            "Hard precision floor for medium-confidence threshold selection. "
            "Format: 'dataset_name:min_precision' (e.g. 'test_provirus:0.85'). "
            "Repeatable."
        ),
    )

    # Weight tuning
    parser.add_argument("--tune-weights", action="store_true",
        help="If set, iteratively tune w_neg_mge, w_neg_host, and w_pos_virus over provided ranges to optimize mean AUPRC across test datasets.")
    parser.add_argument("--tune-w-neg-mge-range", default=None,
        help="Range for w_neg_mge as 'min:max:step' (or 'min,max,step'). Example: '0.5:10:0.5'.")
    parser.add_argument("--tune-w-neg-host-range", default=None,
        help="Range for w_neg_host as 'min:max:step' (or 'min,max,step'). Example: '0.5:10:0.5'.")
    parser.add_argument("--tune-w-pos-virus-range", default=None,
        help="Range for w_pos_virus as 'min:max:step' (or 'min,max,step'). Example: '0.5:10:0.5'.")
    parser.add_argument("--tune-max-iters", type=int, default=5,
        help="Max coordinate-descent passes over weights (default: 5).")
    parser.add_argument("--tune-outdir", default=None,
        help="If tuning used, write per-run pickles & summary to this directory.")
    parser.add_argument("--tune-parallel-jobs", type=int, default=1,
        help="Parallel jobs for candidate evaluation during weight tuning. Each job gets threads//jobs threads (default: 1).")
    parser.add_argument("--tune-ignore-substr-for-score", default=None,
        help="If set, ignore datasets whose name contains this substring when computing the tuning objective (e.g. 'mge').")
    parser.add_argument("--tune-favor-dataset", default=None,
        help="If set, favor this dataset key (exact match) in the tuning objective mean AUPRC by multiplying its weight.")
    parser.add_argument("--tune-favor-weight", type=float, default=2.0,
        help="Weight multiplier for --tune-favor-dataset (default: 2.0).")

    parser.add_argument("--bayes-candidates-in-parallel", type=int, default=1,
        help="How many BayesSearch candidates to evaluate concurrently. Keep small; default 1 (fold-parallel).")
    parser.add_argument("--reliability-bins", type=int, default=15,
        help="Number of bins for reliability (calibration) computation (default: 15).")
    parser.add_argument("--reliability-strategy", choices=["uniform", "quantile"], default="uniform",
        help="Binning strategy for reliability: 'uniform' equal-width, 'quantile' equal-count (default: uniform).")
    parser.add_argument("-o", "--output", required=True,
        help="Path to save results pickle file.")

    args = parser.parse_args()
    main(args)