"""Resolve the annotate and de-novo database directories from a user-provided
--db-dir.

The annotate and de-novo databases are distributed separately, each extracting
to its own named directory. A --db-dir may therefore point either directly at a
database directory or at a parent that contains it as an immediate subdirectory
(e.g. when both databases are downloaded side by side). These helpers locate the
correct directory by content marker so both standalone modules and end-to-end
work with either layout.
"""

from __future__ import annotations

import glob
import os

# Files that uniquely identify each database directory by content.
_ANNOTATE_MARKERS = ("KEGG_cutoffs.tsv", "FOAM_cutoffs.tsv", "KEGG.h3f", "KEGG.hmm")


def _immediate_subdirs(d):
    try:
        return [os.path.join(d, s) for s in sorted(os.listdir(d)) if os.path.isdir(os.path.join(d, s))]
    except OSError:
        return []


def _has_any_file(d, names):
    return any(os.path.isfile(os.path.join(d, n)) for n in names)


def _has_glob(d, pattern):
    return bool(glob.glob(os.path.join(d, pattern)))


def resolve_annotate_db(db_dir):
    """Return the directory containing the annotate HMM/cutoff files.

    Checks db_dir itself, then its immediate subdirectories. Falls back to
    db_dir unchanged so downstream code raises a clear file-not-found error.
    """
    db_dir = os.path.abspath(db_dir)
    if _has_any_file(db_dir, _ANNOTATE_MARKERS):
        return db_dir
    for sub in _immediate_subdirs(db_dir):
        if _has_any_file(sub, _ANNOTATE_MARKERS):
            return sub
    return db_dir


def resolve_denovo_db(db_dir):
    """Return the directory containing the de-novo model checkpoint.

    Checks db_dir, a legacy 'checkamg_pst' subdir, then immediate
    subdirectories. Falls back to db_dir unchanged.
    """
    db_dir = os.path.abspath(db_dir)
    candidates = [db_dir, os.path.join(db_dir, "checkamg_pst")] + _immediate_subdirs(db_dir)
    seen = set()
    for d in candidates:
        if d in seen:
            continue
        seen.add(d)
        if os.path.isdir(d) and _has_glob(d, "*.ckpt"):
            return d
    return db_dir
