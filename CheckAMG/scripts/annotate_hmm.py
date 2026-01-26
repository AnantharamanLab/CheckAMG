#!/usr/bin/env python3

import os
import sys
import resource
import logging
os.environ["POLARS_MAX_THREADS"] = str(snakemake.threads)
import polars as pl
from pathlib import Path
from pyhmmer import easel, plan7, hmmer
import json
import math
import uuid
from datetime import datetime
from pyfastatools import Parser, write_fasta
from tqdm import tqdm

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
    format="%(asctime)s | %(levelname)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.FileHandler(log_file, mode="a"),
        logging.StreamHandler(sys.stdout),
    ],
)
logger = logging.getLogger()

print("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n                Step 5/11: Assign functions to proteins                 \n========================================================================\n")

def atomic_write_text(out_path: Path, write_fn, mode="w", encoding="utf-8"):
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with open(tmp_path, mode, encoding=encoding) as f:
            write_fn(f)
            f.flush()
            os.fsync(f.fileno())
        os.replace(tmp_path, out_path)
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except Exception:
                pass

def write_done_marker(done_path: Path, meta: dict):
    def _w(f):
        f.write(json.dumps(meta, sort_keys=True) + "\n")

    atomic_write_text(done_path, _w)

def done_path_for(path: Path) -> Path:
    return Path(str(path) + ".done")

def is_complete(path: Path) -> bool:
    dp = done_path_for(path)
    if not path.exists():
        return False
    if not dp.exists():
        return False
    try:
        if path.stat().st_size == 0:
            return False
    except Exception:
        return False
    return True

def remove_incomplete_outputs(paths):
    for p in paths:
        p = Path(p)
        dp = done_path_for(p)
        tp = p.with_suffix(p.suffix + ".tmp")
        if tp.exists():
            try:
                tp.unlink()
            except Exception:
                pass
        if p.exists() and not dp.exists():
            try:
                p.unlink()
            except Exception:
                pass
        if dp.exists() and not p.exists():
            try:
                dp.unlink()
            except Exception:
                pass

def extract_query_info(hits, db_path):
    s = str(db_path)
    if "Pfam" in s or "pfam" in s:
        hmm_id = hits.query.accession.decode()
    elif "FOAM" in s or "foam" in s:
        hmm_id = hits.query.accession.decode()
    elif "eggNOG" in s or "eggnog" in s:
        hmm_id = hits.query.name.decode().split(".")[0]
    else:
        query_name = hits.query.name.decode()
        if ".wlink.txt.mafft" in query_name:
            hmm_id = query_name.split(".")[1]
        else:
            hmm_id = (
                query_name.replace("_alignment", "")
                .replace(".mafft", "")
                .replace(".txt", "")
                .replace(".hmm", "")
                .replace("_protein.alignment", "")
            )
    return hmm_id

def aggregate_sequences(protein_dir):
    all_sequences = []
    protein_dir = Path(protein_dir)
    for fasta_file in protein_dir.rglob("*"):
        if fasta_file.suffix.lower() in (".faa", ".fasta"):
            all_sequences.extend(Parser(str(fasta_file)).all())
    return all_sequences

def split_aggregated_sequences(all_sequences, chunk_size):
    for i in range(0, len(all_sequences), chunk_size):
        yield all_sequences[i:i + chunk_size]

def determine_chunk_size(n_sequences, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8):
    total_bytes = n_sequences * est_bytes_per_seq
    allowed_bytes = max_chunk_fraction * mem_limit * (1024**3)
    n_chunks = max(1, math.ceil(total_bytes / allowed_bytes))
    return math.ceil(n_sequences / n_chunks)

def find_resumable_tmp_dir(wdir: Path) -> Path | None:
    wdir = Path(wdir)
    candidates = sorted(
        [p for p in wdir.glob("hmmsearch_tmp_*") if p.is_dir()],
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    for d in candidates:
        if (d / "manifest.json").exists():
            return d
        if any(d.glob("seqbatch_*.faa")):
            return d
    return None

def load_manifest(tmp_dir: Path) -> dict | None:
    mpath = tmp_dir / "manifest.json"
    if not mpath.exists():
        return None
    try:
        return json.loads(mpath.read_text())
    except Exception:
        return None

def write_manifest(tmp_dir: Path, manifest: dict):
    mpath = tmp_dir / "manifest.json"

    def _w(f):
        f.write(json.dumps(manifest, sort_keys=True, indent=2) + "\n")

    atomic_write_text(mpath, _w)

def assign_db(db_path):
    s = str(db_path)
    if "KEGG" in s or "kegg" in s or "kofam" in s:
        return "KEGG"
    elif "FOAM" in s or "foam" in s:
        return "FOAM"
    elif "Pfam" in s or "pfam" in s:
        return "Pfam"
    elif "dbcan" in s or "dbCAN" in s or "dbCan" in s:
        return "dbCAN"
    elif "METABOLIC_custom" in s or "metabolic_custom" in s:
        return "METABOLIC"
    elif "CAMPER" in s or "camper" in s:
        return "CAMPER"
    elif "VOG" in s or "vog" in s:
        return "VOG"
    elif "eggNOG" in s or "eggnog" in s:
        return "eggNOG"
    elif "PHROG" in s or "phrog" in s:
        return "PHROG"
    elif "user_custom" in s:
        return "user_custom"
    else:
        return None
    
def _nan_to_none(v):
    if v is None:
        return None
    if isinstance(v, float) and math.isnan(v):
        return None
    return float(v)


# Load KEGG thresholds
KEGG_THRESHOLDS = {}
KEGG_THRESHOLDS_PATH = snakemake.params.kegg_cutoff_file
if Path(KEGG_THRESHOLDS_PATH).exists():
    kegg_df = pl.read_csv(
        KEGG_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "threshold": pl.Float64},
        separator="\t",
    )
    KEGG_THRESHOLDS = dict(zip(kegg_df["id"].to_list(), kegg_df["threshold"].to_list()))
else:
    logger.warning(f"KEGG thresholds file not found at {KEGG_THRESHOLDS_PATH}. KEGG thresholds will not be used to filter HMMsearch results!")

# Load FOAM thresholds
FOAM_THRESHOLDS = {}
FOAM_THRESHOLDS_PATH = snakemake.params.foam_cutoff_file
if Path(FOAM_THRESHOLDS_PATH).exists():
    foam_df = pl.read_csv(
        FOAM_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "cutoff_full": pl.Float64, "cutoff_domain": pl.Float64},
        separator="\t",
    )
    ids = foam_df["id"].to_list()
    fulls = foam_df["cutoff_full"].to_list()
    doms = foam_df["cutoff_domain"].to_list()
    FOAM_THRESHOLDS = {i: (_nan_to_none(f), _nan_to_none(d)) for i, f, d in zip(ids, fulls, doms)}
else:
    logger.warning(f"FOAM thresholds file not found at {FOAM_THRESHOLDS_PATH}. FOAM thresholds will not be used to filter HMMsearch results!")

# Load CAMPER thresholds
CAMPER_THRESHOLDS = {}
CAMPER_THRESHOLDS_PATH = snakemake.params.camper_cutoff_file
if Path(CAMPER_THRESHOLDS_PATH).exists():
    camper_df = pl.read_csv(
        CAMPER_THRESHOLDS_PATH,
        schema_overrides={"id": pl.Utf8, "cutoff_full": pl.Float64, "cutoff_domain": pl.Float64},
        separator="\t",
    )
    ids = camper_df["id"].to_list()
    fulls = camper_df["cutoff_full"].to_list()
    doms = camper_df["cutoff_domain"].to_list()
    CAMPER_THRESHOLDS = {i: (_nan_to_none(f), _nan_to_none(d)) for i, f, d in zip(ids, fulls, doms)}
else:
    logger.warning(f"CAMPER thresholds file not found at {CAMPER_THRESHOLDS_PATH}. CAMPER thresholds will not be used to filter HMMsearch results!")


def get_kegg_threshold(hmm_id):
    return KEGG_THRESHOLDS.get(hmm_id, None)

def get_foam_threshold(hmm_id):
    return FOAM_THRESHOLDS.get(hmm_id, (None, None))

def get_camper_threshold(hmm_id):
    return CAMPER_THRESHOLDS.get(hmm_id, (None, None))

def standardize_and_filter_hmm_results(tsv_path, hmm_path, out_full_path, out_best_path):
    db = assign_db(hmm_path)

    best_by_seq = {}

    def _is_alignment_row(alignment_type):
        return alignment_type in ("full", "domain")

    def _better(a, b):
        return a < b

    tsv_path = Path(tsv_path)
    out_full_path = Path(out_full_path)
    out_best_path = Path(out_best_path)

    remove_incomplete_outputs([out_full_path, out_best_path])

    def _write_full(f_out):
        f_out.write(
            "hmm_id\tdb\tsequence\tevalue\tscore\talignment_type\tcoverage_sequence\tcoverage_hmm\tkeep\tnote\n"
        )
        with open(tsv_path, "r", encoding="utf-8", errors="replace") as f_in:
            for line in f_in:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 9:
                    continue

                sequence, hmm_id, evalue, score, alignment_type, seq_cov, hmm_cov, keep, note = parts

                try:
                    evalue_f = float(evalue)
                    score_f = float(score)
                    seq_cov_f = float(seq_cov)
                    hmm_cov_f = float(hmm_cov)
                except ValueError:
                    continue

                keep_bool = str(keep).strip().lower() == "true"
                keep_txt = "true" if keep_bool else "false"

                f_out.write(
                    f"{hmm_id}\t{db}\t{sequence}\t{evalue_f}\t{score_f:.6f}\t{alignment_type}\t"
                    f"{seq_cov_f:.3f}\t{hmm_cov_f:.3f}\t{keep_txt}\t{note}\n"
                )

                if not keep_bool:
                    continue
                if not _is_alignment_row(alignment_type):
                    continue

                cand_key = (evalue_f, -score_f, -hmm_cov_f, -seq_cov_f, hmm_id, alignment_type, keep_bool, note)
                prev = best_by_seq.get(sequence)
                if prev is None or _better(cand_key, prev[0]):
                    best_by_seq[sequence] = (
                        cand_key,
                        (hmm_id, db, sequence, evalue_f, score_f, alignment_type, seq_cov_f, hmm_cov_f),
                    )

    atomic_write_text(out_full_path, _write_full)
    write_done_marker(done_path_for(out_full_path), {"tsv": str(out_full_path), "db": db, "status": "ok"})

    def _write_best(f_out):
        f_out.write("hmm_id\tdb\tsequence\tevalue\tscore\talignment_type\tcoverage_sequence\tcoverage_hmm\n")
        for _, row in best_by_seq.items():
            hmm_id, db, seq, evalue_f, score_f, alignment_type, seq_cov_f, hmm_cov_f = row[1]
            f_out.write(
                f"{hmm_id}\t{db}\t{seq}\t{evalue_f}\t{score_f:.6f}\t{alignment_type}\t{seq_cov_f:.3f}\t{hmm_cov_f:.3f}\n"
            )

    atomic_write_text(out_best_path, _write_best)
    write_done_marker(done_path_for(out_best_path), {"tsv": str(out_best_path), "db": db, "status": "ok"})

def hmmsearch_serial(
    outfile: Path,
    batch_fasta: str,
    db_path: str,
    seq_lengths: dict,
    min_coverage: float,
    min_score: float,
    min_bitscore_fraction: float,
    evalue: float,
    cpus: int,
):
    outfile = Path(outfile)
    remove_incomplete_outputs([outfile])

    alphabet = easel.Alphabet.amino()
    hmm_list = list(plan7.HMMFile(db_path))
    db = assign_db(db_path)

    def _merge_len(intervals):
        if not intervals:
            return 0
        intervals = sorted(intervals)
        total = 0
        cur_s, cur_e = intervals[0]
        for s, e in intervals[1:]:
            if s <= cur_e + 1:
                if e > cur_e:
                    cur_e = e
            else:
                total += cur_e - cur_s + 1
                cur_s, cur_e = s, e
        total += cur_e - cur_s + 1
        return total

    tmp_out = outfile.with_suffix(outfile.suffix + ".tmp")

    n_rows = 0
    with open(tmp_out, "w", encoding="utf-8") as out, easel.SequenceFile(
        batch_fasta, digital=True, alphabet=alphabet
    ) as seqs:
        for hits in hmmer.hmmsearch(queries=hmm_list, sequences=seqs, E=0.01, cpus=cpus):
            hmm = hits.query
            hmm_id = extract_query_info(hits, db_path)

            for hit in hits:
                hit_name = hit.name.decode()

                domains = list(hit.domains.reported)
                if not domains:
                    continue

                target_intervals = []
                hmm_intervals = []
                hmm_length = None
                for dom in domains:
                    a = dom.alignment
                    target_intervals.append((a.target_from, a.target_to))
                    hmm_intervals.append((a.hmm_from, a.hmm_to))
                    if hmm_length is None:
                        hmm_length = a.hmm_length

                target_aln_len = _merge_len(target_intervals)
                hmm_aln_len = _merge_len(hmm_intervals)

                seq_len = seq_lengths.get(hit_name, 0)
                seq_cov_full = (target_aln_len / seq_len) if seq_len else 0.0
                hmm_cov_full = (hmm_aln_len / hmm_length) if (hmm_length and hmm_length > 0) else 0.0

                # Default to using sequence-level metrics (one output row per hit)
                alignment_type = "full"
                score = hit.score
                reported_evalue = hit.evalue
                seq_coverage = seq_cov_full
                profile_coverage = hmm_cov_full

                # Default keep/note
                keep = True
                note = ""

                # Pfam: apply GA cutoffs where available, otherwise default
                if db == "Pfam":
                    note += "Pfam;"
                    if profile_coverage < min_coverage:
                        note += "coverage_below_minimum;"
                        keep = False
                    else:
                        note += "coverage_above_minimum;"
                    if hmm.cutoffs.gathering is not None:
                        note += "has_valid_GA;"
                        if score < hmm.cutoffs.gathering1:
                            note += "score_below_GA;"
                            keep = False
                        else:
                            note += "score_above_GA;"
                    else:
                        note += "no_GA_found;"
                        if score < min_score:
                            note += "score_fails_default_filter;"
                            keep = False
                        else:
                            note += "score_passes_default_filter;"

                # KEGG: check threshold + heuristic first, then default
                elif db == "KEGG":
                    note += "KEGG;"
                    if profile_coverage < min_coverage:
                        note += "coverage_below_minimum;"
                        keep = False
                    else:
                        note += "coverage_above_minimum;"
                    kegg_thresh = get_kegg_threshold(hmm_id)
                    if kegg_thresh is not None:
                        note += "has_valid_threshold;"
                        if score < kegg_thresh:
                            note += "score_below_threshold;"
                            if reported_evalue > evalue:
                                note += "evalue_fails_heuristic;"
                                keep = False
                            else:
                                note += "evalue_passes_heuristic;"
                            if score < min_bitscore_fraction * kegg_thresh:
                                note += "score_fails_heuristic;"
                                keep = False
                            else:
                                note += "score_passes_heuristic;"
                        else:
                            note += "score_above_threshold;"
                    else:
                        note += "no_threshold_found;"
                        if score < min_score:
                            note += "score_fails_default_filter;"
                            keep = False
                        else:
                            note += "score_passes_default_filter;"

                # FOAM: check full threshold first, then domain threshold, then default
                elif db == "FOAM":
                    note += "FOAM;"
                    foam_full, foam_dom = get_foam_threshold(hmm_id)

                    if foam_full is not None:
                        note += "has_valid_full_threshold;"
                        if profile_coverage < min_coverage:
                            note += "coverage_below_minimum;"
                            keep = False
                        else:
                            note += "coverage_above_minimum;"
                        if score < foam_full:
                            note += "score_below_full_threshold;"
                            if reported_evalue > evalue:
                                note += "evalue_fails_heuristic;"
                                keep = False
                            else:
                                note += "evalue_passes_heuristic;"
                            if score < min_bitscore_fraction * foam_full:
                                note += "score_fails_heuristic;"
                                keep = False
                            else:
                                note += "score_passes_heuristic;"
                        else:
                            note += "score_above_full_threshold;"

                    elif foam_dom is not None:
                        note += "has_valid_domain_threshold;"
                        note += "summary_row;domain_rows_emitted;"
                        any_domain_keep = False

                        for dom in domains:
                            a = dom.alignment
                            dom_target_len = a.target_to - a.target_from + 1
                            dom_hmm_len = a.hmm_to - a.hmm_from + 1

                            dom_seq_cov = (dom_target_len / seq_len) if seq_len else 0.0
                            dom_hmm_cov = (dom_hmm_len / a.hmm_length) if (a.hmm_length and a.hmm_length > 0) else 0.0

                            dom_score = dom.score
                            dom_evalue = getattr(dom, "i_evalue", None)
                            if dom_evalue is None:
                                dom_evalue = hit.evalue

                            dom_keep = True
                            dom_note = "FOAM;has_valid_domain_threshold;domain_row;"

                            if dom_hmm_cov < min_coverage:
                                dom_note += "coverage_below_minimum;"
                                dom_keep = False
                            else:
                                dom_note += "coverage_above_minimum;"

                            if dom_score < foam_dom:
                                dom_note += "score_below_domain_threshold;"
                                if dom_evalue > evalue:
                                    dom_note += "evalue_fails_heuristic;"
                                    dom_keep = False
                                else:
                                    dom_note += "evalue_passes_heuristic;"
                                if dom_score < min_bitscore_fraction * foam_dom:
                                    dom_note += "score_fails_heuristic;"
                                    dom_keep = False
                                else:
                                    dom_note += "score_passes_heuristic;"
                            else:
                                dom_note += "score_above_domain_threshold;"

                            out.write(
                                f"{hit_name}\t{hmm_id}\t{dom_evalue:.1E}\t{dom_score:.6f}\t"
                                f"domain\t{dom_seq_cov}\t{dom_hmm_cov}\t{dom_keep}\t{dom_note}\n"
                            )
                            if dom_keep:
                                any_domain_keep = True

                        keep = any_domain_keep
                        if keep:
                            note += "any_domain_kept;"
                        else:
                            note += "no_domains_kept;"

                        alignment_type = "summary"

                    else:
                        note += "no_threshold_found;"
                        if profile_coverage < min_coverage:
                            note += "coverage_below_minimum;"
                            keep = False
                        else:
                            note += "coverage_above_minimum;"
                        if score < min_score:
                            note += "score_fails_default_filter;"
                            keep = False
                        else:
                            note += "score_passes_default_filter;"

                # CAMPER: check full threshold first, then domain threshold, then default
                elif db == "CAMPER":
                    note += "CAMPER;"
                    camper_full, camper_dom = get_camper_threshold(hmm_id)

                    if camper_full is not None:
                        note += "has_valid_full_threshold;"
                        if profile_coverage < min_coverage:
                            note += "coverage_below_minimum;"
                            keep = False
                        else:
                            note += "coverage_above_minimum;"
                        if score < camper_full:
                            note += "score_below_full_threshold;"
                            if reported_evalue > evalue:
                                note += "evalue_fails_heuristic;"
                                keep = False
                            else:
                                note += "evalue_passes_heuristic;"
                            if score < min_bitscore_fraction * camper_full:
                                note += "score_fails_heuristic;"
                                keep = False
                            else:
                                note += "score_passes_heuristic;"
                        else:
                            note += "score_above_full_threshold;"

                    elif camper_dom is not None:
                        note += "has_valid_domain_threshold;"
                        note += "summary_row;domain_rows_emitted;"
                        any_domain_keep = False

                        for dom in domains:
                            a = dom.alignment
                            dom_target_len = a.target_to - a.target_from + 1
                            dom_hmm_len = a.hmm_to - a.hmm_from + 1

                            dom_seq_cov = (dom_target_len / seq_len) if seq_len else 0.0
                            dom_hmm_cov = (dom_hmm_len / a.hmm_length) if (a.hmm_length and a.hmm_length > 0) else 0.0

                            dom_score = dom.score
                            dom_evalue = getattr(dom, "i_evalue", None)
                            if dom_evalue is None:
                                dom_evalue = hit.evalue

                            dom_keep = True
                            dom_note = "CAMPER;has_valid_domain_threshold;domain_row;"

                            if dom_hmm_cov < min_coverage:
                                dom_note += "coverage_below_minimum;"
                                dom_keep = False
                            else:
                                dom_note += "coverage_above_minimum;"

                            if dom_score < camper_dom:
                                dom_note += "score_below_domain_threshold;"
                                if dom_evalue > evalue:
                                    dom_note += "evalue_fails_heuristic;"
                                    dom_keep = False
                                else:
                                    dom_note += "evalue_passes_heuristic;"
                                if dom_score < min_bitscore_fraction * camper_dom:
                                    dom_note += "score_fails_heuristic;"
                                    dom_keep = False
                                else:
                                    dom_note += "score_passes_heuristic;"
                            else:
                                dom_note += "score_above_domain_threshold;"

                            out.write(
                                f"{hit_name}\t{hmm_id}\t{dom_evalue:.1E}\t{dom_score:.6f}\t"
                                f"domain\t{dom_seq_cov}\t{dom_hmm_cov}\t{dom_keep}\t{dom_note}\n"
                            )
                            if dom_keep:
                                any_domain_keep = True

                        keep = any_domain_keep
                        if keep:
                            note += "any_domain_kept;"
                        else:
                            note += "no_domains_kept;"

                        alignment_type = "summary"

                    else:
                        note += "no_threshold_found;"
                        if profile_coverage < min_coverage:
                            note += "coverage_below_minimum;"
                            keep = False
                        else:
                            note += "coverage_above_minimum;"
                        if score < min_score:
                            note += "score_fails_default_filter;"
                            keep = False
                        else:
                            note += "score_passes_default_filter;"

                # METABOLIC GA where available, otherwise default
                elif db == "METABOLIC":
                    note += "METABOLIC;"
                    if profile_coverage < min_coverage:
                        note += "coverage_below_minimum;"
                        keep = False
                    else:
                        note += "coverage_above_minimum;"
                    if hmm.cutoffs.gathering is not None:
                        note += "has_valid_GA;"
                        if score < hmm.cutoffs.gathering1:
                            note += "score_below_GA;"
                            keep = False
                        else:
                            note += "score_above_GA;"
                    else:
                        note += "no_GA_found;"
                        if score < min_score:
                            note += "score_fails_default_filter;"
                            keep = False
                        else:
                            note += "score_passes_default_filter;"

                # Default fallback for other databases (like PHROG and dbCAN)
                else:
                    note += f"{db};"
                    if profile_coverage < min_coverage:
                        note += "coverage_below_minimum;"
                        keep = False
                    else:
                        note += "coverage_above_minimum;"
                    if score < min_score:
                        note += "score_fails_default_filter;"
                        keep = False
                    else:
                        note += "score_passes_default_filter;"

                out.write(
                    f"{hit_name}\t{hmm_id}\t{reported_evalue:.1E}\t{score:.6f}\t{alignment_type}\t"
                    f"{seq_coverage}\t{profile_coverage}\t{keep}\t{note}\n"
                )
                n_rows += 1

        out.flush()
        os.fsync(out.fileno())

    os.replace(tmp_out, outfile)
    write_done_marker(done_path_for(outfile), {"tsv": str(outfile), "db": db, "rows": n_rows, "status": "ok"})
    return str(outfile)

def main():
    protein_dir = snakemake.params.protein_dir
    wdir = Path(snakemake.params.wdir)
    hmm_vscores = snakemake.params.hmm_vscores
    cov_fraction = snakemake.params.cov_fraction
    db_dir = snakemake.params.db_dir
    output = Path(snakemake.params.vscores)
    all_hmm_results = Path(snakemake.params.all_hmm_results)
    filtered_hmm_results = Path(snakemake.params.filtered_hmm_results)
    num_threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    minscore = snakemake.params.min_bitscore
    min_bitscore_fraction = snakemake.params.min_bitscore_fraction_heuristic
    evalue = snakemake.params.max_evalue
    keep_full_results = snakemake.params.keep_full_hmm_results

    logger.info("Protein HMM searches starting...")

    set_memory_limit(mem_limit)
    logger.debug(f"Memory limit set to {mem_limit} GB")

    priority_order = ["KEGG", "FOAM", "PHROG", "VOG", "Pfam", "eggNOG", "dbCAN", "CAMPER", "METABOLIC", "user_custom"]
    hmm_paths = sorted(
        [Path(db_dir) / f for f in os.listdir(db_dir) if f.endswith((".H3M", ".h3m"))],
        key=lambda x: priority_order.index(assign_db(x)) if assign_db(x) in priority_order else float("inf"),
    )

    resumable = find_resumable_tmp_dir(wdir)
    tmp_dir = None
    manifest = None

    if resumable is not None:
        manifest = load_manifest(resumable)
        if manifest is None:
            logger.info(f"Found previous tmp dir (no manifest): {resumable}")
            tmp_dir = resumable
        else:
            same_inputs = (
                manifest.get("protein_dir") == str(protein_dir)
                and manifest.get("db_dir") == str(db_dir)
                and manifest.get("cov_fraction") == float(cov_fraction)
                and manifest.get("minscore") == float(minscore)
                and manifest.get("min_bitscore_fraction") == float(min_bitscore_fraction)
                and manifest.get("evalue") == float(evalue)
            )
            if same_inputs:
                logger.info(f"Resuming from tmp dir: {resumable}")
                tmp_dir = resumable
            else:
                logger.info(f"Found tmp dir but manifest does not match current inputs: {resumable}")

    if tmp_dir is None:
        run_id = f"{datetime.now().strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
        tmp_dir = wdir / f"hmmsearch_tmp_{run_id}"
        tmp_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Creating tmp dir: {tmp_dir}")

    aggregated = aggregate_sequences(protein_dir)
    seq_lengths = {rec.header.name: len(rec.seq) for rec in aggregated}
    N = len(aggregated)
    
    batch_files = []
    batch_keys = []

    manifest = load_manifest(tmp_dir)
    if manifest is not None and "batch_keys" in manifest and "batch_files" in manifest:
        batch_keys = list(manifest["batch_keys"])
        batch_files = list(manifest["batch_files"])
        missing_batches = [bf for bf in batch_files if not Path(bf).exists()]
        if missing_batches:
            logger.warning("Manifest found but some batch FASTAs are missing; regenerating batches.")
            batch_files = []
            batch_keys = []

    if not batch_keys:
        chunk_size = determine_chunk_size(N, mem_limit, est_bytes_per_seq=32768, max_chunk_fraction=0.8)
        for idx, batch in enumerate(split_aggregated_sequences(aggregated, chunk_size)):
            batch_key = f"seqbatch_{idx}"
            batch_fasta = tmp_dir / f"{batch_key}.faa"

            if not batch_fasta.exists() or batch_fasta.stat().st_size == 0:
                remove_incomplete_outputs([batch_fasta])
                with open(batch_fasta, "w", encoding="utf-8") as f:
                    for rec in batch:
                        write_fasta(rec, f)

            batch_files.append(str(batch_fasta))
            batch_keys.append(batch_key)

        logger.info(f"Splitting {N:,} sequences into {len(batch_keys):,} batches of {chunk_size:,}")

        manifest = {
            "created": datetime.now().isoformat(timespec="seconds"),
            "protein_dir": str(protein_dir),
            "db_dir": str(db_dir),
            "n_sequences": int(N),
            "cov_fraction": float(cov_fraction),
            "minscore": float(minscore),
            "min_bitscore_fraction": float(min_bitscore_fraction),
            "evalue": float(evalue),
            "batch_keys": batch_keys,
            "batch_files": batch_files,
            "hmm_paths": [str(p) for p in hmm_paths],
        }
        write_manifest(tmp_dir, manifest)
    else:
        logger.info(f"Using existing batches from manifest: {len(batch_keys):,} batches")

    result_paths = []
    for db_path in hmm_paths:
        db_name = assign_db(db_path)
        logger.info(f"Running HMMsearches against {db_name} profile HMMs...")

        if len(batch_keys) > 1:
            for batch_key, batch_fasta in tqdm(
                list(zip(batch_keys, batch_files)),
                total=len(batch_keys),
                desc=f"HMMsearches ({db_name})",
                unit="batch",
            ):
                outfile = tmp_dir / f"{db_path.stem}_{batch_key}_search.tsv"
                if is_complete(outfile):
                    logger.info(f"Skipping (complete): {outfile.name}")
                    result_paths.append((str(outfile), db_path))
                    continue

                remove_incomplete_outputs([outfile])
                out_path = hmmsearch_serial(
                    outfile=outfile,
                    batch_fasta=batch_fasta,
                    db_path=str(db_path),
                    seq_lengths=seq_lengths,
                    min_coverage=cov_fraction,
                    min_score=minscore,
                    min_bitscore_fraction=min_bitscore_fraction,
                    evalue=evalue,
                    cpus=num_threads,
                )
                result_paths.append((out_path, db_path))
        else:
            for batch_key, batch_fasta in zip(batch_keys, batch_files):
                outfile = tmp_dir / f"{db_path.stem}_{batch_key}_search.tsv"
                if is_complete(outfile):
                    logger.info(f"Skipping (complete): {outfile.name}")
                    result_paths.append((str(outfile), db_path))
                    continue

                remove_incomplete_outputs([outfile])
                out_path = hmmsearch_serial(
                    outfile=outfile,
                    batch_fasta=batch_fasta,
                    db_path=str(db_path),
                    seq_lengths=seq_lengths,
                    min_coverage=cov_fraction,
                    min_score=minscore,
                    min_bitscore_fraction=min_bitscore_fraction,
                    evalue=evalue,
                    cpus=num_threads,
                )
                result_paths.append((out_path, db_path))

    logger.info("Filtering HMMsearch results...")
    full_paths = []
    best_paths = []
    for result_path, db_path in result_paths:
        db = assign_db(db_path)
        logger.info(f"Filtering results from {db}...")
        result_path = Path(result_path)
        full_path = Path(str(result_path).replace("_search.tsv", "_full.tsv"))
        best_path = Path(str(result_path).replace("_search.tsv", "_best_kept.tsv"))

        full_done = is_complete(full_path)
        best_done = is_complete(best_path)

        if full_done and best_done:
            logger.info(f"Skipping standardize (complete): {full_path.name} and {best_path.name}")
        else:
            logger.debug(f"Standardizing {result_path} to {full_path} and {best_path}")
            standardize_and_filter_hmm_results(result_path, db_path, full_path, best_path)

        full_paths.append(str(full_path))
        best_paths.append(str(best_path))
    
    logger.info("Aggregating final HMMsearch results...")
    # Load + concat full and best
    schema_overrides_full = {
        "hmm_id": pl.Utf8,
        "db": pl.Utf8,
        "sequence": pl.Utf8,
        "evalue": pl.Float64,
        "score": pl.Float64,
        "alignment_type": pl.Utf8,
        "coverage_sequence": pl.Float64,
        "coverage_hmm": pl.Float64,
        "keep": pl.Boolean,
        "note": pl.Utf8,
    }

    schema_overrides_best = {
        "hmm_id": pl.Utf8,
        "db": pl.Utf8,
        "sequence": pl.Utf8,
        "evalue": pl.Float64,
        "score": pl.Float64,
        "alignment_type": pl.Utf8,
        "coverage_sequence": pl.Float64,
        "coverage_hmm": pl.Float64,
    }

    if keep_full_results:
        logger.debug("Loading and concatenating full HMMsearch results...")
        dfs_full = []
        for p in full_paths:
            logger.debug(f"Loading full HMMsearch results from {p}")
            try:
                dfs_full.append(pl.read_csv(p, separator="\t", schema_overrides=schema_overrides_full))
            except Exception as e:
                logger.warning(f"Failed to load {p}: {e}")
        logger.debug("Concatenating full HMMsearch results...")
        combined_full_df = pl.concat(dfs_full) if dfs_full else pl.DataFrame(schema=schema_overrides_full)
        combined_full_df.write_csv(all_hmm_results, separator="\t")

    logger.debug("Loading and concatenating best-kept HMMsearch results...")
    dfs_best = []
    for p in best_paths:
        logger.debug(f"Loading best-kept HMMsearch results from {p}")
        try:
            dfs_best.append(pl.read_csv(p, separator="\t", schema_overrides=schema_overrides_best))
        except Exception as e:
            logger.warning(f"Failed to load {p}: {e}")
    logger.debug("Concatenating best-kept HMMsearch results...")
    combined_best_df = pl.concat(dfs_best) if dfs_best else pl.DataFrame(schema=schema_overrides_best)
    combined_best_df.write_csv(filtered_hmm_results, separator="\t")

    # Use best-kept for vscore assignment downstream
    logger.debug("Loading HMM V-scores...")
    vscores_df = pl.read_csv(
        hmm_vscores,
        schema_overrides={"id": pl.Utf8, "V-score": pl.Float64, "VL-score": pl.Float64, "db": pl.Categorical, "name": pl.Utf8},
        separator="\t",
    ).with_columns(
        [
            pl.when(pl.col("db") == "Pfam")
            .then(pl.col("id").str.replace(r"\.\d+$", ""))
            .otherwise(pl.col("id"))
            .alias("id_norm")
        ]
    )

    logger.debug("Merging best HMM results with V-scores...")
    combined_best_df = combined_best_df.with_columns(
        [
            pl.when(pl.col("db") == "Pfam")
            .then(pl.col("hmm_id").str.replace(r"\.\d+$", ""))
            .otherwise(pl.col("hmm_id"))
            .alias("hmm_id_norm")
        ]
    )

    merged_df = combined_best_df.rename({"hmm_id": "id"}).join(
        vscores_df, left_on="hmm_id_norm", right_on="id_norm", how="left"
    ).with_columns(
        [pl.col("id").alias("hmm_id")]
    )

    merged_df = merged_df.filter(pl.col("V-score").is_not_null())

    cols_to_drop = ["name", "db_right", "id", "id_norm", "hmm_id_norm"]
    for col in cols_to_drop:
        if col in merged_df.columns:
            merged_df = merged_df.drop(col)

    merged_df = merged_df.sort(["sequence", "score", "V-score", "db"])
    merged_df.write_csv(output, separator="\t")

    logger.debug("Cleaning up temporary files and HMMsearch directory...")
    for f in tmp_dir.iterdir():
        if f.name == "manifest.json":
            continue
        try:
            if f.is_file():
                f.unlink()
        except Exception:
            pass
    try:
        (tmp_dir / "manifest.json").unlink()
    except Exception:
        pass
    try:
        tmp_dir.rmdir()
    except Exception:
        pass

    logger.info("Protein HMM searches completed.")
    if keep_full_results:
        logger.info(f"Full results: {all_hmm_results}")
    logger.info(f"Best hits per protein: {filtered_hmm_results}")

if __name__ == "__main__":
    main()
