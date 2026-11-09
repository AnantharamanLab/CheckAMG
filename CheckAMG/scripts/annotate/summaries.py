import os
import textwrap
from datetime import datetime, timezone
from importlib.metadata import version, PackageNotFoundError

import polars as pl

KEYS = ["Protein", "Contig", "Genome"]
AVG_TYPE = {"metabolic": "AMG", "physiological": "APG", "regulatory": "AReG"}
CONFIDENCES = ["high", "medium", "low"]

EXC_LABELS = {
    "essential": "an essential or core (non-auxiliary) function",
    "glucan": "a glycoside hydrolase or glycosyltransferase function",
    "nucleotide": "a nucleotide metabolism function",
    "methyl": "a methyltransferase function",
    "lipid": "a lipid metabolism function",
}


def checkamg_version() -> str:
    try:
        return version("checkamg")
    except PackageNotFoundError:
        return "unknown"


def read_table(path) -> pl.DataFrame:
    path = str(path)
    if path.endswith(".parquet"):
        return pl.read_parquet(path)
    return pl.read_csv(path, separator="\t", infer_schema_length=100000)


def write_table(df: pl.DataFrame, path) -> None:
    path = str(path)
    if path.endswith(".parquet"):
        df.write_parquet(path)
    else:
        df.write_csv(path, separator="\t")


def read_counts(path: str) -> dict:
    if not path or not os.path.exists(path):
        return {}
    df = pl.read_csv(path, separator="\t", schema={"metric": pl.Utf8, "value": pl.Utf8})
    return dict(zip(df["metric"].to_list(), df["value"].to_list()))


def _normalize_ref_id(id_expr: pl.Expr) -> pl.Expr:
    id_expr = id_expr.cast(pl.Utf8, strict=False).str.strip_chars()
    return pl.when(id_expr.str.starts_with("PF")).then(id_expr.str.replace(r"\.\d+$", "", literal=False)).otherwise(id_expr)


def max_amg_weight_per_protein(df: pl.DataFrame, amg_ref: pl.DataFrame) -> pl.DataFrame:
    weights = (
        amg_ref
        .with_columns(_normalize_ref_id(pl.col("id")).alias("id"), pl.col("AMG_weight").cast(pl.Float64, strict=False))
        .sort("AMG_weight", descending=True, nulls_last=True)
        .group_by("id", maintain_order=True).agg(pl.col("AMG_weight").first(), pl.col("name").first().alias("_ref_name"))
    )
    id_cols = {
        "KEGG_hmm_id": None,
        "FOAM_hmm_id": None,
        "Pfam_hmm_id": r"\.\d+$",
        "dbCAN_hmm_id": r"_(.*)",
        "METABOLIC_hmm_id": None,
        "CAMPER_hmm_id": None,
        "PHROG_hmm_id": None,
    }
    parts = []
    for col, pattern in id_cols.items():
        if col not in df.columns:
            continue
        expr = pl.col(col).cast(pl.Utf8, strict=False)
        if pattern:
            expr = expr.str.replace(pattern, "", literal=False)
        parts.append(df.select(KEYS + [expr.alias("id")]).drop_nulls("id"))
    if not parts:
        return df.select(KEYS).with_columns(pl.lit(None, dtype=pl.Float64).alias("AMG Weight"), pl.lit(None, dtype=pl.Utf8).alias("_ref_name")).unique()
    return (
        pl.concat(parts)
        .join(weights, on="id", how="inner")
        .sort("AMG_weight", descending=True, nulls_last=True)
        .group_by(KEYS, maintain_order=True).agg(pl.col("AMG_weight").first().alias("AMG Weight"), pl.col("_ref_name").first())
    )


def readable_reason(token: str, min_amg_weight: float, amg_weight) -> str:
    tok = token.strip()
    if tok == "below_min_amg_weight":
        if amg_weight is not None:
            return f"AMG weight {amg_weight:.3f} is below the minimum of {min_amg_weight}"
        return f"AMG weight is below the minimum of {min_amg_weight}"
    if tok.startswith("hard|"):
        parts = tok.split("|")
        label = EXC_LABELS.get(parts[1] if len(parts) > 1 else "", "a filtered function")
        hit = " ".join(parts[2:4]) if len(parts) >= 4 else ""
        return f"Annotation matches {label}" + (f" ({hit})" if hit else "")
    if tok == "context|nonviral_region":
        return "Outside a strict viral region (--filter-ambig-regions)"
    if tok.startswith("context|avg_array"):
        n = tok.rsplit("_", 1)[-1]
        if n.isdigit():
            return f"In a run of {n} or more adjacent AVG candidates (--filter-avg-arrays)"
        return "In a run of adjacent AVG candidates (--filter-avg-arrays)"
    return tok


def build_filtered_candidates(audits: dict, final_df: pl.DataFrame, amg_ref: pl.DataFrame, min_amg_weight: float) -> pl.DataFrame:
    avg_keys = final_df.filter(pl.col("Protein Classification") != "unclassified").select(KEYS)
    conf = final_df.select(KEYS + ["Protein Viral Origin Confidence"])
    frames = []
    for avg_type, audit in audits.items():
        if audit is None or audit.height == 0:
            continue
        removed = audit.filter(pl.col("removed").cast(pl.Boolean, strict=False).fill_null(False))
        if removed.height == 0:
            continue
        removed = removed.unique(subset=KEYS, keep="first")
        if avg_type == "AMG":
            removed = (
                removed.join(max_amg_weight_per_protein(removed, amg_ref), on=KEYS, how="left")
                .with_columns(pl.coalesce([pl.col("Function"), pl.col("_ref_name")]).alias("Function"))
                .drop("_ref_name")
            )
        else:
            removed = removed.with_columns(pl.lit(None, dtype=pl.Float64).alias("AMG Weight"))
        reasons = [
            "; ".join(
                readable_reason(t, min_amg_weight, w)
                for t in str(r or "").split(",") if t.strip()
            ) or "Unknown"
            for r, w in zip(removed["remove_reason"].to_list(), removed["AMG Weight"].to_list())
        ]
        frames.append(
            removed.select(KEYS + ["Function", "AMG Weight"])
            .with_columns(
                pl.lit(avg_type).alias("Candidate AVG Type"),
                pl.Series("Reason Filtered", reasons, dtype=pl.Utf8),
            )
        )
    cols = KEYS + ["Candidate AVG Type", "Function", "Protein Viral Origin Confidence", "AMG Weight", "Reason Filtered"]
    if not frames:
        return pl.DataFrame(schema={c: (pl.Float64 if c == "AMG Weight" else pl.Utf8) for c in cols})
    out = (
        pl.concat(frames, how="diagonal_relaxed")
        .join(avg_keys, on=KEYS, how="anti")
        .join(conf, on=KEYS, how="left")
        .select(cols)
        .with_columns(pl.col("AMG Weight").round(3))
        .sort(["Candidate AVG Type", "Genome", "Contig", "Protein"])
    )
    return out


CAT_TOP = "Functional Categories (Level 1)"
CAT_PATHS = "Functional Category Paths (L1 > L2 > L3)"


def collapse_categories(cat_df: pl.DataFrame) -> pl.DataFrame:
    levels = ["category_L1", "category_L2", "category_L3"]
    cat_df = cat_df.with_columns(
        [pl.when(pl.col(c).str.strip_chars() == "").then(None).otherwise(pl.col(c)).alias(c) for c in levels]
    ).with_columns(
        pl.when(pl.col("category_L1").is_not_null())
        .then(pl.concat_str(levels, separator=" > ", ignore_nulls=True))
        .alias("_path")
    )
    collapsed = cat_df.group_by(KEYS).agg(
        pl.col("category_L1").drop_nulls().unique().sort().str.join("; ").alias(CAT_TOP),
        pl.col("_path").drop_nulls().unique().sort().str.join("; ").alias(CAT_PATHS),
    )
    return collapsed.with_columns(
        [pl.when(pl.col(c) == "").then(pl.lit(None, dtype=pl.Utf8)).otherwise(pl.col(c)).alias(c) for c in [CAT_TOP, CAT_PATHS]]
    )


def build_avg_predictions(final_df: pl.DataFrame, curated: dict, categories: dict) -> pl.DataFrame:
    frames = []
    for classification, avg_type in AVG_TYPE.items():
        base = final_df.filter(pl.col("Protein Classification") == classification)
        cur = curated[classification]
        support_cols = [c for c in ["Supporting_HMM_db", "Supporting_HMM_id", "Supporting_HMM_bitscore", "AMG_weight"] if c in cur.columns]
        support = cur.select(KEYS + support_cols).unique(subset=KEYS, keep="first")
        cats = collapse_categories(categories[classification])
        frames.append(
            base.join(support, on=KEYS, how="left")
            .join(cats, on=KEYS, how="left")
            .with_columns(pl.lit(avg_type).alias("AVG Type"))
        )
    df = pl.concat(frames, how="diagonal_relaxed")
    if "AMG_weight" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("AMG_weight"))
    for c in ["Supporting_HMM_db", "Supporting_HMM_id"]:
        if c not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8).alias(c))
    if "Supporting_HMM_bitscore" not in df.columns:
        df = df.with_columns(pl.lit(None, dtype=pl.Float64).alias("Supporting_HMM_bitscore"))
    df = df.rename({
        "Protein Viral Origin Confidence": "Viral Origin Confidence",
        "Protein Viral Origin Probability": "Viral Origin Probability",
        "AMG_weight": "AMG Weight",
        "Supporting_HMM_db": "Supporting HMM Database",
        "Supporting_HMM_id": "Supporting HMM ID",
        "Supporting_HMM_bitscore": "Supporting HMM Bit Score",
    })
    conf_order = {"high": 0, "medium": 1, "low": 2}
    type_order = {"AMG": 0, "APG": 1, "AReG": 2}
    return (
        df.select([
            "Protein", "Contig", "Genome",
            "AVG Type", "Function",
            "Viral Origin Confidence", "Viral Origin Probability",
            "AMG Weight",
            CAT_TOP, CAT_PATHS,
            "Supporting HMM Database", "Supporting HMM ID", "Supporting HMM Bit Score",
            "Protein in Strict Viral Region",
        ])
        .with_columns(
            pl.col("AVG Type").replace_strict(type_order, return_dtype=pl.Int32).alias("_t"),
            pl.col("Viral Origin Confidence").replace_strict(conf_order, default=3, return_dtype=pl.Int32).alias("_c"),
        )
        .sort(["_t", "_c", "Genome", "Contig", "Protein"])
        .drop(["_t", "_c"])
        .with_columns(
            pl.col("Viral Origin Probability").cast(pl.Float64).round(3),
            pl.col("AMG Weight").cast(pl.Float64).round(3),
            pl.col("Supporting HMM Bit Score").cast(pl.Float64).round(1),
        )
    )


def build_genome_summary(final_df: pl.DataFrame, avg_pred: pl.DataFrame) -> pl.DataFrame:
    base = final_df.group_by("Genome").agg(
        pl.col("Contig").n_unique().alias("Contigs"),
        pl.col("Protein").count().alias("Proteins"),
    )
    count_exprs = [
        ((pl.col("AVG Type") == t) & (pl.col("Viral Origin Confidence") == c)).sum().cast(pl.Int64).alias(f"{t}s ({c})")
        for t in ["AMG", "APG", "AReG"] for c in CONFIDENCES
    ]
    counts = avg_pred.group_by("Genome").agg(count_exprs) if avg_pred.height else None
    high_fns = (
        avg_pred.filter(pl.col("Viral Origin Confidence") == "high")
        .group_by("Genome")
        .agg(
            (pl.col("AVG Type") + pl.lit(": ") + pl.col("Function").fill_null("unknown function"))
            .unique().sort().str.join("; ").alias("High-Confidence AVG Functions")
        )
    )
    out = base
    if counts is not None:
        out = out.join(counts, on="Genome", how="left")
    else:
        out = out.with_columns([pl.lit(0, dtype=pl.Int64).alias(f"{t}s ({c})") for t in ["AMG", "APG", "AReG"] for c in CONFIDENCES])
    count_cols = [f"{t}s ({c})" for t in ["AMG", "APG", "AReG"] for c in CONFIDENCES]
    out = out.with_columns([pl.col(c).fill_null(0) for c in count_cols])
    out = out.join(high_fns, on="Genome", how="left")
    return out.select(["Genome", "Contigs", "Proteins"] + count_cols + ["High-Confidence AVG Functions"]).sort("Genome")


RUN_PARAMETER_LABELS = [
    ("min_len", "Minimum contig length (bp)"),
    ("min_cds", "Minimum ORFs per contig"),
    ("annotation_percent_threshold", "Minimum annotated fraction for context analysis"),
    ("cov_fraction", "Minimum HMM and query coverage"),
    ("max_evalue", "Fallback maximum E-value"),
    ("min_bitscore", "Fallback minimum bit score"),
    ("min_bitscore_fraction_heuristic", "Bit score fraction heuristic"),
    ("window_size", "VL-score window size (bp)"),
    ("minimum_flank_vscore", "Minimum flank V-score"),
    ("minimum_window_avg_vlscore", "Minimum window average VL-score"),
    ("use_hallmark", "Use hallmark genes for flanks"),
    ("min_amg_weight", "Minimum AMG weight"),
    ("filter_presets", "Filter presets"),
    ("filter_nonviral_regions", "Filter predictions outside strict viral regions"),
    ("filter_avg_arrays", "Filter AVG arrays"),
    ("avg_array_len_limit", "AVG array length limit"),
]


def build_run_summary(run_params: dict, counts_length: dict, counts_cds: dict, final_df: pl.DataFrame, avg_pred: pl.DataFrame, n_filtered: int) -> pl.DataFrame:
    rows = []

    def add(section, metric, value):
        rows.append((section, metric, "" if value is None else str(value)))

    add("run", "CheckAMG version", checkamg_version())
    add("run", "Module", "annotate")
    add("run", "Date (UTC)", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M"))
    add("run", "Annotate database", os.path.basename(str(run_params["db_dir"]).rstrip("/")))
    add("run", "Input type", run_params["input_type"])

    for key, label in RUN_PARAMETER_LABELS:
        add("parameters", label, run_params[key])

    if counts_length:
        add("inputs", "Input genomes", counts_length.get("input_genomes"))
        add("inputs", "Input contigs", counts_length.get("input_contigs"))
        add("inputs", "Contigs passing minimum length", counts_length.get("contigs_passing_min_length"))
        add("inputs", "Genomes passing minimum length", counts_length.get("genomes_passing_min_length"))
    if counts_cds:
        if not counts_length:
            add("inputs", "Input genomes", counts_cds.get("genomes_before_min_orf"))
            add("inputs", "Input contigs", counts_cds.get("contigs_before_min_orf"))
        add("inputs", "Contigs passing minimum ORFs", counts_cds.get("contigs_passing_min_orf"))
        add("inputs", "Genomes passing minimum ORFs", counts_cds.get("genomes_passing_min_orf"))
    add("inputs", "Genomes analyzed", final_df["Genome"].n_unique())
    add("inputs", "Contigs analyzed", final_df["Contig"].n_unique())
    add("inputs", "Proteins analyzed", final_df.height)

    for t in ["AMG", "APG", "AReG"]:
        sub = avg_pred.filter(pl.col("AVG Type") == t)
        add("predictions", f"{t}s (all confidence levels)", sub.height)
        for c in CONFIDENCES:
            add("predictions", f"{t}s ({c} confidence)", sub.filter(pl.col("Viral Origin Confidence") == c).height)
    add("predictions", "AVG candidates removed by curation", n_filtered)

    return pl.DataFrame(rows, schema=["Section", "Metric", "Value"], orient="row")


def _fmt(n) -> str:
    try:
        return f"{int(n):,}"
    except (TypeError, ValueError):
        return "NA"


def _annotate_run_block(final_df: pl.DataFrame, avg_pred: pl.DataFrame, n_filtered: int, counts_length: dict, counts_cds: dict) -> str:
    n_input = counts_length.get("input_genomes") or counts_cds.get("genomes_before_min_orf")
    by_type = ", ".join(f"{_fmt(avg_pred.filter(pl.col('AVG Type') == t).height)} {t}s" for t in ["AMG", "APG", "AReG"])
    by_conf = ", ".join(f"{c} {_fmt(avg_pred.filter(pl.col('Viral Origin Confidence') == c).height)}" for c in CONFIDENCES)
    lines = [
        f"Genomes analyzed:   {_fmt(final_df['Genome'].n_unique())}" + (f" (of {_fmt(n_input)} input genomes)" if n_input else ""),
        f"Proteins analyzed:  {_fmt(final_df.height)}",
        f"AVGs predicted:     {_fmt(avg_pred.height)} ({by_type})",
        f"  by viral-origin confidence: {by_conf}",
        f"AVG candidates removed by curation: {_fmt(n_filtered)} (see detailed/filtered_candidates.tsv)",
    ]
    return "\n".join(lines)


README_TEXT = """CheckAMG annotate results
=========================

CheckAMG version: {version}

This run
--------
{run_block}

Start here
----------
1. avg_predictions.tsv    One row per predicted auxiliary viral gene (AVG). This is the main result.
2. genome_summary.tsv     One row per genome, with AVG counts by type and viral-origin confidence.
3. faa/                   Amino-acid sequences of predicted AVGs: AMGs_*, APGs_*, AReGs_*, split by
                          viral-origin confidence (high, medium, low) plus *_all.faa.
4. detailed/              Supporting tables:
                            all_proteins.tsv         One row per analyzed protein, including proteins
                                                     that are not AVGs.
                            filtered_candidates.tsv  Proteins that matched an AVG reference but were
                                                     removed during curation, with the reason (one row
                                                     per protein and candidate AVG type).
                            run_summary.tsv          What was run (version, database, parameters), how
                                                     many inputs were analyzed or removed by the
                                                     length/ORF filters, and how many AVGs were found.
                          plus the full tables used to build the files above (annotations, genome
                          context features, curated tables, long-format categories, filter audits).

If this is part of a CheckAMG end-to-end run, the combined results of annotate and de-novo are in
03_aggregate/ (see 03_aggregate/README.txt).

AVG types
---------
AMG   auxiliary metabolic gene
APG   auxiliary physiology gene
AReG  auxiliary regulatory gene

Choosing a viral-origin confidence level
----------------------------------------
Every AVG is reported with one of three levels. Choose which to keep based on what your input is:
  high     Recommended for mixed-community metagenomes, or when input composition is unknown.
  medium   Can be added to high for viromes, confident predictions from virus-prediction tools,
           and isolate or complete viral genomes.
  low      Not recommended for any analysis. These proteins lack sufficient support for viral
           origin and may be encoded by host or other non-viral sequence.
Counting every row of avg_predictions.tsv without filtering by confidence overstates the number of
virus-encoded AVGs.

Columns in avg_predictions.tsv
------------------------------
Protein, Contig, Genome     Identifiers. For multi-contig bins/vMAGs, Genome is the bin name.
AVG Type                    AMG, APG, or AReG.
Function                    Description of the reference family that supports the AVG call. Among
                            the protein's HMM hits to reference families of its AVG type (for AMGs,
                            only families at or above the minimum AMG weight), the hit with the
                            highest bit score is used. A protein that qualifies as more than one
                            AVG type is reported once in this table and in detailed/all_proteins.tsv,
                            with priority AMG, then APG, then AReG (it still appears in each
                            matching FASTA file in faa/).
Viral Origin Confidence     high, medium, or low likelihood that the protein is virus-encoded,
                            from the genome-context LightGBM classifier. Refers to viral origin,
                            not to the functional annotation.
Viral Origin Probability    LightGBM probability (0-1) underlying the confidence level.
AMG Weight                  AMGs only. AMG weight (0-1) of the reference family that supplied
                            Function. It combines how specifically the family is associated with
                            metabolism (its metabolic ratio, from KEGG BRITE, FOAM, GO, and EC
                            classifications) with how rare it is in viral genomes (its VL-score).
                            Higher values mean more specifically metabolic and less common in
                            viruses. AMGs require at least the minimum AMG weight
                            ({min_amg_weight} in this run, default 0.6).
Functional Categories (Level 1)
                            Top-level functional categories of the reference family that
                            supplied Function, separated by "; ". Useful for filtering.
Functional Category Paths (L1 > L2 > L3)
                            Every category assignment of that family as a full path, separated
                            by "; ". A family can belong to several pathways. This lists
                            potential pathways and does not mean the gene acts in all (or any)
                            of them. For counts per category, use the long-format
                            detailed/*_gene_categories.tsv files (one row per assignment).
Supporting HMM Database     Database of the HMM hit that supports the AVG call.
Supporting HMM ID           ID of that HMM.
Supporting HMM Bit Score    Bit score of that hit.
Protein in Strict Viral Region
                            True if the protein lies in a strict viral region, as identified by
                            CheckAMG annotate's V- and VL-score-based viral region algorithm (see
                            "What are strict viral regions?" in the wiki FAQ). This is an optional
                            stringency flag. Viral Origin Confidence is the primary viral-origin
                            call. Use this flag only for extra-stringent filtering (the
                            --filter-ambig-regions option applies it automatically).

Columns in detailed/all_proteins.tsv
------------------------------------
Protein Classification      metabolic (AMG), physiological (APG), regulatory (AReG), or
                            unclassified.
Annotation Status           AVG; AVG candidate, filtered (see detailed/filtered_candidates.tsv);
                            Annotated, not an AVG; or Unannotated.
Protein Viral Origin Confidence / Probability
                            As above, for every protein.
Protein in Strict Viral Region
                            As above.
Function                    For AVGs, chosen as described for avg_predictions.tsv. For all other
                            proteins, the best-scoring HMM annotation across all databases.
<database> columns          Best hit per database (KEGG, FOAM, Pfam, dbCAN/CAZy, METABOLIC,
                            CAMPER, PHROG).
Best Scoring HMM, Best Scoring HMM Annotation, Best Scoring HMM Origin
                            Highest-scoring hit across all databases. This can differ from the
                            hit supporting an AVG call.

Notes on other detailed/ tables
-------------------------------
*_score columns are HMM bit scores. *_coverage columns are the fraction of the HMM profile
covered by the alignment. *_V-score and *_VL-score columns are viral-likeness scores of the
annotation (higher = more virus-like). Distances are in base pairs.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
FAQ: https://github.com/AnantharamanLab/CheckAMG/wiki/FAQ
"""


def write_readme(path: str, min_amg_weight: float, final_df: pl.DataFrame, avg_pred: pl.DataFrame, n_filtered: int, counts_length: dict, counts_cds: dict) -> None:
    block = _annotate_run_block(final_df, avg_pred, n_filtered, counts_length, counts_cds)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(README_TEXT.format(version=checkamg_version(), run_block=block, min_amg_weight=min_amg_weight))


FINAL_CONFIDENCES = ["very high", "high", "medium", "low"]
EVIDENCE_CLASSES = ["annotate and de-novo", "annotate only", "de-novo only"]
SUPPORT_COLS = ["AMG Weight", "Supporting HMM Database", "Supporting HMM ID", "Supporting HMM Bit Score"]


def build_aggregate_avg_predictions(all_df: pl.DataFrame, annotate_avg: pl.DataFrame = None) -> pl.DataFrame:
    df = all_df.filter(pl.col("Evidence Class") != "not supported")
    if annotate_avg is not None:
        support = annotate_avg.select(KEYS + [c for c in SUPPORT_COLS if c in annotate_avg.columns]).unique(subset=KEYS, keep="first")
        df = df.join(support, on=KEYS, how="left")
    for c in SUPPORT_COLS:
        if c not in df.columns:
            df = df.with_columns(pl.lit(None, dtype=pl.Utf8 if c in ["Supporting HMM Database", "Supporting HMM ID"] else pl.Float64).alias(c))
    type_order = {"AMG": 0, "APG": 1, "AReG": 2}
    conf_order = {c: i for i, c in enumerate(FINAL_CONFIDENCES)}
    return (
        df.with_columns(
            pl.col("AVG Type").replace_strict(type_order, default=3, return_dtype=pl.Int32).alias("_t"),
            pl.col("Final AVG Confidence Level").replace_strict(conf_order, default=4, return_dtype=pl.Int32).alias("_c"),
            pl.col("Final AVG Probability").cast(pl.Float64).round(3),
            pl.col("AMG Weight").cast(pl.Float64).round(3),
            pl.col("Supporting HMM Bit Score").cast(pl.Float64).round(1),
        )
        .sort(["_t", "_c", "Genome", "Contig", "Protein"])
        .select([
            "Protein", "Contig", "Genome",
            "Evidence Class", "AVG Type", "Function",
            "Final AVG Confidence Level", "Final AVG Probability",
            "Viral Origin Confidence (annotate)",
            "AMG Weight",
            CAT_TOP, CAT_PATHS,
            "Supporting HMM Database", "Supporting HMM ID", "Supporting HMM Bit Score",
            "Protein in Strict Viral Region",
        ])
    )


def build_aggregate_genome_summary(all_df: pl.DataFrame, avg_pred: pl.DataFrame) -> pl.DataFrame:
    base = all_df.group_by("Genome").agg(
        pl.col("Contig").n_unique().alias("Contigs"),
        pl.col("Protein").count().alias("Proteins"),
    )
    conf_cols = [f"AVGs ({c} final confidence)" for c in FINAL_CONFIDENCES]
    type_cols = ["AMGs", "APGs", "AReGs", "AVGs (no annotate type)"]
    ev_cols = [f"AVGs ({e})" for e in EVIDENCE_CLASSES]
    exprs = (
        [pl.len().cast(pl.Int64).alias("AVGs")]
        + [(pl.col("Final AVG Confidence Level") == c).sum().cast(pl.Int64).alias(n) for c, n in zip(FINAL_CONFIDENCES, conf_cols)]
        + [(pl.col("AVG Type") == t).sum().cast(pl.Int64).alias(f"{t}s") for t in ["AMG", "APG", "AReG"]]
        + [pl.col("AVG Type").is_null().sum().cast(pl.Int64).alias("AVGs (no annotate type)")]
        + [(pl.col("Evidence Class") == e).sum().cast(pl.Int64).alias(n) for e, n in zip(EVIDENCE_CLASSES, ev_cols)]
    )
    count_cols = ["AVGs"] + conf_cols + type_cols + ev_cols
    high_fns = (
        avg_pred.filter(pl.col("Final AVG Confidence Level").is_in(["very high", "high"]) & pl.col("Function").is_not_null())
        .group_by("Genome")
        .agg(
            (pl.col("AVG Type").fill_null("AVG") + pl.lit(": ") + pl.col("Function"))
            .unique().sort().str.join("; ").alias("High-Confidence AVG Functions")
        )
    )
    out = base
    if avg_pred.height:
        out = out.join(avg_pred.group_by("Genome").agg(exprs), on="Genome", how="left")
    else:
        out = out.with_columns([pl.lit(0, dtype=pl.Int64).alias(c) for c in count_cols])
    out = out.with_columns([pl.col(c).fill_null(0) for c in count_cols]).join(high_fns, on="Genome", how="left")
    return out.select(["Genome", "Contigs", "Proteins"] + count_cols + ["High-Confidence AVG Functions"]).sort("Genome")


def _aggregate_run_block(all_df: pl.DataFrame, avg_pred: pl.DataFrame) -> str:
    lines = [
        f"Proteins:  {_fmt(all_df.height)}",
        f"AVGs supported by annotate, de-novo, or both: {_fmt(avg_pred.height)}",
    ]
    for e in EVIDENCE_CLASSES:
        lines.append(f"  {e}: {_fmt(avg_pred.filter(pl.col('Evidence Class') == e).height)}")
    by_conf = ", ".join(f"{c} {_fmt(avg_pred.filter(pl.col('Final AVG Confidence Level') == c).height)}" for c in FINAL_CONFIDENCES)
    lines.append(f"  by final AVG confidence: {by_conf}")
    return "\n".join(lines)


AGGREGATE_README_TEXT = """CheckAMG aggregate results
==========================

CheckAMG version: {version}

This run
--------
{run_block}

Start here
----------
1. avg_predictions.tsv    One row per auxiliary viral gene (AVG) supported by annotate, de-novo, or
                          both. This is the main result.
2. genome_summary.tsv     One row per genome, with AVG counts by final AVG confidence, AVG type,
                          and evidence class.
3. detailed/              Supporting tables:
                            all_proteins.tsv                 One row per protein, including proteins
                                                             not supported as AVGs.
                            aggregated_results_detailed.tsv  Every probability and confidence level
                                                             from both modules, for every protein.
                            run_summary.tsv                  Inputs, support definitions, and counts.

AVG types
---------
AMG   auxiliary metabolic gene
APG   auxiliary physiology gene
AReG  auxiliary regulatory gene

Which confidence to use
-----------------------
Final AVG Confidence Level  From de-novo. Combines the probability that the protein is auxiliary
                            (AVG-like) with the probability that it is viral: very high, high,
                            medium, or low. When it is available, use it to judge AVG predictions.
Viral Origin Confidence (annotate)
                            From annotate's genome-context classifier: high, medium, or low. This
                            is supporting evidence for viral origin.

Evidence Class
--------------
{evidence_block}

Columns in avg_predictions.tsv
------------------------------
Protein, Contig, Genome     Identifiers. For multi-contig bins/vMAGs, Genome is the bin name.
Evidence Class              See above.
AVG Type                    AMG, APG, or AReG, from annotate. Empty if annotate did not classify the
                            protein as an AVG.
Function                    From annotate. For AVGs, the description of the reference family that
                            supports the AVG call. Otherwise, the best-scoring HMM annotation, or empty
                            if the protein has no annotation.
Final AVG Confidence Level  See above.
Final AVG Probability       de-novo probability (0-1) underlying the final AVG confidence level.
Viral Origin Confidence (annotate)
                            See above.
AMG Weight, Functional Categories (Level 1), Functional Category Paths (L1 > L2 > L3),
Supporting HMM Database, Supporting HMM ID, Supporting HMM Bit Score, Protein in Strict Viral Region
                            From annotate, as described in the annotate README.txt. Empty for proteins
                            without an annotate AVG type. Functional categories are "Unknown" for
                            proteins with no annotation and "Other" for annotated proteins with no
                            assigned category.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
FAQ: https://github.com/AnantharamanLab/CheckAMG/wiki/FAQ
"""


def _join_levels(levels) -> str:
    levels = list(levels)
    return levels[0] if len(levels) == 1 else ", ".join(levels[:-1]) + (", or " if len(levels) > 2 else " or ") + levels[-1]


def _evidence_block(annotate_supported_conf, denovo_supported_conf) -> str:
    entries = [
        ("annotate and de-novo", f"Annotate classifies the protein as an AVG with {_join_levels(annotate_supported_conf)} viral-origin confidence, and de-novo assigns a final AVG confidence of {_join_levels(denovo_supported_conf)}. These have the most support."),
        ("annotate only", "Supported by annotate but not by de-novo (Final AVG Confidence Level is below those levels, or empty if de-novo made no prediction for the protein)."),
        ("de-novo only", "Supported by de-novo but not by annotate. This includes proteins that annotate classified as AVGs at a lower viral-origin confidence, and proteins with no annotate AVG type (AVG Type is empty), which annotation-based prediction alone would miss."),
    ]
    lines = ["Every row of avg_predictions.tsv is supported by at least one module:"]
    for label, text in entries:
        lines.append(textwrap.fill(text, width=100, initial_indent=f"  {label:<22}", subsequent_indent=" " * 24))
    return "\n".join(lines)


def write_aggregate_readme(path: str, all_df: pl.DataFrame, avg_pred: pl.DataFrame, annotate_supported_conf, denovo_supported_conf) -> None:
    block = _evidence_block(annotate_supported_conf, denovo_supported_conf)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(AGGREGATE_README_TEXT.format(version=checkamg_version(), run_block=_aggregate_run_block(all_df, avg_pred), evidence_block=block))


END_TO_END_README_TEXT = """CheckAMG end-to-end results
===========================

CheckAMG version: {version}

Start with 03_aggregate/avg_predictions.tsv. It lists every auxiliary viral gene (AVG) supported by
annotate, de-novo, or both, and 03_aggregate/README.txt explains its columns.

01_annotate/    Output of CheckAMG annotate: annotation-based AVG predictions at every viral-origin
                confidence level (see 01_annotate/results/README.txt).
02_denovo/      Output of CheckAMG de-novo: protein language model predictions for every protein
                (predictions.tsv).
03_aggregate/   Combined results of annotate and de-novo. Start here.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
"""


def write_end_to_end_readme(path: str) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(END_TO_END_README_TEXT.format(version=checkamg_version()))
