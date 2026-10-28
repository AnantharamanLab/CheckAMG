import os
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

MIN_AMG_WEIGHT = 0.6
ALLOWED_FILTERS = {"glucan", "nucleotide", "methyl", "lipid"}
HMM_SOURCES = ["KEGG", "FOAM", "Pfam", "dbCAN", "METABOLIC", "CAMPER", "PHROG"]
LOG_REASONS = {"3+ AMGs in a row": "array_3plus", "AMG weight below minimum": "below_min_amg_weight",
               "other non-auxiliary function": "keyword_flag_filter"}


def _norm(expr):
    import polars as pl

    return pl.when(expr.str.starts_with("PF")).then(expr.str.replace(r"\.\d+$", "")).otherwise(expr)


def checkamg_reference(files_dir):
    """Id sets CheckAMG annotate filters on, read from its shipped reference tables."""
    import polars as pl

    files_dir = Path(files_dir)
    amgs = pl.read_csv(files_dir / "AMGs.tsv", separator="\t", infer_schema_length=0).with_columns(
        _norm(pl.col("id")).alias("id"), pl.col("AMG_weight").cast(pl.Float64, strict=False).fill_null(0.0))
    read_ids = lambda name: set(pl.read_csv(files_dir / name, separator="\t", infer_schema_length=0).with_columns(_norm(pl.col("id")).alias("id"))["id"].to_list())
    filters = pl.read_csv(files_dir / "AMG_filters.tsv", separator="\t", infer_schema_length=0).with_columns(_norm(pl.col("id")).alias("id"))
    active_cols = [c for c in filters.columns if c.startswith("filter_") and c[7:] not in ALLOWED_FILTERS]
    return {"all_amg": set(amgs["id"].to_list()), "pass_amg": set(amgs.filter(pl.col("AMG_weight") >= MIN_AMG_WEIGHT)["id"].to_list()),
            "phys": read_ids("APGs.tsv"), "reg": read_ids("AReGs.tsv"),
            "active": set(filters.filter(pl.any_horizontal([pl.col(c) == "True" for c in active_cols]))["id"].to_list())}


def _reconstruct(job):
    sample_dir, sample, wanted, reference, threads = job
    os.environ["POLARS_MAX_THREADS"] = str(threads)
    import polars as pl

    results = Path(sample_dir, "results")
    ga = pl.read_parquet(results / "gene_annotations.parquet", columns=["Protein", "Contig", "Genome"] + [f"{s}_hmm_id" for s in HMM_SOURCES])
    gc = pl.read_parquet(results / "genes_genomic_context.parquet", columns=["protein", "gene_number"]).rename({"protein": "Protein"}).unique("Protein")
    fr = pl.read_parquet(results / "final_results.parquet", columns=["Protein", "Protein Classification"]).unique("Protein")
    curated = set(pl.read_parquet(results / "metabolic_genes_curated.parquet", columns=["Protein"])["Protein"].to_list())
    ids = [_norm(pl.col("Pfam_hmm_id")) if s == "Pfam" else (pl.col("dbCAN_hmm_id").str.replace(r"_(.*)", "") if s == "dbCAN" else pl.col(f"{s}_hmm_id"))
           for s in HMM_SOURCES]
    ga = ga.join(gc, on="Protein", how="left").with_columns(pl.concat_list(ids).list.drop_nulls().alias("ids"))
    isin = lambda key: pl.col("ids").list.eval(pl.element().is_in(list(reference[key]))).list.any()
    ga = ga.with_columns(isin("all_amg").alias("cand"), isin("pass_amg").alias("wpass"), isin("phys").alias("is_phys"), isin("reg").alias("is_reg"),
                         isin("active").alias("flag_active"), (pl.col("ids").list.len() == 0).alias("no_hit"))
    ga = ga.with_columns((pl.col("wpass") | pl.col("is_phys") | pl.col("is_reg")).alias("member")).sort(["Genome", "Contig", "gene_number"])
    ga = ga.with_columns(pl.col("member").rle_id().over(["Genome", "Contig"]).alias("_rid"))
    ga = ga.with_columns((pl.col("member") & (pl.len().over(["Genome", "Contig", "_rid"]) >= 3)).alias("in_array_ge3"))
    ga = ga.with_columns(
        pl.when(~pl.col("cand")).then(pl.lit(None, dtype=pl.Utf8))
        .when(~pl.col("wpass")).then(pl.lit("below_min_amg_weight"))
        .when(pl.col("flag_active")).then(pl.lit("keyword_flag_filter"))
        .when(pl.col("in_array_ge3")).then(pl.lit("array_3plus"))
        .otherwise(pl.lit("kept")).alias("recon"))
    log_text = Path(sample_dir, "CheckAMG_annotate.log").read_text(encoding="utf-8", errors="replace")
    logged = {m.group(2): int(m.group(1).replace(",", "")) for m in re.finditer(r"Filtered ([\d,]+) AMGs for reason: (.+)", log_text)}
    counts = dict(ga.group_by("recon").len().rows())
    check = {"sample": sample, "recon_kept": counts.get("kept", 0), "curated": len(curated),
             "kept_set_equal": set(ga.filter(pl.col("recon") == "kept")["Protein"].to_list()) == curated}
    for text, reason in LOG_REASONS.items():
        check[f"log_{reason}"] = logged.get(text)
        check[f"recon_{reason}"] = counts.get(reason, 0)
    genes = ga.filter(pl.col("Protein").is_in(list(wanted))).join(fr, on="Protein", how="left").select(
        pl.lit(sample).alias("sample"), pl.col("Protein").alias("gene"), "no_hit", "cand", "recon",
        pl.col("Protein Classification").alias("checkamg_protein_classification"))
    return genes, check


def noncall_reason(called):
    """Reason CheckAMG did not call a gene, null where called is true, in the order CheckAMG annotate applies its filters."""
    import polars as pl

    return (pl.when(called).then(pl.lit(None, dtype=pl.Utf8))
            .when(pl.col("checkamg_protein_classification").is_null()).then(pl.lit("not_in_checkamg_input"))
            .when(pl.col("no_hit")).then(pl.lit("no_checkamg_hmm_hit"))
            .when(~pl.col("cand")).then(pl.lit("checkamg_hits_not_in_AMG_reference"))
            .when(pl.col("recon") == "below_min_amg_weight").then(pl.lit("in_AMG_reference_below_min_amg_weight"))
            .when(pl.col("recon") == "keyword_flag_filter").then(pl.lit("removed_by_keyword_flag_filter"))
            .when(pl.col("recon") == "array_3plus").then(pl.lit("removed_by_array_of_3plus_context_filter"))
            .otherwise(pl.lit("other_unexplained")))


def reconstruct_checkamg_filters(checkamg_dir, samples, files_dir, workers, threads_per_worker=1):
    """Filter outcome of the wanted genes per (subset, sample, wanted) job, rebuilt read-only from CheckAMG outputs and checked against its log."""
    import polars as pl

    reference = checkamg_reference(files_dir)
    jobs = [(str(Path(checkamg_dir, subset, sample)), sample, set(wanted), reference, threads_per_worker) for subset, sample, wanted in samples]
    with ProcessPoolExecutor(max_workers=max(1, min(workers, len(jobs))), mp_context=get_context("forkserver")) as pool:
        parts = list(pool.map(_reconstruct, jobs))
    return pl.concat([genes for genes, _ in parts]), pl.DataFrame([check for _, check in parts])
