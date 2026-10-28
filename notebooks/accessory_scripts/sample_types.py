import polars as pl

GENE_KEY = ["gene", "sample", "source", "ecosystem"]
BENCHMARK_SUBSETS = ("complete_virus_genomes", "metagenomes", "viromes")
N_BENCHMARK_SAMPLES = 21


def benchmark_gene_id_mapping(path, columns=None):
    """gene_id_mapping.parquet restricted to the benchmark samples, the subsets with VirSorter2 outputs."""
    mapping = pl.read_parquet(path, columns=columns).filter(pl.col("subset").is_in(BENCHMARK_SUBSETS))
    n_samples = mapping["sample"].n_unique()
    assert n_samples == N_BENCHMARK_SAMPLES, f"{n_samples} benchmark samples in {path}, expected {N_BENCHMARK_SAMPLES}"
    return mapping


def assert_benchmark_samples(frame, samples, label):
    outside = _collect(frame.select("sample").unique().filter(~pl.col("sample").is_in(list(samples))))
    assert outside.height == 0, f"{label}: genes of samples outside the benchmark: {sorted(outside['sample'].to_list())}"


# Sample type by sample name, overriding the directory placement of these six soil samples
VIROMES = ("soil_T42_15_1_55", "soil_T42_15_2_56", "soil_T42_15_3_57")
MIXED_METAGENOMES = ("soil_T42_15_1_43", "soil_T42_15_2_44", "soil_T42_15_3_45")


def sample_type_corrections():
    return {**{sample: "Viromes" for sample in VIROMES},
            **{sample: "Mixed metagenomes" for sample in MIXED_METAGENOMES}}


def apply_sample_type_corrections(frame, corrections):
    """Replace the display sample type of the corrected samples, keyed by sample name rather than directory placement."""
    return frame.with_columns(
        pl.col("sample").replace_strict(corrections, default=None, return_dtype=pl.Utf8).fill_null(pl.col("source")).alias("source"))


def _collect(frame):
    return frame.collect() if isinstance(frame, pl.LazyFrame) else frame


def assert_sample_types(frame, corrections, label):
    pairs = _collect(frame.filter(pl.col("sample").is_in(list(corrections))).select("sample", "source").unique())
    assert pairs.height > 0, f"{label}: no rows for the corrected samples"
    wrong = pairs.filter(pl.col("source") != pl.col("sample").replace_strict(corrections, default=None, return_dtype=pl.Utf8))
    assert wrong.height == 0, f"{label}: as-labeled sample types remain for {sorted(wrong.rows())}"


def assert_key_join_complete(left, right, key, label, corrections):
    """Every corrected-sample row of left that right holds on the key without source also matches on the full key."""
    loose = [c for c in key if c != "source"]
    soil = _collect(left.filter(pl.col("sample").is_in(list(corrections))).select(key).unique())
    right_keys = _collect(right.filter(pl.col("sample").is_in(list(corrections))).select(key).unique())
    assert soil.height > 0, f"{label}: no rows for the corrected samples"
    n_full = soil.join(right_keys, on=key, how="semi").height
    n_loose = soil.join(right_keys.select(loose).unique(), on=loose, how="semi").height
    assert n_full == n_loose, f"{label}: {n_loose - n_full} corrected-sample keys match on {loose} but not on {key}"
    return n_full
