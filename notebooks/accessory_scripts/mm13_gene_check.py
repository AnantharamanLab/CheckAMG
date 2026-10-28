import hashlib
import json
import multiprocessing
from concurrent.futures import ProcessPoolExecutor

import pyarrow as pa
import pyarrow.parquet as pq
from Bio.SeqIO.FastaIO import SimpleFastaParser

N_EXAMPLES = 5
BENCHMARK_SUBSETS = ("complete_virus_genomes", "metagenomes", "viromes")


def _sequence_digests(path):
    with open(path) as handle:
        for title, sequence in SimpleFastaParser(handle):
            # Line wrapping and a terminal stop differ between the DRAM-V and reformatted writers
            yield title.split(None, 1)[0], hashlib.md5(sequence.rstrip("*").upper().encode()).digest()


def compare_sample(mapping_path, sample, published_faa, mm13_faa):
    pa.set_cpu_count(1)
    mapping = pq.read_table(mapping_path, columns=["old_gene", "new_gene"], filters=[("sample", "==", sample)])
    old_to_new = dict(zip(mapping.column("old_gene").to_pylist(), mapping.column("new_gene").to_pylist()))
    published = dict(_sequence_digests(published_faa))
    seen, unmapped, not_in_published, differ = set(), 0, [], []
    for old_gene, digest in _sequence_digests(mm13_faa):
        # Genes without a GFF match kept their DRAM-V name in the published reformatted file
        gene = old_to_new.get(old_gene)
        if gene is None:
            unmapped += 1
            gene = old_gene
        seen.add(gene)
        expected = published.get(gene)
        if expected is None:
            not_in_published.append(gene)
        elif expected != digest:
            differ.append(gene)
    not_in_mm13 = [gene for gene in published if gene not in seen]
    return {
        "sample": sample,
        "n_published": len(published),
        "n_mm13": len(seen),
        "n_unmapped_by_old_id": unmapped,
        "n_sequence_differs": len(differ),
        "n_not_in_published": len(not_in_published),
        "n_not_in_mm13": len(not_in_mm13),
        "examples_sequence_differs": differ[:N_EXAMPLES],
        "examples_not_in_published": not_in_published[:N_EXAMPLES],
        "examples_not_in_mm13": not_in_mm13[:N_EXAMPLES],
        "identical": not (differ or not_in_published or not_in_mm13),
    }


def check_gene_set(mapping_path, genes_reformatted_dir, mm13_dir, cache_path, workers):
    """Per-sample comparison of mm13 genes.faa with the published reformatted genes, cached by input fingerprints."""
    published_samples = sorted((p.parent.parent.name, p.parent.name) for subset in BENCHMARK_SUBSETS
                               for p in genes_reformatted_dir.joinpath(subset).glob("*/genes_reformatted.faa"))
    mm13_samples = sorted((p.parent.parent.name, p.parent.name) for subset in BENCHMARK_SUBSETS
                          for p in mm13_dir.joinpath(subset).glob("*/genes.faa"))
    if mm13_samples != published_samples:
        raise RuntimeError(
            f"mm13 genes.faa is missing for {sorted(set(published_samples) - set(mm13_samples))} and unexpected for "
            f"{sorted(set(mm13_samples) - set(published_samples))}, finish the mm13 re-run first")
    published_faa = [genes_reformatted_dir.joinpath(subset, sample, "genes_reformatted.faa") for subset, sample in mm13_samples]
    mm13_faa = [mm13_dir.joinpath(subset, sample, "genes.faa") for subset, sample in mm13_samples]
    key = [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in [mapping_path, *published_faa, *mm13_faa]]
    cached = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if cached.get("key") == key:
        return mm13_samples, cached["results"]
    with ProcessPoolExecutor(max_workers=min(workers, len(mm13_samples)),
                             mp_context=multiprocessing.get_context("forkserver")) as pool:
        results = list(pool.map(compare_sample, [mapping_path] * len(mm13_samples),
                                [sample for _, sample in mm13_samples], published_faa, mm13_faa))
    cache_path.write_text(json.dumps({"key": key, "results": results}, indent=1, default=str))
    return mm13_samples, results
