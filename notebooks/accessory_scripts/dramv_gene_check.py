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


def compare_sample(mapping_path, sample, reformatted_faa, dramv_faa):
    pa.set_cpu_count(1)
    mapping = pq.read_table(mapping_path, columns=["old_gene", "new_gene"], filters=[("sample", "==", sample)])
    old_to_new = dict(zip(mapping.column("old_gene").to_pylist(), mapping.column("new_gene").to_pylist()))
    reformatted = dict(_sequence_digests(reformatted_faa))
    seen, unmapped, not_in_reformatted, differ = set(), 0, [], []
    for old_gene, digest in _sequence_digests(dramv_faa):
        # Genes without a GFF match kept their DRAM-V name in the reformatted file
        gene = old_to_new.get(old_gene)
        if gene is None:
            unmapped += 1
            gene = old_gene
        seen.add(gene)
        expected = reformatted.get(gene)
        if expected is None:
            not_in_reformatted.append(gene)
        elif expected != digest:
            differ.append(gene)
    not_in_dramv = [gene for gene in reformatted if gene not in seen]
    return {
        "sample": sample,
        "n_reformatted": len(reformatted),
        "n_dramv": len(seen),
        "n_unmapped_by_old_id": unmapped,
        "n_sequence_differs": len(differ),
        "n_not_in_reformatted": len(not_in_reformatted),
        "n_not_in_dramv": len(not_in_dramv),
        "examples_sequence_differs": differ[:N_EXAMPLES],
        "examples_not_in_reformatted": not_in_reformatted[:N_EXAMPLES],
        "examples_not_in_dramv": not_in_dramv[:N_EXAMPLES],
        "identical": not (differ or not_in_reformatted or not_in_dramv),
    }


def check_gene_set(mapping_path, genes_reformatted_dir, dramv_dir, cache_path, workers):
    """Per-sample comparison of DRAM-V genes.faa with the reformatted genes, cached by input fingerprints."""
    reformatted_samples = sorted((p.parent.parent.name, p.parent.name) for subset in BENCHMARK_SUBSETS
                                 for p in genes_reformatted_dir.joinpath(subset).glob("*/genes_reformatted.faa"))
    dramv_samples = sorted((p.parent.parent.name, p.parent.name) for subset in BENCHMARK_SUBSETS
                           for p in dramv_dir.joinpath(subset).glob("*/genes.faa"))
    if dramv_samples != reformatted_samples:
        raise RuntimeError(
            f"DRAM-V genes.faa is missing for {sorted(set(reformatted_samples) - set(dramv_samples))} and unexpected for "
            f"{sorted(set(dramv_samples) - set(reformatted_samples))}, finish the DRAM-V run first")
    reformatted_faa = [genes_reformatted_dir.joinpath(subset, sample, "genes_reformatted.faa") for subset, sample in dramv_samples]
    dramv_faa = [dramv_dir.joinpath(subset, sample, "genes.faa") for subset, sample in dramv_samples]
    key = [[str(p), p.stat().st_size, p.stat().st_mtime_ns] for p in [mapping_path, *reformatted_faa, *dramv_faa]]
    cached = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    if cached.get("key") == key:
        return dramv_samples, cached["results"]
    with ProcessPoolExecutor(max_workers=min(workers, len(dramv_samples)),
                             mp_context=multiprocessing.get_context("forkserver")) as pool:
        results = list(pool.map(compare_sample, [mapping_path] * len(dramv_samples),
                                [sample for _, sample in dramv_samples], reformatted_faa, dramv_faa))
    cache_path.write_text(json.dumps({"key": key, "results": results}, indent=1, default=str))
    return dramv_samples, results
