import os
import re
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import get_context
from pathlib import Path

CAT_SUFFIX = re.compile(r"-cat_\d+$")


def fasta_lengths(path):
    lengths, name, n = {}, None, 0
    with open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                if name is not None:
                    lengths[name] = n
                name, n = CAT_SUFFIX.sub("", line[1:].split()[0].decode()), 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def for_dramv_fasta(virsorter_dir, subset, sample):
    found = sorted(Path(virsorter_dir, subset, sample, "for-dramv").glob("*.fa"))
    if len(found) != 1:
        raise FileNotFoundError(f"expected one for-dramv FASTA for {subset}/{sample}, found {len(found)}")
    return found[0]


def _sample_lengths(job):
    virsorter_dir, subset, sample = job
    lengths = fasta_lengths(for_dramv_fasta(virsorter_dir, subset, sample))
    return subset, sample, list(lengths), list(lengths.values())


def true_contig_lengths(virsorter_dir, pairs, workers):
    """Sequence length of every VirSorter2 for-dramv record, keyed by subset, sample and header minus -cat_N."""
    import polars as pl

    jobs = [(str(virsorter_dir), subset, sample) for subset, sample in pairs]
    os.environ.setdefault("POLARS_MAX_THREADS", "1")
    with ProcessPoolExecutor(max_workers=max(1, min(workers, len(jobs))), mp_context=get_context("forkserver")) as pool:
        parts = list(pool.map(_sample_lengths, jobs))
    return pl.concat([
        pl.DataFrame({"subset": subset, "sample": sample, "scaffold": names, "true_len": lens},
                     schema={"subset": pl.Utf8, "sample": pl.Utf8, "scaffold": pl.Utf8, "true_len": pl.Int64})
        for subset, sample, names, lens in parts
    ])
