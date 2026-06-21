# CheckAMG/scripts/common/fasta_io.py

import os
import gzip
import shlex
from pathlib import Path
from pyfastatools import Parser

_NUCL_EXTS = (".fna", ".fasta", ".fa")
_PROT_EXTS = (".faa", ".fasta", ".fa")


class _SimpleHeader:
    __slots__ = ("name", "desc")

    def __init__(self, name, desc=""):
        self.name = name
        self.desc = desc


class _SimpleRecord:
    __slots__ = ("header", "seq")

    def __init__(self, name, seq, desc=""):
        self.header = _SimpleHeader(name, desc)
        self.seq = seq


def is_gz(path):
    return str(path).endswith(".gz")


def _strip_gz(path):
    s = str(path)
    return s[:-3] if s.endswith(".gz") else s


def is_nucl_fasta_path(path):
    s = _strip_gz(path)
    return any(s.endswith(e) for e in _NUCL_EXTS)


def is_prot_fasta_path(path):
    s = _strip_gz(path)
    return any(s.endswith(e) for e in _PROT_EXTS)


def is_strict_faa_path(path):
    return _strip_gz(path).endswith(".faa")


def _open_text(path):
    if is_gz(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def _iter_fasta_records_text(path):
    name = None
    desc = ""
    seq_chunks = []
    with _open_text(path) as fh:
        for line in fh:
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield _SimpleRecord(name, "".join(seq_chunks), desc)
                header = line[1:].rstrip("\n\r")
                if not header:
                    name = ""
                    desc = ""
                else:
                    parts = header.split(None, 1)
                    name = parts[0]
                    desc = parts[1] if len(parts) > 1 else ""
                seq_chunks = []
            else:
                seq_chunks.append(line.strip())
        if name is not None:
            yield _SimpleRecord(name, "".join(seq_chunks), desc)


def iter_fasta_records(path):
    """
    Iterate FASTA records (gz or plaintext). Records expose `.header.name`,
    `.header.desc`, and `.seq`.
    """
    if is_gz(path):
        return _iter_fasta_records_text(path)
    return Parser(path)


def iter_fasta_header_lines(path):
    """Yield raw header text (everything after `>`, no trailing newline). gz aware."""
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].rstrip("\n\r")


def write_fasta_custom(record, handle, small_line=75, big_chunk=100000):
    """
    Write a FASTA record. Wraps short sequences at `small_line`; very large
    sequences are written in `big_chunk` blocks to avoid excessive line counts.
    """
    seq = record.seq
    length = len(seq)
    name = str(record.header.name)
    desc = getattr(record.header, "desc", "") or ""
    handle.write(f">{name} {desc}\n" if desc else f">{name}\n")
    if length <= 1_000_000:
        for i in range(0, length, small_line):
            handle.write(seq[i:i + small_line] + "\n")
    else:
        for i in range(0, length, big_chunk):
            handle.write(seq[i:i + big_chunk] + "\n")


def output_basename_for_bin(bin_path, target_ext=None):
    """
    Build a non-gzipped output basename for a bin FASTA. If `target_ext` is
    given (e.g. ".fna"), replace `.fasta`/`.fa` with it; `.faa` is kept as-is
    unless explicitly retargeted.
    """
    out_base = os.path.basename(str(bin_path))
    if out_base.endswith(".gz"):
        out_base = out_base[:-3]
    if target_ext:
        for ext in (".fasta", ".fa"):
            if out_base.endswith(ext):
                out_base = out_base[:-len(ext)] + target_ext
                break
    return out_base


def _clean_path(p):
    if p is None:
        return None
    s = str(p) if not isinstance(p, Path) else str(p)
    s = os.path.expandvars(os.path.expanduser(s)).strip()
    return s or None


def list_fastas_in_dir(dir_path, predicate):
    if not os.path.isdir(dir_path):
        return []
    return sorted(
        os.path.join(dir_path, f)
        for f in os.listdir(dir_path)
        if predicate(f) and os.path.isfile(os.path.join(dir_path, f))
    )


_PRED_BY_KIND = {
    "nucl": is_nucl_fasta_path,
    "prot": is_prot_fasta_path,
    "faa": is_strict_faa_path,
}


def _expand_bin_input(bin_input, kind, logger, *, strict=True):
    """
    Resolve a bin input which may be:
      - a list/tuple of paths
      - a single directory (returns matching files inside)
      - a single FASTA file
      - a whitespace-separated string of paths/dirs
    Returns a deduplicated, ordered list of files matching `kind`.
    With `strict=False`, a single non-existent path silently yields [] (used
    for workflow-managed paths that may not have been written yet).
    """
    if kind not in _PRED_BY_KIND:
        raise ValueError(f"Unknown bin kind: {kind}")
    pred = _PRED_BY_KIND[kind]
    if not bin_input:
        return []

    if isinstance(bin_input, (list, tuple)):
        candidates = [_clean_path(x) for x in bin_input if x]
    else:
        s = _clean_path(bin_input)
        if s is None:
            return []
        logger.debug(f"bin_input ({kind}) raw: {s!r}")
        if os.path.isdir(s):
            return list_fastas_in_dir(s, pred)
        if os.path.isfile(s) and pred(s):
            return [s]
        if not any(ch.isspace() for ch in s):
            if strict:
                raise FileNotFoundError(
                    f"Bin input is not an existing directory or matching FASTA file ({kind}): {s!r}"
                )
            logger.debug(f"Bin input path missing or unmatched ({kind}): {s!r}")
            return []
        try:
            tokens = shlex.split(s)
        except ValueError:
            tokens = s.split()
        candidates = [_clean_path(x) for x in tokens]

    bin_files = []
    for c in candidates:
        if not c:
            continue
        if os.path.isdir(c):
            bin_files.extend(list_fastas_in_dir(c, pred))
        elif os.path.isfile(c) and pred(c):
            bin_files.append(c)
        else:
            if strict:
                raise FileNotFoundError(
                    f"Bin input is not an existing directory or matching FASTA file ({kind}): {c!r}"
                )
            logger.debug(f"Skipping unresolvable bin candidate: {c!r}")

    seen = set()
    return [p for p in bin_files if not (p in seen or seen.add(p))]


def _collect(single_path, bin_input, kind, logger, *, strict, label):
    single_clean = _clean_path(single_path)
    has_single = bool(single_clean) and os.path.isfile(single_clean)
    if single_clean and not has_single:
        if strict:
            raise FileNotFoundError(f"Single-contig {label} FASTA not found: {single_clean!r}")
        logger.debug(f"Skipping missing single-contig {label} FASTA: {single_clean!r}")

    bin_files = _expand_bin_input(bin_input, kind=kind, logger=logger, strict=strict)
    has_bin_files = len(bin_files) > 0

    if not has_single and not has_bin_files:
        raise FileNotFoundError(
            f"No valid input {label} FASTAs found (single contigs and bins both missing)."
        )

    input_files = ([single_clean] if has_single else []) + bin_files
    logger.debug(f"Detected single-contig {label} FASTA: {has_single}")
    logger.debug(f"Detected bin {label} FASTAs: {len(bin_files):,}")
    if bin_files:
        logger.debug(f"First 5 bin {label} FASTAs: " + ", ".join(bin_files[:5]))
    return input_files, bin_files, has_single, has_bin_files


def collect_nucleotide_inputs(input_single_contigs, input_bins, logger):
    """
    Resolve user-supplied nucleotide FASTAs (.fna/.fasta/.fa, optionally .gz).
    Bins may be passed as a directory, a list, or a whitespace-separated string.
    """
    return _collect(
        input_single_contigs, input_bins,
        kind="nucl", logger=logger, strict=True, label="nucleotide",
    )


def collect_protein_inputs(input_single_contig_prots, input_bin_prots, logger):
    """
    Resolve user-supplied protein FASTAs (.faa/.fasta/.fa, optionally .gz).
    Bins may be passed as a directory, a list, or a whitespace-separated string.
    """
    return _collect(
        input_single_contig_prots, input_bin_prots,
        kind="prot", logger=logger, strict=True, label="protein",
    )


def collect_pyrodigal_outputs(single_contig_prots_path, bin_proteins_dir, logger):
    """
    Collect pyrodigal-gv outputs: a single-contig .faa file (optional) and a
    directory of .faa bin files (optional). Workflow-managed paths can be
    absent when the upstream step decided not to write them, so missing paths
    are tolerated; only `.faa` (and `.faa.gz`) files are picked up.
    """
    return _collect(
        single_contig_prots_path, bin_proteins_dir,
        kind="faa", logger=logger, strict=False, label="protein",
    )
