from pathlib import Path
from typing import Iterator

import gc
import os
import shutil
import time
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)

import faiss
import torch
from pst.embed import main as esm2_embed
from pst.utils import graphify

from CheckAMG.scripts.denovo.constants import ESM2_EMBED_CHUNK_SIZE, ESM2_EMBED_RETRIES
from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: Compute ESM2 embeddings from input proteins", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")

import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="lightning")
warnings.filterwarnings("ignore", message=".*LeafSpec.*")

# Silence Lightning's promotional "Tip: ... install litlogger ..." message
# without touching its other rank_zero_info() output (Seed set, GPU available,
# TPU available, etc.). The tip is emitted by the trainer's logger_connector
# via logging.getLogger("lightning.pytorch.utilities.rank_zero").
import logging as _logging
class _DropLitLoggerTip(_logging.Filter):
    def filter(self, record):
        try:
            msg = record.getMessage()
        except Exception:
            return True
        return "litlogger" not in msg.lower()
_logging.getLogger("lightning.pytorch.utilities.rank_zero").addFilter(_DropLitLoggerTip())
_logging.getLogger("pytorch_lightning.utilities.rank_zero").addFilter(_DropLitLoggerTip())

# Per-chunk embedding would otherwise repeat Lightning's "Seed set to N",
# "GPU available", "TPU available", "Using bfloat16 AMP" and
# "LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [...]" banner for every chunk; silence
# those startup loggers so the consolidated per-chunk completion line below is
# the only output.
for _name in (
    "lightning.fabric.utilities.seed",
    "lightning.fabric.utilities.rank_zero",
    "lightning.pytorch.utilities.rank_zero",
    "lightning.pytorch.accelerators.cuda",
    "pytorch_lightning.utilities.rank_zero",
    "pytorch_lightning.accelerators.cuda",
):
    _logging.getLogger(_name).setLevel(_logging.WARNING)

# pst.embed.main.embed() constructs a fresh Lightning Trainer per chunk with
# hardcoded kwargs. To present a SINGLE progress bar covering all chunks, we
# register a shared Callback (see _GlobalEmbedProgressCallback below) and patch
# Trainer.__init__ to (a) suppress Lightning's default per-chunk progress bar
# and model summary, and (b) inject our shared callback into every Trainer.
import lightning as _L
_global_embed_pbar: dict = {"cb": None}

_orig_trainer_init = _L.Trainer.__init__
def _quiet_trainer_init(self, *args, **kwargs):
    kwargs.setdefault("enable_model_summary", False)
    cb = _global_embed_pbar.get("cb")
    if cb is not None:
        kwargs["enable_progress_bar"] = False
        existing = list(kwargs.get("callbacks") or [])
        if cb not in existing:
            existing.append(cb)
        kwargs["callbacks"] = existing
    return _orig_trainer_init(self, *args, **kwargs)
_L.Trainer.__init__ = _quiet_trainer_init

# Tensor Cores on recent NVIDIA GPUs (e.g. H100/H200) benefit from TF32 matmuls.
# Setting this also silences Lightning's float32_matmul_precision warning.
try:
    torch.set_float32_matmul_precision("high")
except Exception:
    pass

ESM2_CKPT_FILES = (
    "esm2_t30_150M_UR50D.pt",
    "esm2_t30_150M_UR50D-contact-regression.pt",
)


def _build_torch_hub_from_ckpt_dir(esm2_ckpt_dir: Path, work_dir: Path) -> Path:
    """Build a torch_hub directory whose 'checkpoints/' subdir contains the
    user-supplied ESM2 .pt files (symlinked). Returns the torch_hub path
    suitable for passing as ModelArgs.torch_hub.
    """
    missing = [f for f in ESM2_CKPT_FILES if not (esm2_ckpt_dir / f).is_file()]
    if missing:
        raise FileNotFoundError(
            f"Provided ESM2 checkpoint directory '{esm2_ckpt_dir}' is missing required files: {missing}"
        )
    torch_hub = work_dir / ".esm2_torch_hub"
    checkpoints = torch_hub / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    for fname in ESM2_CKPT_FILES:
        src = (esm2_ckpt_dir / fname).resolve()
        dst = checkpoints / fname
        if dst.is_symlink() or dst.exists():
            try:
                if dst.resolve() == src:
                    continue
            except OSError:
                pass
            dst.unlink()
        dst.symlink_to(src)
    return torch_hub


class _GlobalEmbedProgressCallback(_L.Callback):
    """A single tqdm bar shared across every chunk's Trainer.predict() call,
    tracking proteins embedded out of the total across all chunks.

    Between retries of the same chunk (e.g. after a transient CUDA fault) the
    bar is rolled back to a "committed" baseline so re-embedding a chunk does
    not double-count its proteins.
    """

    def __init__(self, total: int, initial: int = 0, desc: str = "ESM2 embedding"):
        super().__init__()
        self._total = total
        self._initial = initial
        self._committed = initial
        self._desc = desc
        self._bar = None

    def _ensure_bar(self):
        if self._bar is None:
            from tqdm.auto import tqdm
            self._bar = tqdm(
                total=self._total,
                initial=self._initial,
                desc=self._desc,
                unit=" proteins",
                unit_scale=False,
                smoothing=0.1,
                mininterval=1.0,
                dynamic_ncols=True,
            )

    def begin_chunk(self):
        """Mark the bar position at chunk entry as the rollback baseline."""
        self._ensure_bar()
        self._committed = self._bar.n

    def reset_to_committed(self):
        """Roll the bar back to the most recent committed baseline (used between
        retry attempts within a chunk)."""
        if self._bar is None:
            return
        self._bar.n = self._committed
        self._bar.refresh()

    def on_predict_batch_end(self, trainer, pl_module, outputs, batch, batch_idx, dataloader_idx=0):
        # SequenceDataset batches are (labels, seqs, tokens); len(labels) is the
        # number of proteins embedded in this step.
        self._ensure_bar()
        try:
            n = len(batch[0])
        except Exception:
            n = 1
        self._bar.update(n)

    def close(self):
        if self._bar is not None:
            self._bar.close()
            self._bar = None


def _build_model_cfg(esm2_ckpt_dir: Path | None, work_dir: Path):
    """Construct the ModelArgs for ESM2, pointing torch_hub at user-supplied
    checkpoints when provided, otherwise relying on (and warning about) the
    default torch hub cache."""
    if esm2_ckpt_dir:
        esm2_ckpt_dir = Path(esm2_ckpt_dir)
        torch_hub = _build_torch_hub_from_ckpt_dir(esm2_ckpt_dir, work_dir)
        logger.info(f"Using user-provided ESM2 checkpoints from {esm2_ckpt_dir} (torch_hub={torch_hub}).")
        return esm2_embed.ModelArgs(
            esm=esm2_embed.ESM2Models.esm2_t30_150M,
            torch_hub=torch_hub,
        )

    default_ckpts = [
        os.path.expanduser(f"~/.cache/torch/hub/checkpoints/{f}") for f in ESM2_CKPT_FILES
    ]
    if not all(os.path.exists(p) for p in default_ckpts):
        logger.warning(
            f"Required files {', '.join(ESM2_CKPT_FILES)} are not present in "
            "~/.cache/torch/hub/checkpoints. They will be downloaded, which may take some "
            "time on the first run. To avoid downloading, pass --esm2-ckpt-dir pointing to "
            "a directory containing these files."
        )
    else:
        logger.debug("Using ESM2 checkpoints already cached in ~/.cache/torch/hub/checkpoints.")
    return esm2_embed.ModelArgs(esm=esm2_embed.ESM2Models.esm2_t30_150M)


def _iter_fasta_records(fasta_file: Path) -> Iterator[tuple[str, list[str]]]:
    """Yield (header_line, sequence_lines) records, preserving raw line content
    so chunks are byte-identical slices of the input."""
    with fasta_file.open() as fp:
        header: str | None = None
        seq_lines: list[str] = []
        for line in fp:
            if line.startswith(">"):
                if header is not None:
                    yield header, seq_lines
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            yield header, seq_lines


def _scaffold_of(header: str) -> str:
    """Extract the scaffold id from a prodigal-style header '>scaffold_ptn # ...'."""
    ptn = header[1:].split()[0]
    return ptn.rsplit("_", 1)[0]


def _split_fasta_into_chunks(
    fasta_file: Path, chunks_dir: Path, chunk_size: int
) -> tuple[list[Path], list[int]]:
    """Split the FASTA into chunk files of up to ``chunk_size`` proteins each,
    never splitting a scaffold across chunks (proteins of a scaffold are
    contiguous in prodigal output, and graphify groups proteins by scaffold).

    Returns ``(chunk_paths, chunk_counts)`` so callers can compute total /
    resumed protein counts without re-reading the chunk files.

    Splitting is deterministic, so re-running overwrites chunk files identically
    and is safe alongside the resume logic in :func:`_embed_and_graphify_chunk`.
    """
    chunks_dir.mkdir(parents=True, exist_ok=True)

    chunk_paths: list[Path] = []
    chunk_counts: list[int] = []
    chunk_idx = 0
    proteins_in_chunk = 0
    current_scaffold: str | None = None
    out_fp = None

    try:
        for header, seq_lines in _iter_fasta_records(fasta_file):
            scaffold = _scaffold_of(header)
            # Only roll over to a new chunk at a scaffold boundary.
            start_new = out_fp is None or (
                proteins_in_chunk >= chunk_size and scaffold != current_scaffold
            )
            if start_new:
                if out_fp is not None:
                    out_fp.close()
                    chunk_counts.append(proteins_in_chunk)
                chunk_idx += 1
                path = chunks_dir / f"chunk_{chunk_idx:05d}.faa"
                out_fp = path.open("w")
                chunk_paths.append(path)
                proteins_in_chunk = 0

            out_fp.write(header)
            out_fp.writelines(seq_lines)
            proteins_in_chunk += 1
            current_scaffold = scaffold
    finally:
        if out_fp is not None:
            out_fp.close()
            chunk_counts.append(proteins_in_chunk)

    return chunk_paths, chunk_counts


def _try_clear_cuda():
    gc.collect()
    try:
        if torch.cuda.is_available():
            torch.cuda.synchronize()
            torch.cuda.empty_cache()
    except Exception:
        # After a hard CUDA fault the context may be unusable; ignore.
        pass


def _run_esm2_on_chunk(chunk_fasta: Path, model_cfg, devices: int, tmpdir: Path, attempts: int):
    """Embed a single chunk FASTA with ESM2, retrying on failure. Each attempt
    rebuilds the trainer/model (a fresh GPU context) and starts from a clean
    tmpdir so a partially-written 'predictions.h5' is never reused."""
    pbar = _global_embed_pbar.get("cb")
    if pbar is not None:
        pbar.begin_chunk()
    last_exc: Exception | None = None
    for attempt in range(1, attempts + 1):
        if tmpdir.exists():
            shutil.rmtree(tmpdir, ignore_errors=True)
        if pbar is not None and attempt > 1:
            pbar.reset_to_committed()
        try:
            esm2_embed.embed(
                input=chunk_fasta,
                outdir=tmpdir,
                model_cfg=model_cfg,
                trainer_cfg=esm2_embed.TrainerArgs(
                    devices=devices,
                    precision="bf16-mixed" if torch.cuda.is_available() else 32,
                ),
            )
            return
        except Exception as e:  # retry on any embedding failure (e.g. transient CUDA faults)
            last_exc = e
            logger.warning(
                f"ESM2 embedding attempt {attempt}/{attempts} for chunk {chunk_fasta.name} failed: {e}"
            )
            _try_clear_cuda()
            if attempt < attempts:
                time.sleep(5)
    raise RuntimeError(
        f"ESM2 embedding failed for chunk {chunk_fasta.name} after {attempts} attempt(s)"
    ) from last_exc


def _embed_and_graphify_chunk(
    idx: int,
    total: int,
    chunk_fasta: Path,
    model_cfg,
    devices: int,
    chunks_dir: Path,
    attempts: int,
) -> Path:
    """Embed + graphify one chunk, returning its graph-formatted h5. If that
    output already exists (a previous run completed this chunk), it is reused."""
    chunk_graphfmt = chunks_dir / f"{chunk_fasta.stem}.graphfmt.h5"
    if chunk_graphfmt.exists():
        logger.debug(f"Chunk {idx}/{total} ({chunk_fasta.name}) already embedded at {chunk_graphfmt}, resuming.")
        return chunk_graphfmt

    tmpdir = chunks_dir / f"tmp_esm2_embed_{chunk_fasta.stem}"
    logger.debug(f"Embedding chunk {idx}/{total} ({chunk_fasta.name}) using ESM2 with {devices} device(s)...")
    _run_esm2_on_chunk(chunk_fasta, model_cfg, devices, tmpdir, attempts)

    result_files = [p for p in tmpdir.glob("*.h5") if p.is_file()]
    if len(result_files) != 1:
        raise RuntimeError(
            f"Expected exactly one ESM2 result .h5 in {tmpdir} after embedding, found: {result_files}"
        )
    result_h5 = result_files[0]

    logger.debug(f"Graphifying chunk {idx}/{total} embeddings...")
    try:
        # to_graph_format validates the output and unlinks it on failure, so a
        # surviving chunk_graphfmt is always complete and valid.
        graphify.to_graph_format(
            io=graphify.IOArgs(
                file=result_h5,
                fasta_file=chunk_fasta,
                output=chunk_graphfmt,
            ),
            optional=graphify.OptionalArgs(),
        )
    except Exception as e:
        logger.error(f"Failed to graphify chunk {idx}/{total} ({chunk_fasta.name}): {e}")
        raise RuntimeError(f"Failed to graphify chunk {idx}/{total} ({chunk_fasta.name})") from e

    # The intermediate result .h5 lives inside tmpdir; drop the whole tmpdir.
    shutil.rmtree(tmpdir, ignore_errors=True)
    return chunk_graphfmt


def _fmt_duration(seconds: float) -> str:
    s = max(0, int(seconds))
    h, rem = divmod(s, 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _cleanup_chunks_dir(chunks_dir: Path):
    try:
        for p in chunks_dir.glob("chunk_*.faa"):
            p.unlink()
        for p in chunks_dir.glob("tmp_esm2_embed_*"):
            if p.is_dir():
                shutil.rmtree(p, ignore_errors=True)
        try:
            chunks_dir.rmdir()
        except OSError:
            pass
    except Exception as e:  # cleanup is best-effort
        logger.debug(f"Chunk directory cleanup encountered a non-fatal issue: {e}")


def _esm2_embed(
    fasta_file: Path,
    devices: int,
    outdir: Path,
    graphfmt_file: Path,
    esm2_ckpt_dir: Path | None = None,
    chunk_size: int = ESM2_EMBED_CHUNK_SIZE,
    attempts: int = ESM2_EMBED_RETRIES,
) -> Path:
    if graphfmt_file.exists():
        logger.info(
            f"Final graph-formatted ESM2 embeddings already exist at {graphfmt_file}, "
            "skipping embedding and graphification."
        )
        return graphfmt_file

    model_cfg = _build_model_cfg(esm2_ckpt_dir, outdir)

    chunks_dir = outdir / "esm2_chunks"
    logger.info(
        f"Splitting {fasta_file} into chunks of up to {chunk_size:,} proteins (scaffold-aware) "
        f"under {chunks_dir} ..."
    )
    chunk_fastas, chunk_counts = _split_fasta_into_chunks(fasta_file, chunks_dir, chunk_size)
    if not chunk_fastas:
        raise RuntimeError(f"No protein sequences found in {fasta_file}; nothing to embed.")
    total_chunks = len(chunk_fastas)
    total_proteins = sum(chunk_counts)
    resumed_proteins = sum(
        n for cf, n in zip(chunk_fastas, chunk_counts)
        if (chunks_dir / f"{cf.stem}.graphfmt.h5").exists()
    )
    already_done = sum(
        1 for cf in chunk_fastas if (chunks_dir / f"{cf.stem}.graphfmt.h5").exists()
    )
    if already_done:
        logger.info(
            f"Input split into {total_chunks} chunk(s) ({total_proteins:,} proteins); "
            f"resuming from prior run ({already_done}/{total_chunks} chunks already embedded)."
        )
    else:
        logger.info(f"Input split into {total_chunks} chunk(s) ({total_proteins:,} proteins).")

    # Register a single tqdm bar that spans every chunk's Trainer.predict call.
    # The Trainer monkey-patch picks it up via the _global_embed_pbar holder.
    _global_embed_pbar["cb"] = _GlobalEmbedProgressCallback(
        total=total_proteins,
        initial=resumed_proteins,
    )

    # The shared tqdm bar (set up just above) shows per-chunk progress; per-chunk
    # completion is logged at DEBUG only to keep the INFO stream clean.
    chunk_graphfmts: list[Path] = []
    try:
        for idx, chunk_fasta in enumerate(chunk_fastas, start=1):
            t0 = time.perf_counter()
            chunk_graphfmts.append(
                _embed_and_graphify_chunk(
                    idx, total_chunks, chunk_fasta, model_cfg, devices, chunks_dir, attempts
                )
            )
            logger.debug(
                f"Chunk {idx}/{total_chunks} done in {_fmt_duration(time.perf_counter() - t0)}."
            )
    finally:
        cb = _global_embed_pbar.get("cb")
        if cb is not None:
            cb.close()
        _global_embed_pbar["cb"] = None

    if len(chunk_graphfmts) == 1:
        # Single chunk: move directly to the final output (identical to the
        # pre-chunking single-file behavior).
        os.replace(chunk_graphfmts[0], graphfmt_file)
    else:
        logger.info(f"Merging {len(chunk_graphfmts)} chunk graph files into {graphfmt_file} ...")
        tmp_merged = graphfmt_file.with_name(graphfmt_file.name + ".merging.h5")
        if tmp_merged.exists():
            tmp_merged.unlink()
        # merge_graph_files validates the merged output and unlinks the inputs
        # on success; write to a temp path then atomically rename so a partial
        # merge is never mistaken for a finished output.
        graphify.merge_graph_files(chunk_graphfmts, tmp_merged)
        os.replace(tmp_merged, graphfmt_file)

    _cleanup_chunks_dir(chunks_dir)
    logger.info(f"Graphified ESM2 embeddings saved to {graphfmt_file}.")
    return graphfmt_file


def _handle_cpu_accelerator(accelerator: str | None, devices: int):
    if accelerator == "cpu" or (
        accelerator == "auto" and not torch.cuda.is_available()
    ):
        torch.set_num_threads(devices)
        faiss.omp_set_num_threads(devices)

def main():
    logger.info("ESM2 embedding and graphification starting...")
    accelerator = str(snakemake.params.accelerator)
    devices = int(snakemake.params.devices)
    logger.debug(f"Accelerator={accelerator}, devices={devices}")
    _handle_cpu_accelerator(accelerator, devices)
    
    fasta_file = Path(snakemake.params.query_fasta_file)
    output_dir = Path(snakemake.params.output_dir)
    output_esm2_graph_file = Path(snakemake.params.esm2_graph_file)
    esm2_ckpt_dir_param = getattr(snakemake.params, "esm2_ckpt_dir", None)
    esm2_ckpt_dir = Path(esm2_ckpt_dir_param) if esm2_ckpt_dir_param else None
    chunk_size_param = getattr(snakemake.params, "esm2_chunk_size", None)
    chunk_size = int(chunk_size_param) if chunk_size_param else ESM2_EMBED_CHUNK_SIZE

    logger.debug(
        f"Input FASTA file: {fasta_file}, Output directory: {output_dir}, "
        f"Output ESM2 graph file: {output_esm2_graph_file}, ESM2 checkpoint dir: {esm2_ckpt_dir}, "
        f"Chunk size: {chunk_size}"
    )
    _esm2_embed(
        fasta_file,
        devices,
        output_dir,
        output_esm2_graph_file,
        esm2_ckpt_dir=esm2_ckpt_dir,
        chunk_size=chunk_size,
    )

    logger.info("ESM2 embedding and graphification completed.")

if __name__ == "__main__":
    main()
