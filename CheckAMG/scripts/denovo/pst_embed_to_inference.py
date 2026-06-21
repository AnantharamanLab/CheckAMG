import math
from pathlib import Path
from typing import NamedTuple

import os
os.environ["NUMEXPR_MAX_THREADS"] = str(snakemake.threads)

import faiss
import numpy as np
import tables as tb
import torch
from pst.typing import OptTensor
from pyfastatools import Parser
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm

from CheckAMG.scripts.denovo.constants import (
    EMBED_FIELD,
    LABEL_FIELD,
    NUM_CLASSES,
)
from CheckAMG.scripts.denovo.data import GenomeDataModule
from CheckAMG.scripts.denovo.embed import embed
from CheckAMG.scripts.denovo.model import CheckAMGPST
from CheckAMG.scripts.denovo.utils import load_model

from CheckAMG.scripts.common.snakemake_logging import init_snakemake_logger, emit_step_banner
from CheckAMG.scripts.common.mem_limit import set_memory_limit

log_file = snakemake.params.log
logger = init_snakemake_logger(snakemake, redirect_streams=False)
emit_step_banner(title=f"Step {snakemake.params.step[0]}/{snakemake.params.step[1]}: PST model inference", log_path=log_file, append_raw=False)
set_memory_limit(snakemake.resources.mem)
logger.debug(f"Memory limit set to {snakemake.resources.mem} GB")


class Decomposed(NamedTuple):
    avg_like: torch.Tensor
    viral: torch.Tensor

    def __iter__(self):
        return zip(self.avg_like, self.viral)


# Only sets the progress-bar granularity for the add phase; chunked add is byte-identical to a single add, but below ~50k per-call overhead grows sharply.
INDEX_ADD_CHUNK_SIZE = 50_000


def build_index(
    data: torch.Tensor, n_cells: int | None = None, n_threads: int | None = None
) -> faiss.IndexIVFFlat:
    # CPU mode passes the job's allotted thread budget; GPU mode passes None to keep FAISS at all cores.
    if n_threads is not None:
        n_threads = max(1, int(n_threads))
        faiss.omp_set_num_threads(n_threads)

    n, dim = data.shape

    max_cells = max(n // 39, 1)
    if n_cells is None:
        sqrt_cells = int(math.sqrt(n)) * 4

        # for small datasets, max_cells is smaller than sqrt_cells
        n_cells = min(max_cells, sqrt_cells)
    else:
        n_cells = min(n_cells, max_cells)

    quantizer = faiss.IndexFlatL2(dim)
    index = faiss.IndexIVFFlat(quantizer, dim, n_cells)

    logger.debug(f"Building FAISS index over {n} vectors using up to {faiss.omp_get_max_threads()} threads")
    # faiss-cpu bundles a single-threaded BLAS, so its default sgemm path runs k-means assignment on one core. Forcing the non-BLAS path parallelizes assignment over FAISS's OpenMP pool and is deterministic across thread counts.
    prev_blas_threshold = faiss.cvar.distance_compute_blas_threshold
    faiss.cvar.distance_compute_blas_threshold = n + 1
    try:
        cp = index.cp
        # Subsample the training set once, replicating FAISS's seeded rand_perm, so the per-iteration loop trains on the same points as a one-shot index.train() without re-subsampling (and re-copying) every iteration.
        n_sub = n_cells * cp.max_points_per_centroid
        if n_sub < n:
            perm = np.empty(n, dtype=np.int32)
            faiss.rand_perm(faiss.swig_ptr(perm), n, cp.seed)
            train_vectors = data[torch.from_numpy(perm[:n_sub].astype(np.int64))].contiguous()
        else:
            train_vectors = data

        # Run k-means one iteration at a time so the bar advances per iteration; byte-identical to a single index.train() (seeded init + warm start reproduce the internal loop).
        clustering = faiss.Clustering(dim, n_cells, cp)
        clustering.niter = 1
        quantizer.reset()
        for _ in tqdm(range(cp.niter), desc="Clustering training vectors", unit="iter"):
            clustering.train(train_vectors, quantizer)  # type: ignore
        index.is_trained = True

        for start in tqdm(
            range(0, n, INDEX_ADD_CHUNK_SIZE),
            desc="Adding training vectors to index",
            unit="chunk",
        ):
            index.add(data[start : start + INDEX_ADD_CHUNK_SIZE])  # type: ignore
    finally:
        faiss.cvar.distance_compute_blas_threshold = prev_blas_threshold
    return index


def collate(batchlist: list[tuple[torch.Tensor]]) -> torch.Tensor:
    return torch.stack([x[0] for x in batchlist])


def knn_search(
    data: torch.Tensor,
    index: faiss.Index,
    k: int,
    batch_size: int = 4096,
    verbose: bool = True,
) -> tuple[torch.Tensor, torch.Tensor]:
    # index is for the training data
    # while data is any arbitrary querying data, ie may not be the training data

    dataset = TensorDataset(data)
    dataloader = DataLoader(
        dataset=dataset, batch_size=batch_size, shuffle=False, collate_fn=collate
    )

    # add 1 to ignore self search
    n = data.shape[0]

    dist = torch.empty((n, k), device=data.device)
    nn_idx = torch.full((n, k), -1, dtype=torch.long, device=data.device)

    batch: torch.Tensor
    for batch_idx, batch in enumerate(
        tqdm(dataloader, disable=not verbose, desc="Computing nearest neighbors")
    ):
        start = batch_idx * batch_size
        stop = start + batch_size

        local_dist, local_nn = index.search(batch, k=k)

        dist[start:stop] = torch.from_numpy(local_dist)
        nn_idx[start:stop] = torch.from_numpy(local_nn)

    dist.masked_fill_(nn_idx == -1, torch.inf)

    return dist, nn_idx


def compute_neighbor_weights(dist: torch.Tensor) -> torch.Tensor:
    std = dist[dist.isfinite()].std()

    # shape: [n, num_neighbors]
    weights = torch.softmax(-dist / std, dim=-1)

    return weights


def bincount_2d(
    input: torch.Tensor, weights: OptTensor = None, minlength: int = 0
) -> torch.Tensor:
    """
    Compute bincount for a 2D long tensor with optional weights.

    This is equivalent to applying torch.bincount to each row independently,
    where each row has its own corresponding weights.

    Args:
        input: 2D tensor of shape (batch_size, seq_len) with non-negative integers
        weights: Optional 2D tensor of same shape as input containing weights
        minlength: Minimum number of bins (default: 0)

    Returns:
        2D tensor of shape (batch_size, num_bins) where num_bins is max(input.max() + 1, minlength)

    Example:
        >>> input = torch.tensor([[0, 1, 1, 2], [0, 0, 1, 3]])
        >>> weights = torch.tensor([[1.0, 2.0, 3.0, 4.0], [0.5, 1.5, 2.5, 3.5]])
        >>> bincount_2d(input, weights)
        tensor([[1.0, 5.0, 4.0, 0.0],
                [2.0, 2.5, 0.0, 3.5]])
    """
    assert input.dim() == 2, "Input must be 2D"
    assert input.dtype in [torch.long, torch.int, torch.int32, torch.int64], (
        "Input must be integer type"
    )

    if weights is not None:
        assert weights.shape == input.shape, "Weights must have same shape as input"

    batch_size = input.shape[0]

    # Determine number of bins
    if input.numel() > 0:
        num_bins = max(int(input.max()) + 1, minlength)
    else:
        num_bins = minlength

    # Initialize output
    if weights is not None:
        output = torch.zeros(
            batch_size, num_bins, dtype=weights.dtype, device=input.device
        )
    else:
        output = torch.zeros(
            batch_size, num_bins, dtype=torch.long, device=input.device
        )

    # Create batch indices for scatter_add
    batch_indices = (
        torch.arange(batch_size, device=input.device).unsqueeze(1).expand_as(input)
    )

    # Flatten for scatter_add
    flat_batch_indices = batch_indices.reshape(-1)
    flat_input = input.reshape(-1)

    if weights is not None:
        flat_weights = weights.reshape(-1)
    else:
        flat_weights = torch.ones_like(flat_input, dtype=torch.long)

    # Use index_put with accumulate=True (equivalent to scatter_add for 2D indexing)
    output.index_put_(
        (flat_batch_indices, flat_input),
        flat_weights,
        accumulate=True,
    )

    return output


def weighted_vote(nn_y: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
    class_probs = bincount_2d(input=nn_y, weights=weights, minlength=NUM_CLASSES)

    # uniform prob for data points without any neighbors
    class_probs[torch.all(nn_y == -1, -1)] = 1 / NUM_CLASSES

    return class_probs


def decompose_class_probs(class_probs: torch.Tensor) -> Decomposed:
    # class_probs: [n, 4 classes]
    # 0 = non AMG-like / nonviral
    # 1 = non AMG-like / viral
    # 2 = AMG-like / non viral
    # 3 = AMG-like / viral

    # so we just need to return the sum of positive probs
    # for the 2 binary labels

    AVG_LIKE_IDX = [2, 3]
    VIRAL_IDX = [1, 3]

    viral_pos_prob = class_probs[:, VIRAL_IDX].sum(-1)
    avg_pos_prob = class_probs[:, AVG_LIKE_IDX].sum(-1)

    return Decomposed(avg_pos_prob, viral_pos_prob)


OptPath = None | Path


def _read_labels(file: Path) -> torch.Tensor:
    with tb.open_file(file.as_posix()) as fp:
        return torch.from_numpy(fp.root[LABEL_FIELD][:])  # type: ignore


def _read_embeddings(file: Path) -> torch.Tensor:
    with tb.open_file(file.as_posix()) as fp:
        return torch.from_numpy(fp.root[EMBED_FIELD][:])  # type: ignore


# TODO: optional move index to gpu
def handle_training_inputs(
    model: CheckAMGPST,
    train_data_file: OptPath,
    train_embed_file: OptPath,
    train_index_file: OptPath,
    train_labels_file: OptPath,
    batch_size: int,
    n_cells: int | None,
    output_dir: Path,
    n_threads: int | None = None,
) -> tuple[faiss.IndexIVFFlat, torch.Tensor]:
    
    if train_index_file is not None:
        logger.info(f"Using precomputed training index from {train_index_file}...")
        if train_labels_file is None:
            logger.error(
                "Training labels file not found at {train_labels_file} while training index file is provided at {train_index_file}"
                "BOTH the training index file and the training labels file must be provided if one is"
                )
            raise ValueError(
                "BOTH the training index file and the training labels file must be provided if one is."
            )
        index = faiss.read_index(str(train_index_file))
        logger.info(f"Using training labels from {train_labels_file}...")
        y = _read_labels(train_labels_file)

    elif train_embed_file is not None:
        logger.debug(f"No precomputed training index found.")
        logger.info(f"Using precomputed training embeddings from {train_embed_file}...")
        logger.info(f"Building index from training embeddings at {train_embed_file}...")
        train_embed = _read_embeddings(train_embed_file)

        logger.info("Building train data search index...")
        index = build_index(train_embed, n_cells=n_cells, n_threads=n_threads)
        # Write the index to disk so it can be reused for inference without needing to recompute the embeddings
        train_index_path = train_embed_file.with_suffix(".index.faiss")
        logger.info(f"Saving training index to {train_index_path}...")
        faiss.write_index(index, train_index_path.as_posix())

        # TODO: should maybe be in the decomposed format?
        logger.info(f"Reading training labels from {train_embed_file}...")
        y = _read_labels(train_embed_file)
        # Write the labels to disk so they can be reused for inference without needing to recompute the embeddings and re-extract the labels from the datamodule
        train_labels_path = train_embed_file.with_suffix(".labels.h5")
        logger.info(f"Saving training labels to {train_labels_path}...")
        with tb.open_file(train_labels_path.as_posix(), mode="w") as fp:
            fp.create_carray(fp.root, LABEL_FIELD, obj=y.numpy())

    elif train_data_file is not None:
        logger.warning(
            f"No precomputed training index or training embeddings found. PST-formatted embedding will be calculated for the training data. This will take some time..."
        )
        datamodule = GenomeDataModule(train_file=train_data_file)

        datamodule.setup()

        train_embed_file = output_dir / f"{train_data_file.stem}.PST-EMBED.h5"
        logger.info(f"Embedding training data -> {train_embed_file}")
        embed(
            model=model,
            dataset=datamodule.train_dataset,
            output=train_embed_file,
            batch_size=batch_size,
        )
        train_embed = _read_embeddings(train_embed_file)

        # Write the index to disk so it can be reused for inference without needing to recompute the embeddings
        logger.info("Building train data search index...")
        index = build_index(train_embed, n_cells=n_cells, n_threads=n_threads)
        train_index_path = train_embed_file.with_suffix(".index.faiss")
        # index = faiss.read_index(str(train_index_path)) # debugging
        logger.info(f"Saving training index to {train_index_path}...")
        faiss.write_index(index, train_index_path.as_posix())

        # Write the labels to disk so they can be reused for inference without needing to recompute the embeddings and re-extract the labels from the datamodule
        logger.info(f"Reading training labels from {train_embed_file}...")
        y = _read_labels(train_embed_file)
        # from CheckAMG.scripts.denovo.utils import _extract_label_from_datamodule # debugging
        # y = _extract_label_from_datamodule(datamodule) # debugging
        # with tb.open_file(train_embed_file, mode="a") as fp: # debugging
        #     atom = tb.Atom.from_dtype(y.numpy().dtype) # debugging
        #     ds = fp.create_carray(fp.root, "label", atom, y.shape) # debugging
        #     ds[:] = y.numpy() # debugging
        train_labels_path = train_embed_file.with_suffix(".labels.h5")
        logger.info(f"Saving training labels to {train_labels_path}...")
        with tb.open_file(train_labels_path.as_posix(), mode="w") as fp:
            fp.create_carray(fp.root, LABEL_FIELD, obj=y.numpy())

    return index, y


def _contig_from_protein(name: str) -> str:
    """Recover the contig a protein is encoded on from its name.

    Pyrodigal names proteins "{contig}_{n}". Contigs with too many CDS are
    chunked by filter_prots.py into "{contig}_SPLIT_{i}_{n}", so strip the
    split marker first to recover the original contig.
    """
    if "_SPLIT_" in name:
        return name.split("_SPLIT_")[0]
    return name.rsplit("_", 1)[0]


def _load_contig_to_genome(file: Path) -> dict[str, str]:
    """Read a Contig<TAB>Genome table (with header) into a dict."""
    mapping = {}
    with open(file) as fh:
        fh.readline() # skip header
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                mapping[parts[0]] = parts[1]
    return mapping


def _confidence_level(prob: float, very_high: float, high: float, medium: float) -> str:
    """Map a probability to a confidence level using the configured cutoffs.

    At or above `very_high` is "very high", at or above `high` is "high", at or
    above `medium` is "medium", and anything below `medium` is "low".
    """
    if prob >= very_high:
        return "very high"
    if prob >= high:
        return "high"
    if prob >= medium:
        return "medium"
    return "low"


def _is_cpu(accelerator: str | None) -> bool:
    return accelerator == "cpu" or (
        accelerator == "auto" and not torch.cuda.is_available()
    )


def _handle_cpu_accelerator(accelerator: str | None, devices: int):
    if _is_cpu(accelerator):
        torch.set_num_threads(devices)
        faiss.omp_set_num_threads(devices)


def main():
    logger.info("PST model inference starting...")

    accelerator = str(snakemake.params.accelerator)
    devices = int(snakemake.params.devices)
    num_probe_cells = snakemake.params.num_probe_cells
    query_file_pst = Path(snakemake.params.query_file_pst) if snakemake.params.query_file_pst else None
    train_data_file = Path(snakemake.params.train_data_file) if snakemake.params.train_data_file else None
    train_embed_file = Path(snakemake.params.train_embed_file) if snakemake.params.train_embed_file else None
    train_index_file = Path(snakemake.params.train_index_file) if snakemake.params.train_index_file else None
    train_labels_file = Path(snakemake.params.train_labels_file) if snakemake.params.train_labels_file else None
    num_index_cells = snakemake.params.num_index_cells
    knn_k = snakemake.params.knn_k
    batch_size = snakemake.params.batch_size
    nn_batch_size = snakemake.params.nn_batch_size
    fasta_file = Path(snakemake.params.query_fasta_file) if snakemake.params.query_fasta_file else None
    output_file = Path(snakemake.params.output_file)
    contig_to_genome_file = Path(snakemake.params.contig_to_genome_file) if snakemake.params.contig_to_genome_file else None
    pst_thresholds = snakemake.params.pst_thresholds
    avg_very_high = pst_thresholds["AVG-like"]["very_high"]
    avg_high = pst_thresholds["AVG-like"]["high"]
    avg_medium = pst_thresholds["AVG-like"]["medium"]
    viral_very_high = pst_thresholds["Viral"]["very_high"]
    viral_high = pst_thresholds["Viral"]["high"]
    viral_medium = pst_thresholds["Viral"]["medium"]
    # Final AVG (AVG-like AND viral) thresholds, applied to the product of the two marginals.
    final_avg_very_high = pst_thresholds["AVG"]["very_high"]
    final_avg_high = pst_thresholds["AVG"]["high"]
    final_avg_medium = pst_thresholds["AVG"]["medium"]
    logger.debug(
        f"PST confidence cutoffs - AVG-like: very_high>={avg_very_high}, high>={avg_high}, medium>={avg_medium}; "
        f"Viral: very_high>={viral_very_high}, high>={viral_high}, medium>={viral_medium}; "
        f"Final AVG: very_high>={final_avg_very_high}, high>={final_avg_high}, medium>={final_avg_medium}"
    )

    logger.debug(f"Accelerator={accelerator}, devices={devices}")
    _handle_cpu_accelerator(accelerator, devices)

    model_ckpt = snakemake.params.model_ckpt
    logger.debug(f"Model checkpoint for PST embedding: {model_ckpt}")
    model = load_model(ckpt=model_ckpt, device=accelerator)

    logger.info("Handling model training inputs...")
    # On CPU, bound FAISS index building to the job's allotted threads; on GPU, leave FAISS free to use all cores.
    index_build_threads = snakemake.threads if _is_cpu(accelerator) else None
    index, train_y = handle_training_inputs(
        model=model,
        train_data_file=train_data_file,
        train_embed_file=train_embed_file,
        train_index_file=train_index_file,
        train_labels_file=train_labels_file,
        batch_size=batch_size,
        n_cells=num_index_cells,
        output_dir=output_file.parent,
        n_threads=index_build_threads,
    )
    
    if num_probe_cells is not None:
        index.nprobe = num_probe_cells

    logger.info("Handling query inputs...")
    query_data = _read_embeddings(query_file_pst)

    dist, nn_idx = knn_search(
        data=query_data, index=index, k=knn_k, batch_size=nn_batch_size
    )
    weights = compute_neighbor_weights(dist)

    nn_y = train_y[nn_idx]
    class_probs = decompose_class_probs(weighted_vote(nn_y, weights))

    if fasta_file is not None:
        seq_iter = (header.name for header in Parser(fasta_file).headers())
    else:
        seq_iter = (f"seq_{i}" for i in range(query_data.shape[0]))

    # The contig-to-genome bridge is only produced by the filter step (fna/faa
    # input modes). When it is absent (esm2/pst_embed modes) the genome falls
    # back to the contig name.
    contig_to_genome = {}
    if contig_to_genome_file is not None and contig_to_genome_file.exists():
        logger.info(f"Loading contig-to-genome mapping from {contig_to_genome_file}...")
        contig_to_genome = _load_contig_to_genome(contig_to_genome_file)
    else:
        logger.debug("No contig-to-genome mapping found; genome defaults to contig name.")

    logger.info(f"Writing PST inference results to {output_file}...")
    with open(output_file, "w") as fp:
        fp.write("Genome\tContig\tProtein\tAVG-like prob\tAVG-like confidence\tViral prob\tViral confidence\tFinal AVG prob\tFinal AVG confidence\n")

        for seqname, (avg_like_prob, viral_prob) in zip(seq_iter, class_probs):
            avg_like_prob = float(avg_like_prob)
            viral_prob = float(viral_prob)
            contig = _contig_from_protein(seqname)
            genome = contig_to_genome.get(contig, contig)
            avg_like_conf = _confidence_level(avg_like_prob, avg_very_high, avg_high, avg_medium)
            viral_conf = _confidence_level(viral_prob, viral_very_high, viral_high, viral_medium)
            # Final AVG = AVG-like AND viral. Probability is the product of the two marginal
            # probabilities; confidence uses the dedicated AND-gate (final AVG) thresholds.
            final_avg_prob = avg_like_prob * viral_prob
            final_avg_conf = _confidence_level(final_avg_prob, final_avg_very_high, final_avg_high, final_avg_medium)
            fp.write(
                f"{genome}\t{contig}\t{seqname}\t"
                f"{avg_like_prob:.4f}\t{avg_like_conf}\t{viral_prob:.4f}\t{viral_conf}\t"
                f"{final_avg_prob:.4f}\t{final_avg_conf}\n"
            )

    logger.info(f"PST model inference completed.")

if __name__ == "__main__":
    main()
