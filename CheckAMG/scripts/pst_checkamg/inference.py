import math
from pathlib import Path
from typing import NamedTuple

import faiss
import pst
import tables as tb
import torch
from pst.embed import main as esm2_embed
from pst.typing import OptTensor
from pst.utils import graphify
from pyfastatools import Parser
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm

from CheckAMG.scripts.pst_checkamg.cli import InferenceArgs
from CheckAMG.scripts.pst_checkamg.constants import (
    CHUNK_SIZE,
    EDGE_STRATEGY,
    EMBED_FIELD,
    FRAGMENT_SIZE,
    LABEL_FIELD,
    NUM_CLASSES,
    THRESHOLD,
)
from CheckAMG.scripts.pst_checkamg.data import GenomeDataModule
from CheckAMG.scripts.pst_checkamg.embed import embed
from CheckAMG.scripts.pst_checkamg.model import CheckAMGPST
from CheckAMG.scripts.pst_checkamg.utils import load_model


class Decomposed(NamedTuple):
    avg_like: torch.Tensor
    viral: torch.Tensor

    def __iter__(self):
        return zip(self.avg_like, self.viral)


def build_index(
    data: torch.Tensor, n_cells: int | None = None
) -> faiss.IndexIVFFlat:
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

    index.train(data)  # type: ignore
    index.add(data)  # type: ignore
    return index


def move_index_to_gpu(index: faiss.Index):
    res = faiss.StandardGpuResources()
    index = faiss.index_cpu_to_gpu(res, 0, index)
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


def decompose_labels(y: torch.Tensor) -> Decomposed:
    decomposer = torch.tensor(
        [
            [False, False],
            [False, True],
            [True, False],
            [True, True],
        ],
        device=y.device,
    )

    y_decomp = decomposer[y]

    return Decomposed(avg_like=y_decomp[:, 0], viral=y_decomp[:, 1])


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
) -> tuple[faiss.IndexIVFFlat, torch.Tensor]:
    if train_index_file is not None:
        if train_labels_file is None:
            raise ValueError(
                "BOTH the training index file and the training labels file must be provided is one is."
            )
        print("Reading precomputed train data index")
        index = faiss.read_index(train_index_file)
        y = _read_labels(train_labels_file)
    elif train_embed_file is not None:
        train_embed = _read_embeddings(train_embed_file)

        print("Building train data search index")
        index = build_index(train_embed, n_cells=n_cells)
        # TODO: save?

        # TODO: should maybe be in the decomposed format?
        y = _read_labels(train_embed_file)
    elif train_data_file is not None:
        datamodule = GenomeDataModule(train_file=train_data_file)

        datamodule.setup()

        train_embed_file = train_data_file.with_suffix(".PST-EMBED.h5")

        print("Embedding training data")
        embed(
            model=model,
            dataset=datamodule.train_dataset,
            output=train_embed_file,
            batch_size=batch_size,
        )

        train_embed = _read_embeddings(train_embed_file)
        index = build_index(train_embed, n_cells=n_cells)
        y: torch.Tensor = datamodule.train_dataset.dataset.y[  # type: ignore
            datamodule.train_dataset.indices  # type: ignore
        ]

    return index, y


def _esm2_embed(
    fasta_file: Path,
    devices: int,
) -> Path:
    outdir = Path("_checkAMGPST")
    esm2_embed.embed(
        input=fasta_file,
        outdir=outdir,
        model_cfg=esm2_embed.ModelArgs(esm=esm2_embed.ESM2Models.esm2_t30_150M),
        trainer_cfg=esm2_embed.TrainerArgs(
            devices=devices,
            precision="bf16-mixed" if torch.cuda.is_available() else 32,
        ),
    )

    output_file = outdir.joinpath("predictions.h5")
    output_file = output_file.rename(Path(fasta_file.name).with_suffix(".h5"))
    outdir.rmdir()

    graphfmt_file = output_file.with_suffix(".graphfmt.h5")

    # now graphify

    graphify.to_graph_format(
        io=graphify.IOArgs(
            file=output_file,
            fasta_file=fasta_file,
            output=graphfmt_file,
        ),
        optional=graphify.OptionalArgs(),
    )

    return graphfmt_file


def handle_query_inputs(
    model: CheckAMGPST,
    fasta_file: OptPath,
    query_file: OptPath,
    query_file_esm2: OptPath,
    batch_size: int,
    devices: int,
) -> torch.Tensor:
    if query_file is not None:
        return _read_embeddings(query_file)
    elif query_file_esm2 is not None:
        dataset = pst.LazyGenomeDataset(
            file=query_file_esm2,
            edge_strategy=EDGE_STRATEGY,
            chunk_size=CHUNK_SIZE,
            threshold=THRESHOLD,
            fragment_size=FRAGMENT_SIZE,
        )

        query_file = query_file_esm2.with_suffix(".PST-EMBED.h5")

        print("Embedding input query data")
        embed(
            model=model,
            dataset=dataset,
            output=query_file,
            batch_size=batch_size,
        )

        return _read_embeddings(query_file)
    elif fasta_file is not None:
        print("Embedding protein sequences with ESM2")
        return handle_query_inputs(
            model=model,
            fasta_file=None,
            query_file=None,
            query_file_esm2=_esm2_embed(fasta_file, devices),
            batch_size=batch_size,
            devices=devices,
        )

    raise ValueError("At least one file input required.")


def _handle_cpu_accelerator(args: InferenceArgs):
    if args.accelerator == "cpu" or (
        args.accelerator == "auto" and not torch.cuda.is_available()
    ):
        torch.set_num_threads(args.devices)
        faiss.omp_set_num_threads(args.devices)


def main():
    args = InferenceArgs.parse_from_cli()
    _handle_cpu_accelerator(args)

    # TODO: device
    model = load_model(ckpt=args.model_ckpt, device=args.accelerator)

    index, train_y = handle_training_inputs(
        model=model,
        train_data_file=args.train_data_file,
        train_embed_file=args.train_embed_file,
        train_index_file=args.train_index_file,
        train_labels_file=args.train_labels_file,
        batch_size=args.batch_size,
        n_cells=args.num_index_cells,
    )

    if args.num_probe_cells is not None:
        index.nprobe = args.num_probe_cells

    query_data = handle_query_inputs(
        model=model,
        fasta_file=args.fasta_file,
        query_file=args.query_file,
        query_file_esm2=args.query_file_esm2,
        batch_size=args.batch_size,
        devices=args.devices,
    )

    dist, nn_idx = knn_search(
        data=query_data, index=index, k=args.knn, batch_size=args.nn_batch_size
    )
    weights = compute_neighbor_weights(dist)

    nn_y = train_y[nn_idx]
    class_probs = decompose_class_probs(weighted_vote(nn_y, weights))

    if args.fasta_file is not None:
        seq_iter = (header.name for header in Parser(args.fasta_file).headers())
    else:
        seq_iter = (f"seq_{i}" for i in range(query_data.shape[0]))

    with open(args.output_file, "w") as fp:
        fp.write("sequence\tAVG-like prob\tViral prob\n")

        for seqname, (avg_like_prob, viral_prob) in zip(seq_iter, class_probs):
            avg_like_prob = float(avg_like_prob)
            viral_prob = float(viral_prob)
            fp.write(f"{seqname}\t{avg_like_prob:.4f}\t{viral_prob:.4f}\n")


if __name__ == "__main__":
    main()
