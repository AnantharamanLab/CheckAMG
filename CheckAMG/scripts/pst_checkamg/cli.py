import argparse
from pathlib import Path
from typing import Optional

import attrs
from attrs import validators
from pst.typing import FilePath
from pst.utils.attrs.dataclass_utils import AttrsDataclassUtilitiesMixin
from pst.utils.attrs.validators import (
    closed_unit_interval,
    file_exists,
    optionally_existing_file,
    positive_float,
    positive_int,
)

from CheckAMG.scripts.pst_checkamg.config import (
    _CONTEXT_SIZE_FIELD,
    _KNN_FIELD,
    _MARGIN_FIELD,
    _MAX_AP_PAIRS_FIELD,
    _N_ATTN_HEADS_FIELD,
    _N_GLOBAL_ATTN_LAYERS_HEAD,
    _NEGATIVE_MINING_FIELD,
    _POSITIVE_MINING_FIELD,
)


def _convert_paths(value: FilePath) -> Path:
    if isinstance(value, str):
        return Path(value)
    return value


_WEIGHTING_FIELD = attrs.field(default="class_freq")


@attrs.define
class Args(AttrsDataclassUtilitiesMixin):
    model_ckpt: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    train_file: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    test_file: Optional[Path] = attrs.field(
        default=None, validator=optionally_existing_file, converter=_convert_paths
    )
    train_embed_output: Optional[Path] = attrs.field(
        default=None, validator=optionally_existing_file, converter=_convert_paths
    )
    batch_size: int = attrs.field(default=16, validator=positive_int)
    lr: float = attrs.field(
        default=1e-3,
        validator=validators.and_(
            positive_float, validators.gt(0.0), validators.le(1.0)
        ),
    )
    logdir: Path = attrs.field(
        default=Path("lightning_logs"), converter=_convert_paths
    )
    max_epochs: int = attrs.field(default=25, validator=positive_int)

    accelerator: str = attrs.field(
        default="auto", validator=validators.in_({"cpu", "gpu", "auto"})
    )

    devices: int = attrs.field(default=1, validator=positive_int)

    dropout: float = attrs.field(default=0.0, validator=closed_unit_interval)

    knn: int = _KNN_FIELD
    margin: float = _MARGIN_FIELD
    n_attn_layers: int = _N_GLOBAL_ATTN_LAYERS_HEAD
    n_attn_heads: int = _N_ATTN_HEADS_FIELD
    max_ap_pairs: int | None = _MAX_AP_PAIRS_FIELD
    verbose: bool = attrs.field(default=False)
    seed: int = attrs.field(default=123, validator=positive_int)

    positive_mining_strategy: str = _POSITIVE_MINING_FIELD
    negative_mining_strategy: str = _NEGATIVE_MINING_FIELD

    class_weighting: str = _WEIGHTING_FIELD
    context_size: int = _CONTEXT_SIZE_FIELD

    @classmethod
    def parse_from_cli(cls):
        parser = argparse.ArgumentParser(
            description="Finetune a PST model to predict whether proteins are (A) viral and (B) auxiliary."
        )

        parser.add_argument(
            "-i",
            "--train-file",
            metavar="FILE",
            type=Path,
            required=True,
            help="Path to graph-formatted .h5 file containing training data.",
        )

        parser.add_argument(
            "-t",
            "--test-file",
            metavar="FILE",
            type=Path,
            help="Path to graph-formatted .h5 file containing test data. If not provided, no test set will be used.",
        )

        parser.add_argument(
            "--train-embed-output",
            metavar="FILE",
            type=Path,
            help="Path to store training protein embeddings. If not specified, these will not be computed after training.",
        )

        parser.add_argument(
            "-m",
            "--model-ckpt",
            metavar="FILE",
            type=Path,
            required=True,
            help="Path to the pretrained PST checkpoint (.ckpt file) to use for finetuning.",
        )

        parser.add_argument(
            "-b",
            "--batch-size",
            metavar="INT",
            type=int,
            default=16,
            help="Batch size in number of scaffolds (default: %(default)s)",
        )

        parser.add_argument(
            "--lr",
            metavar="FLOAT",
            type=float,
            default=1e-3,
            help="Learning rate for the optimizer (default: %(default)s, range: (0.0, 1.0])",
        )

        parser.add_argument(
            "-o",
            "--logdir",
            metavar="DIR",
            type=Path,
            default="lightning_logs",
            help="Directory to save training logs and checkpoints (default: %(default)s)",
        )

        parser.add_argument(
            "-e",
            "--max-epochs",
            metavar="INT",
            type=int,
            default=25,
            help="Maximum number of epochs for training (default: %(default)s)",
        )

        parser.add_argument(
            "-a",
            "--accelerator",
            metavar="STR",
            type=str,
            default="auto",
            choices=["cpu", "gpu", "auto"],
            help="Accelerator to use for training (default: %(default)s, options: 'cpu', 'gpu')",
        )

        parser.add_argument(
            "-d",
            "--devices",
            metavar="INT",
            type=int,
            default=1,
            help="Number of devices to use for training (default: %(default)s). For CPU setups, this is the number of cores/threads. For GPUs, this is the number of GPU devices.",
        )

        parser.add_argument(
            "--dropout",
            metavar="FLOAT",
            type=float,
            default=0.0,
            help="dropout (default: %(default)s, range: [0.0, 1.0])",
        )

        parser.add_argument(
            "--margin",
            metavar="FLOAT",
            type=float,
            default=1.0,
            help="triplet loss margin (default: %(default)s)",
        )

        parser.add_argument(
            "--knn",
            metavar="INT",
            type=int,
            default=5,
            help="number of nearest neighbors for majority voting class label prediction (default: %(default)s)",
        )

        parser.add_argument(
            "--n-attn-layers",
            metavar="INT",
            type=int,
            default=10,
            help="number global attention layers (default: %(default)s)",
        )

        parser.add_argument(
            "--n-attn-heads",
            metavar="INT",
            type=int,
            default=8,
            help="number global attention layers (default: %(default)s)",
        )

        parser.add_argument(
            "--max-ap-pairs",
            metavar="INT|NONE",
            type=int,
            default=250_000,
            help="maximum number of anchor-positive pairs to consider at a time. All pairs with a semihard negative sample will contribute to the loss but will be computed in batches (default: %(default)s)",
        )

        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Add per-step logging to the progress bar instead of only per-epoch logging.",
        )

        parser.add_argument(
            "--seed",
            metavar="INT",
            type=int,
            default=123,
            help="random seed for reproducibility (default: %(default)s)",
        )

        parser.add_argument(
            "--positive-mining-strategy",
            metavar="",
            default="all",
            choices=["all", "hard"],
            help="positive mining strategy (default: %(default)s) [choices: %(choices)s]",
        )

        parser.add_argument(
            "--negative-mining-strategy",
            metavar="",
            default="semihard",
            choices=["hard", "semihard", "both"],
            help="negative mining strategy (default: %(default)s) [choices: %(choices)s]",
        )

        parser.add_argument(
            "--class-weighting",
            metavar="",
            default="class_freq",
            choices=["class_freq", "pair_freq"],
            help="How to scale the loss for each triplet based on the class labels. (default: %(default)s) [choices: 'class_freq' = inverse weighted freq of class labels, 'pair_freq' inverse weighted freq of pairs from the same class]",
        )

        parser.add_argument(
            "--context-size",
            metavar="INT",
            type=int,
            default=4,
            help=(
                "context size to create align short-context views of AVGs with a long context embedding "
                "output from the model. This value is in number of proteins and specifies the max number "
                "of proteins up and downstream of the AVG to use for contextualizing the AVG embedding. "
                "(default: %(default)s)"
            ),
        )

        args = parser.parse_args()
        return cls.from_dict(vars(args))


_optional_file_path_field = attrs.field(
    default=None,
    validator=validators.optional(file_exists),
    converter=_convert_paths,
)

OptPath = Path | None


@attrs.define
class InferenceArgs(AttrsDataclassUtilitiesMixin):
    model_ckpt: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    output_file: Path = attrs.field(converter=_convert_paths)

    # INPUTS
    fasta_file: OptPath = _optional_file_path_field
    query_file: OptPath = _optional_file_path_field
    query_file_esm2: OptPath = _optional_file_path_field

    # TRAIN DATA
    # need any combo of the following
    # 1. full train dataset: compute PST embeddings -> build index, use y
    # 2. train PST embeddings + y: build index, then use y
    # 3. index + y

    train_data_file: OptPath = _optional_file_path_field
    train_embed_file: OptPath = _optional_file_path_field
    train_index_file: OptPath = _optional_file_path_field
    train_labels_file: OptPath = _optional_file_path_field

    # ETC
    batch_size: int = attrs.field(default=16, validator=positive_int)
    nn_batch_size: int = attrs.field(default=4096, validator=positive_int)

    accelerator: str = attrs.field(
        default="auto", validator=validators.in_({"cpu", "gpu", "auto"})
    )
    devices: int = attrs.field(default=1, validator=positive_int)

    num_index_cells: int | None = attrs.field(
        default=None, validator=validators.optional(positive_int)
    )

    num_probe_cells: int | None = attrs.field(
        default=None, validator=validators.optional(positive_int)
    )

    knn: int = _KNN_FIELD

    @classmethod
    def parse_from_cli(cls):
        parser = argparse.ArgumentParser(
            description="Predict viral and auxiliary protein labels based on distance to training reference data."
        )

        query_group = parser.add_argument_group("QUERY INPUTS")

        query_group.add_argument(
            "-i",
            "--query-file",
            metavar="FILE",
            type=Path,
            # required=True,
            help=(
                "Optional path PST embedded query proteins underneath a field called 'ctx_ptn'. If "
                "not provided, either the FASTA file must be provided to fully compute ESM2->PST "
                "embeddings, or the precomputed graph-formatted ESM2 embeddings must be provided, "
                "which get embedded into the PST embedding space before nearest neighbor label "
                "propagation."
            ),
        )

        query_group.add_argument(
            "-e",
            "--query-file-esm2",
            metavar="FILE",
            type=Path,
            # required=True,
            help=(
                "Optional path to precomputed graph-formatted ESM2 embeddings for the query proteins. "
                "Will compute PST embeddings from these before nearest neighbor label propagation."
            ),
        )

        query_group.add_argument(
            "-f",
            "--fasta-file",
            metavar="FILE",
            type=Path,
            help=(
                "Path to FASTA file used to generate the embeddings found in -i/--query-file or "
                "-e/--query-file-esm2. If provided, the sequence names will be written alongside the "
                "predictions. If not, the predictions will be written in the order of embeddings "
                "with a generic sequence identifier. REQUIRED IF: neither -i nor -e are provided"
            ),
        )

        train_group = parser.add_argument_group("TRAIN INPUTS")
        train_group.add_argument(
            "-td",
            "--train-data-file",
            metavar="FILE",
            type=Path,
            help=(
                "Path to graph-formatted .h5 file for training data (used to train model). If "
                "provided, then these data points will be embed using the trained PST model, then "
                "used to build a distance-querying index for nearest neighbor search. TRAINING INPUTS: "
                "Either --train-data-file OR --train-embed-file OR (--train-index-file AND "
                "--train-labels-file) must be provided as reference."
            ),
        )

        train_group.add_argument(
            "-t",
            "--train-embed-file",
            metavar="FILE",
            type=Path,
            help=(
                "Path to train protein PST embeddings WITH the training labels stored as a H5 file. "
                "The PST embeddings must be under 'ctx_ptn' and the training protein labels must be "
                "under 'label'. If provided, this will be used to build a distance-querying index "
                "for nearest neighbor search. TRAINING INPUTS: Either --train-data-file OR "
                "--train-embed-file OR (--train-index-file AND --train-labels-file) must be provided "
                "as reference."
            ),
        )

        train_group.add_argument(
            "-x",
            "--train-index-file",
            metavar="FILE",
            type=Path,
            help=(
                "Path to precomputed train protein index constructed from the training protein PST "
                "embeddings. This is generated by faiss (`faiss.write_index`). TRAINING INPUTS: "
                "Either --train-data-file OR --train-embed-file OR (--train-index-file AND "
                "--train-labels-file) must be provided as reference."
            ),
        )

        train_group.add_argument(
            "-l",
            "--train-labels-file",
            metavar="FILE",
            type=Path,
            help=(
                "Path to train protein labels stored in a H5 file under 'label'. TRAINING INPUTS: "
                "Either --train-data-file OR --train-embed-file OR (--train-index-file AND "
                "--train-labels-file) must be provided as reference."
            ),
        )

        other_io_group = parser.add_argument_group("OTHER IO")

        other_io_group.add_argument(
            "-o",
            "--output-file",
            metavar="FILE",
            type=Path,
            required=True,
            help="REQUIRED: Path to output TSV file that will contain predictions",
        )

        other_io_group.add_argument(
            "-m",
            "--model-ckpt",
            metavar="FILE",
            type=Path,
            required=True,
            help="REQUIRED: Path to the trained CheckAMGPST checkpoint (.ckpt file)",
        )

        misc_group = parser.add_argument_group("MISC")

        misc_group.add_argument(
            "-b",
            "--batch-size",
            metavar="INT",
            type=int,
            default=16,
            help="Batch size in number of scaffolds if PST embeddings need to be computed (default: %(default)s)",
        )

        misc_group.add_argument(
            "--nn-batch-size",
            metavar="INT",
            type=int,
            default=4096,
            help="Batch size in number of proteins when doing kNN searches (default: %(default)s)",
        )

        misc_group.add_argument(
            "-a",
            "--accelerator",
            metavar="STR",
            type=str,
            default="auto",
            choices=["cpu", "gpu", "auto"],
            help="Accelerator to use (default: %(default)s, options: 'cpu', 'gpu')",
        )

        misc_group.add_argument(
            "-d",
            "--devices",
            metavar="INT",
            type=int,
            default=1,
            help="Number of devices to use (default: %(default)s). For CPU setups, this is the number of cores/threads. For GPUs, this is the number of GPU devices.",
        )

        misc_group.add_argument(
            "-k",
            "--knn",
            metavar="INT",
            type=int,
            default=5,
            help="Number of nearest neighbors for distance weighted voting class label prediction (default: %(default)s)",
        )

        misc_group.add_argument(
            "--num-index-cells",
            metavar="INT",
            type=int,
            help=(
                "If a search index needs to be built from the training data, then this is the "
                "number of cells to partition the training data into for inverted lookup. The "
                "default is to determine this number automatically, but for very large dataset (>"
                "10M proteins) may require large amounts of time to create the index. Smaller values "
                "will lead to faster index building time at the expense of slower search time."
            ),
        )

        misc_group.add_argument(
            "--num-probe-cells",
            metavar="INT",
            type=int,
            help=(
                "The number of index cells to visit when searching for nearest neighbors. The "
                "default is to visit 1 cell. Increasing this can slow down search time but "
                "lead to more accurate results with large queries or indexes."
            ),
        )

        args = parser.parse_args()
        return cls.from_dict(vars(args))

    def __attrs_post_init__(self):
        if (
            self.query_file is None
            and self.query_file_esm2 is None
            and self.fasta_file
        ):
            raise ValueError(
                "At least one of query_file, query_file_esm2, and fasta_file must be provided."
            )

        if (
            self.train_data_file is None
            and self.train_embed_file is None
            and self.train_index_file is None
            and self.train_labels_file is None
        ):
            raise ValueError(
                "One of the following set of training data inputs must be provided: (1) The full "
                "graph-formatted training dataset that includes the ESM2 embeddings and the protein "
                "labels. (2) Precomputed PST protein embeddings for the training proteins AND the "
                "protein labels all stored in the same file. (3) A precomputed faiss index file AND "
                "a separate file containing the protein labels for the training proteins."
            )

        if (
            self.train_index_file is not None and self.train_labels_file is None
        ) or (self.train_index_file is None and self.train_labels_file is not None):
            raise ValueError(
                "BOTH the train_index_file and train_labels_file must be provided together."
            )


@attrs.define
class EmbedArgs(AttrsDataclassUtilitiesMixin):
    model_ckpt: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    query_file: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    output_file: Path = attrs.field(converter=_convert_paths)

    batch_size: int = attrs.field(default=16, validator=positive_int)

    accelerator: str = attrs.field(
        default="auto", validator=validators.in_({"cpu", "gpu", "auto"})
    )

    devices: int = attrs.field(default=1, validator=positive_int)
    verbose: bool = True

    @classmethod
    def parse_from_cli(cls):
        parser = argparse.ArgumentParser(
            description="Finetune a PST model to predict whether proteins are (A) viral and (B) auxiliary."
        )

        parser.add_argument(
            "-i",
            "--query-file",
            metavar="FILE",
            type=Path,
            required=True,
            help="Path to graph-formatted .h5 file containing query data.",
        )

        parser.add_argument(
            "-o",
            "--output-file",
            metavar="FILE",
            type=Path,
            help="Path to output TSV file that will contain predictions",
        )

        parser.add_argument(
            "-m",
            "--model-ckpt",
            metavar="FILE",
            type=Path,
            required=True,
            help="Path to the trained CheckAMGPST checkpoint (.ckpt file)",
        )

        parser.add_argument(
            "-b",
            "--batch-size",
            metavar="INT",
            type=int,
            default=16,
            help="Batch size in number of scaffolds (default: %(default)s)",
        )

        parser.add_argument(
            "-a",
            "--accelerator",
            metavar="STR",
            type=str,
            default="auto",
            choices=["cpu", "gpu", "auto"],
            help="Accelerator to use for training (default: %(default)s, options: 'cpu', 'gpu')",
        )

        parser.add_argument(
            "-d",
            "--devices",
            metavar="INT",
            type=int,
            default=1,
            help="Number of devices to use for training (default: %(default)s). For CPU setups, this is the number of cores/threads. For GPUs, this is the number of GPU devices.",
        )

        parser.add_argument(
            "--verbose",
            action="store_true",
            help="Show embedding progress bar.",
        )

        args = parser.parse_args()
        return cls.from_dict(vars(args))
