import argparse
import logging
from pathlib import Path
from typing import Optional
import psutil

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

from CheckAMG.scripts.denovo import denovo
from CheckAMG.scripts.denovo.config import (
    _KNN_FIELD,
    _MARGIN_FIELD,
    _MAX_AP_PAIRS_FIELD,
    _NEGATIVE_MINING_FIELD,
    _POSITIVE_MINING_FIELD,
    _CONTEXT_SIZE_FIELD,
)
from CheckAMG.scripts.common.set_up_snakemake import create_output_dir, run_snakemake
from CheckAMG.scripts.common.args_formatter_logging import (
    CustomHelpFormatter,
    build_rerunnable_command,
)

AVAIL_MEM_GB = psutil.virtual_memory().available / (1024 ** 3)
# DEFAULT_BATCH_SIZE = int((AVAIL_MEM_GB * 0.80))
# DEFAULT_NN_BATCH_SIZE = DEFAULT_BATCH_SIZE * 100
DEFAULT_BATCH_SIZE = 16
DEFAULT_NN_BATCH_SIZE = 4096

for name in ("numexpr", "NumExpr"):
    logging.getLogger(name).setLevel(logging.ERROR)


def _convert_paths(value: FilePath) -> Optional[Path]:
    if value is None:
        return None
    if isinstance(value, str):
        return Path(value)
    return value


_optional_file_path_field = attrs.field(
    default=None,
    validator=validators.optional(file_exists),
    converter=_convert_paths,
)

OptPath = Optional[Path]


@attrs.define
class EmbedArgs(AttrsDataclassUtilitiesMixin):
    """
    Object-based args for embedding mode.
    Fully importable and usable elsewhere.
    """

    model_ckpt: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    query_file: Path = attrs.field(validator=file_exists, converter=_convert_paths)
    output_file: Path = attrs.field(converter=_convert_paths)

    batch_size: int = attrs.field(default=DEFAULT_BATCH_SIZE, validator=positive_int)
    accelerator: str = attrs.field(
        default="auto", validator=validators.in_({"cpu", "gpu", "auto"})
    )
    devices: int = attrs.field(default=1, validator=positive_int)

    verbose: bool = True

    @classmethod
    def from_namespace(cls, ns: argparse.Namespace):
        return cls.from_dict(vars(ns))


@attrs.define
class InferenceArgs(AttrsDataclassUtilitiesMixin):
    """
    Object-based args for de-novo inference.
    Mirrors the CLI exactly but exposes a validated class.
    """

    model_ckpt: Path = attrs.field(validator=file_exists, converter=_convert_paths)

    output: Path = attrs.field(converter=_convert_paths)

    query_file: OptPath = _optional_file_path_field
    query_file_esm2: OptPath = _optional_file_path_field
    query_contigs: OptPath = attrs.field(default=None, converter=_convert_paths)
    query_bins: OptPath = attrs.field(default=None, converter=_convert_paths)
    query_proteins: OptPath = attrs.field(default=None, converter=_convert_paths)
    query_bin_proteins: OptPath = attrs.field(default=None, converter=_convert_paths)

    train_data_file: OptPath = _optional_file_path_field
    train_embed_file: OptPath = _optional_file_path_field
    train_index_file: OptPath = _optional_file_path_field
    train_labels_file: OptPath = _optional_file_path_field

    esm2_ckpt_dir: OptPath = _optional_file_path_field

    batch_size: int = attrs.field(default=DEFAULT_BATCH_SIZE, validator=positive_int)
    nn_batch_size: int = attrs.field(default=DEFAULT_NN_BATCH_SIZE, validator=positive_int)
    knn: int = _KNN_FIELD

    accelerator: str = attrs.field(
        default="auto", validator=validators.in_({"cpu", "gpu", "auto"})
    )
    devices: int = attrs.field(default=1, validator=positive_int)

    num_index_cells: Optional[int] = attrs.field(
        default=None,
        validator=validators.optional(positive_int),
    )

    num_probe_cells: Optional[int] = attrs.field(
        default=None,
        validator=validators.optional(positive_int),
    )

    seed: int = attrs.field(default=10241937, validator=positive_int)

    debug: bool = False

    @classmethod
    def from_namespace(cls, ns: argparse.Namespace):
        return cls.from_dict(vars(ns))


def add_denovo_subcommand(
    parser,
    subparsers,
    scripts_dir,
    default_threads,
    pct_total_cpu,
    available_memory_gb,
    checkamg_version,
):
    denovo_parser = subparsers.add_parser(
        "de-novo",
        help="Predict auxiliary viral genes with a protein-based genome language model, Protein Set Transformer (PST).",
        description="Predict auxiliary viral genes with a protein-based genome language model, Protein Set Transformer (PST).",
        formatter_class=CustomHelpFormatter,
    )
    
    fasta = denovo_parser.add_argument_group("FASTA inputs (slower, required if neither -i nor -e are provided)")
    fasta.add_argument(
        "-c", "--query-contigs", type=str,
        help=(
            "Input nucleotide contigs FASTA (.fna/.fasta; gzipped allowed)."
            ),
        )
    fasta.add_argument(
        "-C", "--query-bins", type=str,
        help=(
            "Folder of binned contig FASTAs (e.g. vMAGs with multiple contigs). "
            "Expects one .fna/.fasta (gzipped allowed) per bin."
            ),
        )
    fasta.add_argument(
        "-p", "--query-proteins", type=str,
        help=(
            "Input amino-acid FASTA from translated contigs (.faa/.fasta; gzipped allowed). "
            "Expected Prodigal headers: >[CONTIG]_[CDS] # START # END # FRAME # ..."
            ),
        )
    fasta.add_argument(
        "-P", "--query-bin-proteins", type=str,
        help=(
            "Folder of amino-acid FASTAs from translated binned contigs (.faa/.fasta; gzipped allowed). "
            "Expects one file per bin, each containing proteins from multiple contigs."
            ),
        )

    query = denovo_parser.add_argument_group("Embedding inputs (optional, faster)")
    query.add_argument(
        "-i", "--query-file", type=str,
        help=(
            "Optional path PST embedded query proteins underneath a field called 'ctx_ptn'. If not "
            "provided, either amino-acid or nucleotide FASTAs are required to compute ESM2->PST "
            "embeddings, or the precomputed graph-formatted ESM2 embeddings must be provided."
            ),
        )
    query.add_argument(
        "-e", "--query-file-esm2", type=str,
        help=(
            "Optional path to precomputed graph-formatted ESM2 embeddings for the query proteins. "
            "Will compute PST embeddings from these before nearest neighbor label propagation."
            ),
        )
    query.add_argument(
        "-f", "--fasta-file", type=str,
        help=(
            "Optional path to the FASTA file used to generate the embeddings found in -i/--query-file "
            "or -e/--query-file-esm2. If provided, the sequence names will be written alongside the "
            "predictions. If not, the predictions will be written in the order of embeddings "
            "with a generic sequence identifier."
        ),
    )

    inference = denovo_parser.add_argument_group("Inference settings")
    inference.add_argument(
        "-b", "--batch-size", dest="batch_size", type=int, default=DEFAULT_BATCH_SIZE,
        help="Batch size in number of scaffolds if PST embeddings need to be computed (default: %(default)s)",
    )
    inference.add_argument(
        "--nn-batch-size", dest="nn_batch_size", type=int, default=DEFAULT_NN_BATCH_SIZE,
        help="Batch size in number of proteins when doing kNN searches (default: %(default)s)",
    )
    inference.add_argument(
        "--num-index-cells", dest="num_index_cells", type=int,
        help=(
            "If a search index needs to be built from the training data, then this is the "
            "number of cells to partition the training data into for inverted lookup. The "
            "default is to determine this number automatically, but for very large dataset (>"
            "10M proteins) may require large amounts of time to create the index. Smaller values "
            "will lead to faster index building time at the expense of slower search time."
        ),
    )
    inference.add_argument(
        "--num-probe-cells", dest="num_probe_cells", type=int,
        help=(
            "The number of index cells to visit when searching for nearest neighbors. The "
            "default is to visit 1 cell. Increasing this can slow down search time but "
            "lead to more accurate results with large queries or indexes."
        ),
    )
    inference.add_argument(
        "-k", "--knn", type=int, default=20,
        help=(
            "Number of nearest neighbors for distance weighted voting class label prediction (default: %(default)s)"
        ),
    )

    reference = denovo_parser.add_argument_group(
        "Reference training data"
        )
    reference.add_argument(
        "-d", "--db-dir", type=str,
        help=(
            "Path to the CheckAMG de-novo database directory (downloaded by 'checkamg download', "
            "separate from the annotate database) containing pre-trained model inputs, or to the parent "
            "directory created by 'checkamg download' that contains it. If not provided, "
            "or using a custom model trained using 'checkamg train', then --train-data-file OR "
            "--train-embed-file OR (--train-index-file AND --train-labels-file), and --model-ckpt must be "
            "provided as reference."
        )
        )
    reference.add_argument(
        "-td", "--train-data-file", type=str,
        help=(
            "Path to graph-formatted .h5 file for training data (used to train model). If "
            "provided, then these data points will be embed using the trained PST model, then "
            "used to build a distance-querying index for nearest neighbor search."
        ),
    )
    reference.add_argument(
        "-te", "--train-embed-file", type=str,
        help=(
            "Path to train protein PST embeddings WITH the training labels stored as a H5 file. "
            "The PST embeddings must be under 'ctx_ptn' and the training protein labels must be "
            "under 'label'. If provided, this will be used to build a distance-querying index "
            "for nearest neighbor search."
        ),
        )
    reference.add_argument(
        "-tx", "--train-index-file", type=str,
        help=(
            "Path to precomputed train protein index constructed from the training protein PST "
            "embeddings. This is generated by faiss (`faiss.write_index`)."
        ),
        )
    reference.add_argument(
        "-tl", "--train-labels-file", type=str,
        help=(
            "Path to train protein labels stored in a H5 file under 'label'."
        ),
        )
    reference.add_argument(
        "-mc", "--model-ckpt", type=str,
        help=(
            "Path to the trained CheckAMG-PST model checkpoint (.ckpt file)."
        )
    )
    reference.add_argument(
        "-ec", "--esm2-ckpt-dir", dest="esm2_ckpt_dir", type=str,
        help=(
            "Optional path to a directory containing pre-downloaded ESM2 checkpoint files "
            "('esm2_t30_150M_UR50D.pt' and 'esm2_t30_150M_UR50D-contact-regression.pt'). "
            "If provided, these will be used instead of downloading them via torch hub. "
            "If not provided, the files will be downloaded to ~/.cache/torch/hub/checkpoints "
            "on first run (unless already cached there)."
        )
    )

    outputs = denovo_parser.add_argument_group("Outputs")
    outputs.add_argument(
        "-o", "--output", required=True, type=str,
        help=(
            "Path to output directory where results, intermediate files, and logs will be written."
        ),
    )

    resources = denovo_parser.add_argument_group("Resources")
    resources.add_argument(
        "-a", "--accelerator", type=str, default="auto", choices=["cpu", "gpu", "auto"],
        help="Accelerator to use (default: %(default)s, options: 'cpu', 'gpu')",
        )
    resources.add_argument(
        "-x", "--devices", type=int, default=1,
        help=(
            "Number of devices to use (default: %(default)s). For CPU setups, this is the number of cores/threads. "
            "For GPUs, multi-GPU is not supported (only the first GPU is used)."
        ),
    )
    resources.add_argument(
        "-m", "--mem", dest="mem", type=int, default=round(available_memory_gb * 0.80),
        help="Memory limit in GB. Default is 80%% of available. CheckAMG de-novo will crash if memory usage exceeds this limit.",
    )

    other = denovo_parser.add_argument_group("Other")
    other.add_argument(
        "--seed", dest="seed", type=int, default=111,
        help="Random seed"
    )
    other.add_argument(
        "--debug", dest="debug", action=argparse.BooleanOptionalAction, default=False,
        help="Enable debug-level logging"
    )


    def _run_denovo(args, parser, this_parser=denovo_parser):
        create_output_dir(args.output)

        config_path, logger, log_file_path, gpu_clamp_warning = denovo.generate_config(args)

        command_full_print = build_rerunnable_command(this_parser, args)

        run_snakemake(
            module_dir_name="denovo",
            module_arg_name="de-novo",
            config_path=config_path,
            scripts_dir=scripts_dir,
            args=args,
            checkamg_version=checkamg_version,
            logger=logger,
            log_file_path=log_file_path,
            command=command_full_print,
            post_command_warning=gpu_clamp_warning,
        )

    denovo_parser.set_defaults(func=_run_denovo)

    return denovo_parser