import argparse
from CheckAMG.scripts.train_pst import train
from CheckAMG.scripts.common.set_up_snakemake import create_output_dir, run_snakemake
from CheckAMG.scripts.common.args_formatter_logging import CustomHelpFormatter
from CheckAMG.scripts.common.args_formatter_logging import build_rerunnable_command

def add_train_subcommand(parser, subparsers, scripts_dir, default_threads, pct_total_cpu, available_memory_gb, checkamg_version):
    train_parser = subparsers.add_parser(
        "train",
        help="Finetune a PST model to predict whether proteins are viral and auxiliary in de-novo mode.",
        description="Finetune a PST model to predict whether proteins are viral and auxiliary in de-novo mode.",
        formatter_class=CustomHelpFormatter,
    )

    inputs = train_parser.add_argument_group("Input train/test data")
    inputs.add_argument(
        "-i",
        "--train-file",
        type=str,
        required=True,
        help="Path to graph-formatted .h5 file containing training data.",
    )
    inputs.add_argument(
        "-te",
        "--test-file",
        type=str,
        help="Path to graph-formatted .h5 file containing test data. If not provided, no test set will be used.",
    )
    inputs.add_argument(
        "-mc",
        "--model-ckpt",
        type=str,
        required=True,
        help="Path to the pretrained PST checkpoint (.ckpt file) to use for finetuning.",
    )

    outputs = train_parser.add_argument_group("Outputs")
    outputs.add_argument(
        "-o",
        "--output",
        type=str,
        help="Directory to save training outputs including logs, checkpoints, and trained model (default: %(default)s)",
        required=True,
    )
    outputs.add_argument(
        "-s",
        "--save-train-embed",
        dest="save_train_embed",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Store training protein embeddings in --output. If false, these will not be computed after training."
        ),
    )

    params = train_parser.add_argument_group("Training parameters")
    params.add_argument(
        "-b",
        "--batch-size",
        type=int,
        default=16,
        help="Batch size in number of scaffolds (default: %(default)s)",
    )
    params.add_argument(
        "--lr",
        type=float,
        default=1e-3,
        help="Learning rate for the optimizer (default: %(default)s, range: (0.0, 0.1])",
    )
    params.add_argument(
        "-e",
        "--max-epochs",
        type=int,
        default=25,
        help="Maximum number of epochs for training (default: %(default)s)",
    )
    params.add_argument(
        "--dropout",
        type=float,
        default=0.0,
        help="dropout (default: %(default)s, range: [0.0, 1.0])",
    )
    params.add_argument(
        "--margin",
        type=float,
        default=1.0,
        help="triplet loss margin (default: %(default)s)",
    )
    params.add_argument(
        "--knn",
        type=int,
        default=5,
        help="number of nearest neighbors for majority voting class label prediction (default: %(default)s)",
    )
    params.add_argument(
        "--max-ap-pairs",
        type=int,
        default=250_000,
        help="maximum number of anchor-positive pairs to consider at a time. All pairs with a semihard negative sample will contribute to the loss but will be computed in batches (default: %(default)s)",
    )
    params.add_argument(
        "--positive-mining-strategy",
        default="esm_easy",
        choices=["all", "hard", "esm_easy"],
        help=(
            "positive mining strategy (default: %(default)s) [choices: %(choices)s]. "
            "'esm_easy' = the nearest same-class protein in the input ESM embedding space; "
            "'hard' = the farthest same-class protein in PST embedding space; "
            "'all' = every same-class pair."
        ),
    )
    params.add_argument(
        "--negative-mining-strategy",
        default="semihard",
        choices=["hard", "semihard", "both"],
        help="negative mining strategy (default: %(default)s) [choices: %(choices)s]",
    )
    params.add_argument(
        "--opposite-class-negatives",
        dest="opposite_class_negatives",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Restrict negative sampling, within the hard/semihard candidate set, to the "
            "most-opposite viral/AVG class present for each anchor (priority: flip viral, then "
            "flip AVG, then flip both). (default: %(default)s)"
        ),
    )
    params.add_argument(
        "--class-weighting",
        default="class_freq",
        choices=["class_freq", "pair_freq"],
        help="How to scale the loss for each triplet based on the class labels. (default: %(default)s) [choices: 'class_freq' = inverse weighted freq of class labels, 'pair_freq' inverse weighted freq of pairs from the same class]",
    )
    params.add_argument(
        "--context-size",
        type=int,
        default=4,
        help=(
            "context size to create align short-context views of AVGs with a long context embedding "
            "output from the model. This value is in number of proteins and specifies the max number "
            "of proteins up and downstream of the AVG to use for contextualizing the AVG embedding. "
            "(default: %(default)s)"
        ),
    )

    resources = train_parser.add_argument_group("Resources")
    resources.add_argument(
        "-a",
        "--accelerator",
        type=str,
        default="auto",
        choices=["cpu", "gpu", "auto"],
        help="Accelerator to use for training (default: %(default)s, options: 'cpu', 'gpu')",
    )
    resources.add_argument(
        "-x",
        "--devices",
        type=int,
        default=1,
        help=(
            "Number of devices to use for training (default: %(default)s). For CPU setups, this is the number of "
            "cores/threads. For GPUs, multi-GPU is not supported (only the first GPU is used)."
        ),
    )
    resources.add_argument(
        "-m", "--mem", dest="mem", type=int, default=round(available_memory_gb * 0.80),
        help="Memory limit in GB. Default is 80%% of available. CheckAMG de-novo will crash if memory usage exceeds this limit.",
    )

    other = train_parser.add_argument_group("Other")
    other.add_argument(
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Add per-step logging to the progress bar instead of only per-epoch logging.",
    )
    other.add_argument(
        "--seed",
        metavar="INT",
        type=int,
        default=111,
        help="random seed for reproducibility (default: %(default)s)",
    )

    def _run_train(args, parser, this_parser=train_parser, _scripts_dir=scripts_dir, _checkamg_version=checkamg_version):
        create_output_dir(args.output)

        config_path, logger, log_file_path, gpu_clamp_warning = train.generate_config(args)

        command_full_print = build_rerunnable_command(this_parser, args)

        run_snakemake(
            module_dir_name="train_pst",
            module_arg_name="train",
            config_path=config_path,
            scripts_dir=_scripts_dir,
            args=args,
            checkamg_version=_checkamg_version,
            logger=logger,
            log_file_path=log_file_path,
            command=command_full_print,
            post_command_warning=gpu_clamp_warning,
        )

    train_parser.set_defaults(func=_run_train)
    return train_parser
