import argparse
from CheckAMG.scripts.aggregate import aggregate
from CheckAMG.scripts.common.set_up_snakemake import create_output_dir, run_snakemake
from CheckAMG.scripts.common.args_formatter_logging import CustomHelpFormatter
from CheckAMG.scripts.common.args_formatter_logging import build_rerunnable_command

def add_aggregate_subcommand(parser, subparsers, scripts_dir, default_threads, pct_total_cpu, available_memory_gb, checkamg_version):
    aggregate_parser = subparsers.add_parser(
        "aggregate",
        help="Aggregate results from the annotate and de-novo modules into one report.",
        description="Aggregate results from the annotate and de-novo modules into one report.",
        formatter_class=CustomHelpFormatter,
    )

    inputs = aggregate_parser.add_argument_group("Input directories")
    inputs.add_argument(
        "-a",
        "--annotate-dir",
        type=str,
        required=True,
        help="Path to the CheckAMG annotate output directory.",
    )
    inputs.add_argument(
        "-d",
        "--denovo-dir",
        type=str,
        required=True,
        help="Path to the CheckAMG de-novo output directory.",
    )

    outputs = aggregate_parser.add_argument_group("Outputs")
    outputs.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Output directory to save the aggregated results.",
    )
    outputs.add_argument(
        "-pq",
        "--save-to-parquet",
        dest="save_to_parquet",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Write the final aggregated results as parquet files instead of TSV.\n"
            "Tables will be smaller files but not human readable without external tools.\n"
            "Recommended for large datasets."
        ),
    )

    resources = aggregate_parser.add_argument_group("Resources")
    resources.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=default_threads,
        help=f"Maximum number of threads allowed. Default is {pct_total_cpu}%% of available.",
    )
    resources.add_argument(
        "-m",
        "--mem",
        dest="mem",
        type=int,
        default=round(available_memory_gb * 0.80),
        help="Memory limit in GB. Default is 80%% of available.",
    )

    other = aggregate_parser.add_argument_group("Other")
    other.add_argument(
        "--debug",
        dest="debug",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Enable debug-level logging.",
    )

    def _run_aggregate(args, parser, this_parser=aggregate_parser, _scripts_dir=scripts_dir, _checkamg_version=checkamg_version):
        create_output_dir(args.output)

        config_path, logger, log_file_path = aggregate.generate_config(args)

        command_full_print = build_rerunnable_command(this_parser, args)

        run_snakemake(
            module_dir_name="aggregate",
            module_arg_name="aggregate",
            config_path=config_path,
            scripts_dir=_scripts_dir,
            args=args,
            checkamg_version=_checkamg_version,
            logger=logger,
            log_file_path=log_file_path,
            command=command_full_print,
        )

    aggregate_parser.set_defaults(func=_run_aggregate)
    return aggregate_parser