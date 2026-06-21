#!/usr/bin/env python3

import argparse
from argparse import Namespace

from CheckAMG.scripts.end_to_end import end_to_end
from CheckAMG.scripts.annotate.cli import (
    _parse_filter_presets,
    _validate_and_resolve_filter_presets,
)
from CheckAMG.scripts.common.args_formatter_logging import (
    CustomHelpFormatter,
    build_rerunnable_command,
)

DEFAULT_BATCH_SIZE = 16
DEFAULT_NN_BATCH_SIZE = 4096


def add_end_to_end_subcommand(parser, subparsers, scripts_dir, default_threads, pct_total_cpu, available_memory_gb, checkamg_version):
    e2e_parser = subparsers.add_parser(
        "end-to-end",
        help="Run annotate, de-novo, and aggregate in tandem into one organized output.",
        description=(
            "Run CheckAMG annotate, then de-novo using the proteins from annotate, and then aggregate, writing each module to its own sub-directory (01_annotate, 02_denovo, 03_aggregate) within one output."
        ),
        formatter_class=CustomHelpFormatter,
    )

    required = e2e_parser.add_argument_group("required arguments")
    required.add_argument(
        "-d", "--db-dir", dest="db_dir", type=str, required=True,
        help="Path to a directory containg both the CheckAMG annotate and CheckAMG de-novo databases.",
    )
    required.add_argument(
        "-o", "--output", dest="output", type=str, required=True,
        help="Output directory for all modules and combined results.",
    )

    inputs = e2e_parser.add_argument_group("input arguments")
    inputs.add_argument(
        "-i", "--input-contigs", dest="input_contigs", type=str, default=None,
        help="Input nucleotide contigs FASTA (.fna/.fasta; gzipped allowed).",
    )
    inputs.add_argument(
        "-I", "--input-bins", dest="input_bins", type=str, default=None,
        help=(
            "Directory containing binned contig FASTAs (e.g. vMAGs with multiple contigs).\n"
            "Expects one .fna/.fasta (gzipped allowed) per bin."
        ),
    )
    inputs.add_argument(
        "-p", "--input-proteins", dest="input_proteins", type=str, default=None,
        help=(
            "Input amino-acid FASTA from translated contigs (.faa/.fasta; gzipped allowed).\n"
            "Expected Prodigal headers: >[CONTIG]_[CDS] # START # END # FRAME # ..."
        ),
    )
    inputs.add_argument(
        "-P", "--input-bin-proteins", dest="input_bin_proteins", type=str, default=None,
        help=(
            "Directory containing amino-acid FASTAs from translated binned contigs (.faa/.fasta; gzipped allowed).\n"
            "Expects one file per bin, each containing proteins from multiple contigs."
        ),
    )
    inputs.add_argument(
        "--input-type", dest="input_type", type=str, choices=["nucl", "prot"], default="nucl",
        help=(
            "Input type: 'nucl' for nucleotide sequences or 'prot' for translated amino-acid sequences.\n"
            "When 'nucl', de-novo reuses the proteins predicted by the annotate module."
        ),
    )

    annotate_opts = e2e_parser.add_argument_group("annotate options")
    annotate_opts.add_argument(
        "-l", "--min-len", dest="min_len", type=int, default=5000,
        help="Minimum length (bp) of input contigs for them to be considered for analysis.",
    )
    annotate_opts.add_argument(
        "-f", "--min-orf", dest="min_orf", type=int, default=4,
        help="Minimum number of ORFs/proteins per contig for it to be considered for analysis.",
    )
    annotate_opts.add_argument(
        "-amg", "--min-amg-weight", dest="min_amg_weight", type=float, default=0.6,
        help="Minimum AMG weight (0.0-1.0) for a putative AMG annotation to be included in the final set.",
    )
    annotate_opts.add_argument(
        "-kf", "--keep-full-hmm-results", dest="keep_full_hmm_results",
        action=argparse.BooleanOptionalAction, default=False,
        help="Write all HMM search results for every hit in each database (large; use with --save-to-parquet).",
    )
    annotate_opts.add_argument(
        "--filter-ambig-regions", dest="filter_ambig_regions",
        action=argparse.BooleanOptionalAction, default=False,
        help="Exclude predictions that fall outside strict viral regions (extra conservative).",
    )
    annotate_opts.add_argument(
        "--filter-avg-arrays", dest="filter_avg_arrays",
        action=argparse.BooleanOptionalAction, default=True,
        help="Exclude AVG predictions that occur in contiguous runs (arrays).",
    )
    annotate_opts.add_argument(
        "--avg-array-limit", dest="avg_array_limit", type=int, default=3,
        help="If --filter-avg-arrays is enabled, exclude runs of AVGs of this length or more.",
    )
    annotate_opts.add_argument(
        "--filter-presets", dest="filter_presets", type=str,
        default="allow_glucan,allow_nucleotide,allow_methyl,allow_lipid",
        help=(
            "Comma-separated preset(s) to control annotation filtering. Valid presets: allow_glucan,\n"
            "allow_nucleotide, allow_methyl, allow_lipid, all_filters, no_filter."
        ),
    )

    denovo_opts = e2e_parser.add_argument_group("de-novo options")
    denovo_opts.add_argument(
        "-b", "--batch-size", dest="batch_size", type=int, default=DEFAULT_BATCH_SIZE,
        help="Batch size in number of scaffolds when computing PST embeddings (default: %(default)s).",
    )
    denovo_opts.add_argument(
        "--nn-batch-size", dest="nn_batch_size", type=int, default=DEFAULT_NN_BATCH_SIZE,
        help="Batch size in number of proteins when doing kNN searches (default: %(default)s).",
    )
    denovo_opts.add_argument(
        "--num-index-cells", dest="num_index_cells", type=int, default=None,
        help="Number of cells to partition the training data into when building a search index.",
    )
    denovo_opts.add_argument(
        "--num-probe-cells", dest="num_probe_cells", type=int, default=None,
        help="Number of index cells to visit when searching for nearest neighbors.",
    )
    denovo_opts.add_argument(
        "-k", "--knn", dest="knn", type=int, default=20,
        help="Number of nearest neighbors for distance-weighted voting (default: %(default)s).",
    )
    denovo_opts.add_argument(
        "--denovo-accelerator", dest="denovo_accelerator", type=str, default="auto",
        choices=["cpu", "gpu", "auto"],
        help=(
            "Accelerator for de-novo only (default: %(default)s). Annotate always runs on CPU.\n"
            "Use 'gpu' to run de-novo on a GPU while annotate uses CPU threads."
        ),
    )
    denovo_opts.add_argument(
        "--denovo-devices", dest="denovo_devices", type=int, default=None,
        help=(
            "Number of devices for de-novo (default: same as --threads when on CPU, capped at 1 on GPU).\n"
            "For CPU this is the number of cores/threads; multi-GPU is not supported."
        ),
    )

    outputs = e2e_parser.add_argument_group("output options")
    outputs.add_argument(
        "-pq", "--save-to-parquet", dest="save_to_parquet",
        action=argparse.BooleanOptionalAction, default=False,
        help="Write intermediate and final tables as parquet instead of TSV.",
    )

    resources = e2e_parser.add_argument_group("resources")
    resources.add_argument(
        "-t", "--threads", dest="threads", type=int, default=default_threads,
        help=f"Maximum number of threads allowed. Default is {pct_total_cpu}%% of available.",
    )
    resources.add_argument(
        "-m", "--mem", dest="mem", type=int, default=round(available_memory_gb * 0.80),
        help="Memory limit in GB. Default is 80%% of available.",
    )
    resources.add_argument(
        "--debug", dest="debug", action=argparse.BooleanOptionalAction, default=False,
        help="Enable debug-level logging.",
    )

    def _run_end_to_end(args, parser, this_parser=e2e_parser, _scripts_dir=scripts_dir, _checkamg_version=checkamg_version):
        if args.input_type == "nucl" and not args.input_contigs and not args.input_bins:
            parser.error("At least one of --input-contigs or --input-bins is required when --input-type is 'nucl'.")
        if args.input_type == "nucl" and (args.input_proteins or args.input_bin_proteins):
            parser.error("Cannot provide --input-proteins or --input-bin-proteins when --input-type is 'nucl'.")
        if args.input_type == "prot" and not args.input_proteins and not args.input_bin_proteins:
            parser.error("At least one of --input-proteins or --input-bin-proteins is required when --input-type is 'prot'.")
        if args.input_type == "prot" and (args.input_contigs or args.input_bins):
            parser.error("Cannot provide --input-contigs or --input-bins when --input-type is 'prot'.")
        if (args.input_contigs and args.input_proteins) or (args.input_bins and args.input_bin_proteins):
            parser.error("Cannot provide both --input-contigs/--input-bins and --input-proteins/--input-bin-proteins.")

        raw_list = _parse_filter_presets(args.filter_presets)
        if not raw_list:
            raw_list = ["allow_glucan", "allow_nucleotide", "allow_methyl", "allow_lipid"]
        elif raw_list == ["all_filters"]:
            raw_list = []
        effective = _validate_and_resolve_filter_presets(raw_list, parser)
        filter_presets = ",".join(effective)

        command = build_rerunnable_command(this_parser, args)

        common = Namespace(output=args.output, debug=args.debug)

        annotate_ns = Namespace(
            output=None,
            db_dir=args.db_dir,
            input_type=args.input_type,
            input_contigs=args.input_contigs,
            input_bins=args.input_bins,
            input_proteins=args.input_proteins,
            input_bin_proteins=args.input_bin_proteins,
            min_len=args.min_len,
            min_orf=args.min_orf,
            min_annot=0.20,
            cov_fraction=0.4,
            evalue=1e-5,
            bit_score=60,
            bitscore_fraction_heuristic=0.5,
            window_size=5000,
            min_flank_Vscore=10.0,
            min_window_avg_VL_score=3.0,
            use_hallmark=False,
            min_amg_weight=args.min_amg_weight,
            filter_ambig_regions=args.filter_ambig_regions,
            filter_avg_arrays=args.filter_avg_arrays,
            avg_array_limit=args.avg_array_limit,
            filter_presets=filter_presets,
            keep_full_hmm_results=args.keep_full_hmm_results,
            save_to_parquet=args.save_to_parquet,
            threads=args.threads,
            mem=args.mem,
            debug=args.debug,
        )

        denovo_devices = args.denovo_devices if args.denovo_devices is not None else args.threads
        denovo_ns = Namespace(
            output=None,
            db_dir=args.db_dir,
            query_contigs=None,
            query_bins=None,
            query_proteins=None,
            query_bin_proteins=None,
            query_file=None,
            query_file_esm2=None,
            fasta_file=None,
            train_data_file=None,
            train_embed_file=None,
            train_index_file=None,
            train_labels_file=None,
            model_ckpt=None,
            esm2_ckpt_dir=None,
            knn=args.knn,
            batch_size=args.batch_size,
            nn_batch_size=args.nn_batch_size,
            num_index_cells=args.num_index_cells,
            num_probe_cells=args.num_probe_cells,
            seed=111,
            accelerator=args.denovo_accelerator,
            devices=denovo_devices,
            mem=args.mem,
            debug=args.debug,
        )

        aggregate_ns = Namespace(
            output=None,
            annotate_dir=None,
            denovo_dir=None,
            save_to_parquet=args.save_to_parquet,
            threads=args.threads,
            mem=args.mem,
            debug=args.debug,
        )

        end_to_end.run_end_to_end(
            common=common,
            annotate_ns=annotate_ns,
            denovo_ns=denovo_ns,
            aggregate_ns=aggregate_ns,
            scripts_dir=_scripts_dir,
            checkamg_version=_checkamg_version,
            command=command,
        )

    e2e_parser.set_defaults(func=_run_end_to_end)
    return e2e_parser
