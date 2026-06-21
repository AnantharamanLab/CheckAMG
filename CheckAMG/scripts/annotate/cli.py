#!/usr/bin/env python3

import argparse
from CheckAMG.scripts.annotate import annotate
from CheckAMG.scripts.common.set_up_snakemake import create_output_dir, run_snakemake
from CheckAMG.scripts.common.args_formatter_logging import CustomHelpFormatter
from CheckAMG.scripts.common.args_formatter_logging import build_rerunnable_command


VALID_FILTER_PRESETS = {
    "all_filters",
    "allow_glucan",
    "allow_nucleotide",
    "allow_methyl",
    "allow_lipid",
    "no_filter",
}


def _parse_filter_presets(preset_str):
    if preset_str is None:
        return []
    parts = [p.strip() for p in preset_str.split(",") if p.strip()]
    return parts


def _validate_and_resolve_filter_presets(presets, parser):
    unknown = [p for p in presets if p not in VALID_FILTER_PRESETS]
    if unknown:
        parser.error(
            f"Unknown --filter-presets value(s): {', '.join(unknown)}. "
            f"Valid options are: {', '.join(sorted(VALID_FILTER_PRESETS))}."
        )

    if len(presets) > 1 and "all_filters" in presets:
        parser.error(
            "Cannot combine 'all_filters' with other filter presets. "
            "Use only 'all_filters' or specify non-default combinations."
        )

    if len(presets) > 1 and "all_filters" in presets:
        parser.error("Cannot combine 'all_filters' with other filter presets. Use only 'all_filters'.")

    return list(presets)

def add_annotate_subcommand(parser, subparsers, scripts_dir, default_threads, pct_total_cpu, available_memory_gb, checkamg_version):
    annotate_parser = subparsers.add_parser(
        "annotate",
        help="Predict and curate auxiliary genes using functional annotations and genomic context.",
        description="Predict and curate auxiliary genes in viral genomes based on functional annotations and genomic context.",
        formatter_class=CustomHelpFormatter,
    )

    required = annotate_parser.add_argument_group("required arguments")
    required.add_argument(
        "-d",
        "--db-dir",
        dest="db_dir",
        type=str,
        required=True,
        help=(
            "Path to the CheckAMG annotate database, or to the parent directory created by "
            "'checkamg download' that contains it."
        ),
    )
    required.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        required=True,
        help="Output directory for all generated files and folders.",
    )

    inputs = annotate_parser.add_argument_group("input arguments")
    inputs.add_argument(
        "-i",
        "--input-contigs",
        dest="input_contigs",
        type=str,
        default=None,
        help=(
            "Input nucleotide contigs FASTA (.fna/.fasta; gzipped allowed)."
        ),
    )
    inputs.add_argument(
        "-I",
        "--input-bins",
        dest="input_bins",
        type=str,
        default=None,
        help=(
            "Folder of binned contig FASTAs (e.g. vMAGs with multiple contigs).\n"
            "Expects one .fna/.fasta (gzipped allowed) per bin."
        ),
    )
    inputs.add_argument(
        "-p",
        "--input-proteins",
        dest="input_proteins",
        type=str,
        default=None,
        help=(
            "Input amino-acid FASTA from translated contigs (.faa/.fasta; gzipped allowed).\n"
            "Expected Prodigal headers: >[CONTIG]_[CDS] # START # END # FRAME # ..."
        ),
    )
    inputs.add_argument(
        "-P",
        "--input-bin-proteins",
        dest="input_bin_proteins",
        type=str,
        default=None,
        help=(
            "Folder of amino-acid FASTAs from translated binned contigs (.faa/.fasta; gzipped allowed).\n"
            "Expects one file per bin, each containing proteins from multiple contigs."
        ),
    )

    inputs.add_argument(
        "--input-type",
        dest="input_type",
        type=str,
        choices=["nucl", "prot"],
        default="nucl",
        help=(
            "Input type: 'nucl' for nucleotide sequences or 'prot' for translated amino-acid sequences.\n"
            "Providing proteins instead of nucleotide sequences skips pyrodigal-gv, and annotations/contextual\n"
            "analyses are performed using the provided proteins. So ensure all proteins from contigs/bins are\n"
            "included and that headers are formatted as expected (see --input-proteins)."
        ),
    )

    thresholds = annotate_parser.add_argument_group("thresholds and HMMsearch settings")
    thresholds.add_argument(
        "-l",
        "--min-len",
        dest="min_len",
        type=int,
        default=5000,
        help="Minimum length (bp) of input contigs for them to be considered for analysis.",
    )
    thresholds.add_argument(
        "-f",
        "--min-orf",
        dest="min_orf",
        type=int,
        default=4,
        help="Minimum number of ORFs/proteins per contig for it to be considered for analysis.",
    )
    thresholds.add_argument(
        "-a",
        "--min-annot",
        dest="min_annot",
        type=float,
        default=0.20,
        help=(
            "Minimum fraction (0.0-1.0) of genes per contig that must receive an annotation\n"
            "to be considered for contextual analysis."
        ),
    )
    thresholds.add_argument(
        "-c",
        "--cov-fraction",
        dest="cov_fraction",
        type=float,
        default=0.4,
        help=(
            "Minimum covered fraction (0.0-1.0) of HMM profiles required to report hits."
        ),
    )
    thresholds.add_argument(
        "-e",
        "--evalue",
        dest="evalue",
        type=float,
        default=1e-5,
        help="Maximum fallback E-value for HMM hits when database-provided cutoffs are unavailable.",
    )
    thresholds.add_argument(
        "-b",
        "--bitscore",
        dest="bit_score",
        type=int,
        default=60,
        help="Minimum fallback bit score for HMM hits when database-provided cutoffs are unavailable.",
    )
    thresholds.add_argument(
        "-bf",
        "--bitscore-fraction-heuristic",
        dest="bitscore_fraction_heuristic",
        type=float,
        default=0.5,
        help=(
            "Retain HMM hits scoring at least this fraction (0.0-1.0) of its database-provided threshold\n"
            "during heuristic filtering."
        ),
    )

    context = annotate_parser.add_argument_group("genomic context settings")
    context.add_argument(
        "-w",
        "--window-size",
        dest="window_size",
        type=int,
        default=5000,
        help="Window size (bp) for local average VL-score calculation.",
    )
    context.add_argument(
        "-v",
        "--min-flank-vscore",
        dest="min_flank_Vscore",
        type=float,
        default=10.0,
        help=(
            "Minimum V-score (0.0-10.0) required in flanking regions to verify viral origin\n"
            "and reduce host-contamination artifacts (higher = more viral-like)."
        ),
    )
    context.add_argument(
        "-vl",
        "--min-window-avg-vlscore",
        dest="min_window_avg_VL_score",
        type=float,
        default=3.0,
        help=(
            "Minimum average VL-score within the specified window size around a gene to be\n"
            "considered a viral region (higher = more viral-like)."
        ),
    )
    context.add_argument(
        "-ha",
        "--use-hallmark",
        dest="use_hallmark",
        default=False,
        action=argparse.BooleanOptionalAction,
        help=(
            "Use viral hallmark genes instead of V-scores when evaluating flanks.\n"
            "Enable to be extra conservative."
        ),
    )

    filters = annotate_parser.add_argument_group("filtering settings")
    filters.add_argument(
        "-amg",
        "--min-amg-weight",
        default=0.6,
        type=float,
        dest="min_amg_weight",
        help=(
            "Minimum AMG weight (0.0-1.0) for a putative AMG annotation to be included in the final\n"
            "set of predicted AMGs. AMG weight is a capped geometric mean of two scores: the ratio of\n"
            "auxiliary-like to non-auxiliary-like metabolic annotations for a given function, and its\n"
            "rarity in viral genomes (VL-score). Higher values (more conservative) indicate that the\n"
            "annotated function is predominantly associated with metabolic pathways and is uncommon in\n"
            "viral genomes, broadly. Lower values (more permissive) indicate functions that are either\n"
            "very common in viral genomes (likely a non-auxiliary role), or predominantly associated\n"
            "with non-metabolic or essential pathways, or both."
        ),
    )
    filters.add_argument(
        "--filter-ambig-regions",
        dest="filter_ambig_regions",
        default=False,
        action=argparse.BooleanOptionalAction,
        help=(
            "Exclude predictions that fall outside strict viral regions (inside ambiguous regions).\n"
            "Strict viral regions are identified from window-average VL-scores and then refined using\n"
            "per-gene V-scores (see --min-window-avg-vlscore and --min-flank-vscore) or viral hallmark\n"
            "genes if --use-hallmark is enabled (stricter, lower recall). When enabled, any prediction\n"
            "not overlapping a strict viral region is filtered out. Disabled by default because it can\n"
            "be too strict when annotation rate is low but other viral origin signals are strong.\n"
            "Enable to be extra conservative."
        ),
    )
    filters.add_argument(
        "--filter-avg-arrays",
        dest="filter_avg_arrays",
        default=True,
        action=argparse.BooleanOptionalAction,
        help="Exclude AVG predictions that occur in contiguous runs (arrays), which suggests\n"
        "non-auxiliary function.",
    )
    filters.add_argument(
        "--avg-array-limit",
        dest="avg_array_limit",
        type=int,
        default=3,
        help="If --filter-avg-arrays is enabled, exclude runs of AVGs of this length or more.",
    )
    filters.add_argument(
        "--filter-presets",
        dest="filter_presets",
        type=str,
        default="allow_glucan,allow_nucleotide,allow_methyl,allow_lipid",
        help=(
            "Comma-separated preset(s) to control annotation filtering, regardless of AMG weight.\n"
            "Valid presets:\n"
            "* allow_glucan (keep glycosyltransferase, glycoside-hydrolase, and related annotations)\n"
            "* allow_nucleotide (keep nucleotide metabolism annotations)\n"
            "* allow_methyl (keep methylase/methyltransferase annotations)\n"
            "* allow_lipid (keep lipopolysaccharide and phospholipid-related annotations)\n"
            "* all_filters (filter annotations in all above categories, strict)\n"
            "* no_filter (disable all filtering including essential viral anntoations, not recommended).\n"
            "Example: --filter-presets allow_glucan,allow_nucleotide."
            )
    )

    outputs = annotate_parser.add_argument_group("output files")
    outputs.add_argument(
        "-kf",
        "--keep-full-hmm-results",
        dest="keep_full_hmm_results",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Write all HMM search results for every hit in each database.\n"
            "By default, only the top hit per protein per database is written to reduce file size.\n"
            "Not recommended for large inputs unless --save-as-parquet is used."
        ),
    )

    outputs.add_argument(
        "-pq",
        "--save-to-parquet",
        dest="save_to_parquet",
        action=argparse.BooleanOptionalAction,
        default=False,
        help=(
            "Write intermediate and final tables as parquet files instead of TSV.\n"
            "Tables will be smaller files but not human readable without external tools.\n"
            "Recommended for large datasets."
        ),
    )

    resources = annotate_parser.add_argument_group("resources")
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
        help="Memory limit in GB. Default is 80%% of available. CheckAMG annotate will crash if memory usage exceeds this limit.",
    )
    resources.add_argument(
        "--debug",
        dest="debug",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Enable debug-level logging.",
    )

    def _run_annotate(args, parser, this_parser=annotate_parser, _scripts_dir=scripts_dir, _checkamg_version=checkamg_version):
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
            raw_list = ["allow_glucan","allow_nucleotide","allow_methyl","allow_lipid"]
        elif raw_list == ["all_filters"]:
            raw_list = []
        effective = _validate_and_resolve_filter_presets(raw_list, parser)
        args.filter_presets = ",".join(effective)

        create_output_dir(args.output)
        config_path, logger, log_file_path = annotate.generate_config(args)

        command_full_print = build_rerunnable_command(this_parser, args)
        run_snakemake(
            module_dir_name="annotate",
            module_arg_name="annotate",
            config_path=config_path,
            scripts_dir=_scripts_dir,
            args=args,
            checkamg_version=_checkamg_version,
            logger=logger,
            log_file_path=log_file_path,
            command=command_full_print,
        )

    annotate_parser.set_defaults(func=_run_annotate)
    return annotate_parser