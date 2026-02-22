#!/usr/bin/env python3

import argparse
import textwrap
import sys
import psutil
from importlib.metadata import version

from CheckAMG.scripts import CheckAMG_annotate, download_dbs
from CheckAMG.scripts.checkAMG_ASCII import ASCII

__version__ = version("checkamg")

available_memory_gb = psutil.virtual_memory().available / (1024 ** 3) # Get available memory in GB
PCT_TOTAL_CPU = 25
default_threads = max(1, int(psutil.cpu_count(logical=True) * (PCT_TOTAL_CPU/100))) # Default to PCT_TOTAL_CPU% of available threads, but at least 1

class CustomHelpFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    def _fill_text(self, text, width, indent):
        text = textwrap.dedent(text).strip()
        return "\n".join(
            textwrap.fill(line, width, initial_indent=indent, subsequent_indent=indent)
            if line and not line.lstrip().startswith("* ")
            else f"{indent}{line.strip()}"
            for line in text.splitlines()
        )

    def _split_lines(self, text, width):
        lines = []
        for part in text.splitlines():
            if not part.strip():
                lines.append("")
            else:
                lines.extend(textwrap.wrap(part, width))
        return lines


VALID_FILTER_PRESETS = {
    "default",
    "allow_glycosyl",
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

    if len(presets) > 1 and "default" in presets:
        parser.error(
            "Cannot combine 'default' with other filter presets. "
            "Use only 'default' or specify non-default combinations."
        )

    if len(presets) > 1 and "no_filter" in presets:
        parser.error("Cannot combine 'no_filter' with other filter presets. Use only 'no_filter'.")

    return list(presets)


def _add_download_subcommand(subparsers):
    download_parser = subparsers.add_parser(
        "download",
        help="Download the databases required by CheckAMG.",
        description=(
            "Download the databases required by CheckAMG.\n "
            "This requires ~40 GB of disk space (or ~21 GB finally, if '--rm-hmm' is provided)."
        ),
        formatter_class=CustomHelpFormatter,
    )

    download_parser.add_argument(
        "-d",
        "--db-dir",
        dest="db_dir",
        type=str,
        required=True,
        help="Directory where the CheckAMG database will be placed.",
    )
    download_parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        default=False,
        help="Force re-download of databases even if they already exist.",
    )
    download_parser.add_argument(
        "-r",
        "--rm-hmm",
        dest="rm_hmm",
        action="store_true",
        default=False,
        help=(
            "Remove human-readable HMM files after downloading to save space.\n"
            "CheckAMG only needs the binary files."
        ),
    )
    download_parser.add_argument(
        "-v",
        "--db-version",
        dest="db_version",
        type=str,
        default=None,
        help="Exact CheckAMG database version identifier to download (overrides latest compatible).",
    )

    def _run_download(args, parser):
        download_dbs.download_db(
            dest=args.db_dir,
            checkamg_version=__version__,
            force=args.force,
            db_version=getattr(args, "db_version", None),
        )
        if args.rm_hmm:
            download_dbs.remove_human_readable_files(dest=args.db_dir)

    download_parser.set_defaults(func=_run_download)
    return download_parser


def _add_annotate_subcommand(subparsers):
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
        help="Path to CheckAMG database files.",
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
        default=0.30,
        help="Minimum covered fraction (0.0-1.0) of HMM profiles required to report hits.",
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
        default=30,
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
        default="default",
        help=(
            "Comma-separated preset(s) controlling functional annotation filtering.\n"
            "Valid presets:\n"
            "* default (recommended)\n"
            "* allow_glycosyl (keep glycosyltransferase, glycoside-hydrolase, and related annotations)\n"
            "* allow_nucleotide (keep nucleotide metabolism annotations)\n"
            "* allow_methyl (keep methylase/methyltransferase annotations)\n"
            "* allow_lipid (keep lipopolysaccharide and phospholipid-related annotations)\n"
            "* no_filter (disable all filtering, not recommended).\n"
            "Example: --filter-presets allow_glycosyl,allow_nucleotide."
            )
    )

    outputs = annotate_parser.add_argument_group("output files")
    outputs.add_argument(
        "-kf",
        "--keep-full-hmm-results",
        dest="keep_full_hmm_results",
        action="store_true",
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
        action="store_true",
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
        help=f"Maximum number of threads allowed. Default is {PCT_TOTAL_CPU}%% of available.",
    )
    resources.add_argument(
        "-m",
        "--mem",
        dest="mem",
        type=int,
        default=round(available_memory_gb * 0.80),
        help="Max memory allowed (GB). Default is 80%% of available.",
    )
    resources.add_argument(
        "--debug",
        dest="debug",
        default=False,
        action=argparse.BooleanOptionalAction,
        help="Enable debug-level logging.",
    )

    def _run_annotate(args, parser):
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
            raw_list = ["default"]
        effective = _validate_and_resolve_filter_presets(raw_list, parser)
        args.filter_presets = ",".join(effective)

        CheckAMG_annotate.create_output_dir(args.output)
        config_path = CheckAMG_annotate.generate_config(args)
        CheckAMG_annotate.run_snakemake(config_path=config_path, args=args, checkamg_version=__version__)

    annotate_parser.set_defaults(func=_run_annotate)
    return annotate_parser


def _add_placeholder_subcommands(subparsers):
    for name, help_text in [
        ("de-novo", "(Not yet implemented) Predict auxiliary genes with an annotation-independent method."),
        ("aggregate", "(Not yet implemented) Aggregate results into a final report."),
        ("end-to-end", "(Not yet implemented) Run annotate, de-novo, and aggregate in tandem."),
    ]:
        p = subparsers.add_parser(
            name,
            help=help_text,
            description="Not yet implemented.",
            formatter_class=CustomHelpFormatter,
        )

        def _run_placeholder(args, parser, cmd=name):
            print(f"CheckAMG {cmd} functionality will be implemented here.")

        p.set_defaults(func=_run_placeholder)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "CheckAMG: Automated discovery and curation of Auxiliary Metabolic Genes (AMGs),\n"
            "          Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs)\n"
            "          encoded in viral genomes."
        ),
        formatter_class=CustomHelpFormatter,
    )
    parser.add_argument("-v", "--version", action="version", version=f"CheckAMG {__version__}")

    subparsers = parser.add_subparsers(title="modules", dest="command", required=True)

    _add_download_subcommand(subparsers)
    _add_annotate_subcommand(subparsers)
    _add_placeholder_subcommands(subparsers)

    if "--version" not in sys.argv and "-v" not in sys.argv:
        print(ASCII)
        sys.stdout.flush()
     
    args = parser.parse_args()
    args._cli_argv = sys.argv[1:] # for clean command logging
    args.func(args, parser)

if __name__ == "__main__":
    main()
