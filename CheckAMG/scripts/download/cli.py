#!/usr/bin/env python3

import argparse

from CheckAMG.scripts.download import download_dbs
from CheckAMG.scripts.common.args_formatter_logging import CustomHelpFormatter
from CheckAMG.scripts.common.args_formatter_logging import build_rerunnable_command

def add_download_subcommand(parser, subparsers, scripts_dir, default_threads, pct_total_cpu, available_memory_gb, checkamg_version):
    download_parser = subparsers.add_parser(
        "download",
        help="Download the databases required by CheckAMG.",
        description=(
            "Download the databases required by CheckAMG.\n\n"
            "About 116 GB of free disk space is required during download. This can be reduced to\n"
            "about 97 GB after downloading finishes if the human-readable HMM files are removed\n"
            "with --rm-hmm. The space required for the database can be further reduced to about \n"
            "21 GB if you do not plan on running checkamg de-novo or checkamg end-to-end and\n"
            "skip downloading the de novo database with --no-download-denovo-db."
        ),
        formatter_class=CustomHelpFormatter,
    )

    download_parser.add_argument(
        "-d",
        "--db-dir",
        dest="db_dir",
        type=str,
        required=True,
        help=(
            "Parent directory where the CheckAMG databases will be placed. The annotate and de-novo\n"
            "databases extract into their own subdirectories here. Pass this same directory to the\n"
            "'--db-dir' of 'checkamg annotate' and 'checkamg denovo'."
        ),
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
    download_parser.add_argument(
        "--download-annotate-db",
        dest="download_annotate",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Download the database required by 'checkamg annotate' (default: enabled).\n"
            "Use '--no-download-annotate-db' to skip it if it already exists."
        ),
    )
    download_parser.add_argument(
        "--download-denovo-db",
        dest="download_denovo",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Download the larger database required by 'checkamg denovo' (default: enabled).\n"
            "Use '--no-download-denovo-db' to skip it if you only need 'checkamg annotate'."
        ),
    )

    def _run_download(args, parser, this_parser=download_parser, _scripts_dir=scripts_dir, _checkamg_version=checkamg_version):
        command_full_print = build_rerunnable_command(this_parser, args)
        download_dbs.main(
            command=command_full_print,
            dest=args.db_dir,
            checkamg_version=checkamg_version,
            force=args.force,
            db_version=getattr(args, "db_version", None),
            rm_hmm=args.rm_hmm,
            download_annotate=args.download_annotate,
            download_denovo=args.download_denovo,
        )

    download_parser.set_defaults(func=_run_download)
    return download_parser

