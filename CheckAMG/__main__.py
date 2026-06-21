#!/usr/bin/env python3

import argparse
import os
import sys
import shlex
import psutil
from importlib.metadata import version

from CheckAMG.scripts.common.checkAMG_ASCII import ASCII
from CheckAMG.scripts.common.args_formatter_logging import CustomHelpFormatter
from CheckAMG.scripts.download.cli import add_download_subcommand
from CheckAMG.scripts.annotate.cli import add_annotate_subcommand
from CheckAMG.scripts.denovo.cli import add_denovo_subcommand
from CheckAMG.scripts.train_pst.cli import add_train_subcommand
from CheckAMG.scripts.aggregate.cli import add_aggregate_subcommand
from CheckAMG.scripts.end_to_end.cli import add_end_to_end_subcommand

__version__ = version("checkamg")

AVAIL_MEM_GB = psutil.virtual_memory().available / (1024 ** 3) # Get available memory in GB
PCT_TOTAL_CPU = 25
DEFAULT_THREADS = max(1, int(psutil.cpu_count(logical=True) * (PCT_TOTAL_CPU/100))) # Default to PCT_TOTAL_CPU% of available threads, but at least 1


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

    add_end_to_end_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
        )
    
    add_annotate_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
        )
    
    add_denovo_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
        )
    
    add_aggregate_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
    )
    
    add_download_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
        )

    add_train_subcommand(
        parser,
        subparsers,
        scripts_dir=os.path.abspath(os.path.dirname(__file__)),
        default_threads=DEFAULT_THREADS,
        pct_total_cpu=PCT_TOTAL_CPU,
        available_memory_gb=AVAIL_MEM_GB,
        checkamg_version=__version__
        )

    if "--version" not in sys.argv and "-v" not in sys.argv:
        print(ASCII)
        sys.stdout.flush()
     
    args = parser.parse_args()
    args._cli_argv = [shlex.quote(x) for x in sys.argv[1:]]
    args.func(args, parser)

if __name__ == "__main__":
    main()
