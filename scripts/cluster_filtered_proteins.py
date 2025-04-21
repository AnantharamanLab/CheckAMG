#!/usr/env/bin/python3

import os
import sys
import subprocess
import shutil
import resource
import platform
import logging

log_level = logging.DEBUG if snakemake.params.debug else logging.INFO
log_file = snakemake.params.log
logging.basicConfig(
    level=log_level,
    format='%(asctime)s | %(levelname)s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    handlers=[
        logging.FileHandler(log_file, mode='a'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger()

print("========================================================================\n       Step 10/22: Align the filtered proteins against each other        \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n       Step 10/22: Align the filtered proteins against each other        \n========================================================================\n")

def set_memory_limit(limit_in_gb):
    current_os = platform.system()
    if current_os == "Linux":
        limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))

def run_mmseqs_command(command, log_file):
    """Run mmseqs command and log output"""
    with open(log_file, "a") as log:
        result = subprocess.run(
            command,
            stdout=log, 
            stderr=subprocess.PIPE, 
            text=True
        )
        if result.returncode != 0:
            log.write(result.stderr)
            log.flush()
            raise subprocess.CalledProcessError(result.returncode, result.args, output=result.stdout, stderr=result.stderr)

def main():
    input_prot_subdir = snakemake.params.input_prot_subdir
    output = snakemake.params.protein_db
    out_db_dir = snakemake.params.out_db_dir
    search_db_dir = snakemake.params.search_parent
    search_db = snakemake.params.search_db
    search_sensitivity = snakemake.params.search_sensitivity
    cov_fraction = snakemake.params.cov_fraction
    search_output = snakemake.params.search_output
    wdir = snakemake.params.wdir
    threads = snakemake.threads
    mem_limit = snakemake.resources.mem
    set_memory_limit(mem_limit)
    
    logger.info("MMseqs2 clustering of filtered proteins starting...")
    
    # Gather input protein files
    if os.path.exists(os.path.join(input_prot_subdir, "vMAG_proteins")) and any(f.endswith('.faa') for f in os.listdir(os.path.join(input_prot_subdir, "vMAG_proteins")) if os.path.isfile(os.path.join(os.path.join(input_prot_subdir, "vMAG_proteins"), f))):
        input_files = [os.path.join(input_prot_subdir, "single_contig_proteins.faa")] + [os.path.join(input_prot_subdir, "vMAG_proteins", faa) for faa in os.listdir(os.path.join(input_prot_subdir, "vMAG_proteins")) if faa.endswith(".faa")]
    else:
        input_files = [os.path.join(input_prot_subdir, "single_contig_proteins.faa")]

    # Step 1: Create MMseqs2 database
    log_file_createdb = os.path.join(wdir, "mmseqs_createdb.log")
    os.makedirs(wdir, exist_ok=True)
    os.makedirs(out_db_dir, exist_ok=True)

    logger.info("Creating MMseqs2 database from input protein files...")
    cmd_createdb = ["mmseqs", "createdb", *input_files, output]
    run_mmseqs_command(cmd_createdb, log_file_createdb)

    # Step 2: MMseqs2 search
    search_log_file = os.path.join(wdir, "mmseqs_search.log")
    logger.info(f"Executing mmseqs search with a maximum sensitivity {search_sensitivity} and a minimum coverage of {cov_fraction}...")
    os.makedirs(search_db_dir, exist_ok=True)
    cmd_search = [
        "mmseqs", "search", output, output, search_db, wdir,
        # "--start-sens", "1", # Faster, less sensitive search, leads to fewer AAI connections and thus more singleton clusters
        "--start-sens", str(round(search_sensitivity * 0.5)), # Slightly slower, slightly more sensitive search, but still faster than just using the maximum sensitivity
        "--sens-steps", "3",
        "-s", str(search_sensitivity),
        "--min-seq-id", "0",
        "-c", str(cov_fraction),
        "--threads", str(threads)
    ]
    run_mmseqs_command(cmd_search, search_log_file)

    # Step 3: Convert MMseqs2 search output to TSV
    convert_log_file = os.path.join(wdir, "mmseqs_convertalis.log")
    logger.info("Converting MMseqs2 search output to TSV format...")
    cmd_convertalis = [
        "mmseqs", "convertalis", output, output, search_db, search_output,
        "--format-output", "query,target,fident,alnlen,qlen,tlen,qcov,tcov,bits",
        "--threads", str(threads)
    ]
    run_mmseqs_command(cmd_convertalis, convert_log_file)
    
    # Step 4: Remove the intermediate search files
    logger.debug(f"Removing mmseqs search files at {search_db_dir}...")
    shutil.rmtree(search_db_dir)
    for item in os.listdir(wdir):
        if not item.startswith("protein_search_output.tsv"):
            if os.path.islink(os.path.abspath(os.path.join(wdir, item))):
                os.unlink(os.path.abspath(os.path.join(wdir, item)))
            elif os.path.isdir(os.path.abspath(os.path.join(wdir, item))):
                shutil.rmtree(os.path.abspath(os.path.join(wdir, item)))
            else:
                os.remove(os.path.abspath(os.path.join(wdir, item)))

    # Step 5: Remove the protein search database
    logger.debug(f"Removing protein search database directory {out_db_dir}...")
    shutil.rmtree(out_db_dir)
    
    logger.info("MMseqs2 clustering of filtered proteins completed.")

if __name__ == "__main__":
    main()
