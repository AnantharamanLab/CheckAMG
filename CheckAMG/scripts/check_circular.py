#!/usr/bin/env python3

# This script is modified from portions of the tool CheckV, originally developed by
# S. Nayfach et al. (2021) in their publication "CheckV: assessing the quality of metagenome-assembled viral genomes" (Nature Biotechnology, https://doi.org/10.1038/s41587-020-00774-7). 
# The original source code can be found at: 
# https://bitbucket.org/berkeleylab/checkv/src/master
# 
# Specifically, code modifications were made based on:
# 1. utility.py - https://bitbucket.org/berkeleylab/checkv/src/master/checkv/utility.py
# 2. complete_genomes.py - https://bitbucket.org/berkeleylab/checkv/src/master/checkv/modules/complete_genomes.py
#
# These modifications allow for assessing sequence completeness based on direct terminal repeats (DTRs)
# and inverted terminal repeats (ITRs) without having to run CheckV in its entirety.
#
# *** License Agreement ***
# CheckV Copyright (c) 2020, The Regents of the University of California, 
# through Lawrence Berkeley National Laboratory (subject to receipt of
# any required approvals from the U.S. Dept. of Energy). All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# (1) Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# (2) Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# (3) Neither the name of the University of California, Lawrence Berkeley
# National Laboratory, U.S. Dept. of Energy nor the names of its contributors
# may be used to endorse or promote products derived from this software
# without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# You are under no obligation whatsoever to provide any bug fixes, patches,
# or upgrades to the features, functionality or performance of the source
# code ("Enhancements") to anyone; however, if you choose to make your
# Enhancements available either publicly, or directly to Lawrence Berkeley
# National Laboratory, without imposing a separate written license agreement
# for such Enhancements, then you hereby grant the following license: a
# non-exclusive, royalty-free perpetual license to install, use, modify,
# prepare derivative works, incorporate into other computer software,
# distribute, and sublicense such enhancements or derivative works thereof,
# in binary and source code form.

import os
import sys
import logging
import re
import collections
from multiprocessing import Pool
from tqdm import tqdm
import resource
import platform
from pyfastatools import Parser

def set_memory_limit(limit_in_gb):
    limit_in_bytes = limit_in_gb * 1024 * 1024 * 1024
    try:
        resource.setrlimit(resource.RLIMIT_AS, (limit_in_bytes, limit_in_bytes))
    except (ValueError, OSError, AttributeError) as e:
        logger.warning(f"Unable to set memory limit. Error: {e}")
    
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

print("========================================================================\n     Step 2/11: Check the circularity of the input genome sequences    \n========================================================================")
with open(log_file, "a") as log:
    log.write("========================================================================\n     Step 2/11: Check the circularity of the input genome sequences    \n========================================================================\n")

# Genome class
class Genome:
    def __init__(self):
        self.seq = None
        self.id = None
        self.length = None
        self.tr = None
        self.flagged = False
        self.reason = None

# TR class
class TR:
    def __init__(self):
        pass

def fetch_dtr(fullseq, min_length=20):
    startseq = fullseq[0:min_length]
    matches = [
        m.start() for m in re.finditer("(?={0})".format(re.escape(startseq)), fullseq)
    ]
    matches = [_ for _ in matches if _ >= len(fullseq) / 2]
    for matchpos in matches:
        endseq = fullseq[matchpos:]
        if fullseq[0 : len(endseq)] == endseq:
            return endseq
    return ""

def reverse_complement(seq):
    trans = str.maketrans("ACTG", "TGAC")
    return seq[::-1].translate(trans)

def fetch_itr(seq, min_len, max_len=1000):
    rev = reverse_complement(seq)
    if seq[:min_len] == rev[:min_len]:
        i = min_len + 1
        while seq[:i] == rev[:i] and i <= max_len:
            i += 1
        return seq[: i - 1]
    return ""

def calculate_kmer_frequency(genome_seq, k=15):
    """Calculate kmer frequency for a genome sequence."""
    kmer_counts = collections.Counter()
    for i in range(len(genome_seq) - k + 1):
        kmer = genome_seq[i:i+k]
        kmer_counts[kmer] += 1
    avg_kmer_freq = sum(kmer_counts.values()) / len(kmer_counts)
    return avg_kmer_freq

def process_genomes_batch(genome_records, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq):
    """Process a batch of genomes."""
    processed_genomes = {}
    for record in genome_records:
        genome = Genome()
        genome.id = record[0]
        genome.seq = str(record[1]).upper()
        genome.length = len(genome.seq)
        genome.kmer_freq = calculate_kmer_frequency(genome.seq, k)
        
        genome.tr = TR()
        dtr = fetch_dtr(genome.seq, tr_min_len)
        itr = fetch_itr(genome.seq, tr_min_len, tr_max_len)
        if len(dtr) < tr_min_len and len(itr) < tr_min_len:
            genome.tr.type = None
        elif len(dtr) >= len(itr):
            genome.tr.type = "DTR"
            genome.tr.seq = dtr
            genome.tr.length = len(dtr)
        else:
            genome.tr.type = "ITR"
            genome.tr.seq = itr
            genome.tr.length = len(itr)

        # Filter terminal repeats and calculate flags
        if genome.tr.type is not None:
            mode_base, mode_count = collections.Counter(genome.tr.seq).most_common(1)[0]
            genome.tr.mode_freq = 1.0 * mode_count / len(genome.tr.seq)
            genome.tr.n_freq = 1.0 * genome.tr.seq.count("N") / len(genome.tr.seq)
            genome.tr.count = genome.seq.count(genome.tr.seq)
            
            # Flag based on user-defined thresholds
            if genome.tr.n_freq > tr_max_ambig:
                genome.flagged = True
                genome.reason = "Too many ambiguous bases in TR"
            elif genome.tr.count > tr_max_count:
                genome.flagged = True
                genome.reason = "Repetitive TR sequence"
            elif genome.tr.mode_freq > tr_max_basefreq:
                genome.flagged = True
                genome.reason = "Low complexity TR"
            else:
                genome.flagged = False
            
        # Check for kmer frequency
        if genome.kmer_freq > kmer_max_freq:
            genome.flagged = True
            genome.reason = "Multiple genome copies detected"

        processed_genomes[genome.id] = genome
    return processed_genomes

def process_genomes_batch_wrapper(batch, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq):
    combined = {}
    for genome_records in batch:
        result = process_genomes_batch(genome_records, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq)
        combined.update(result)
    return combined

def parallel_processing(single_contig_fasta, input_files, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq, num_workers):
    genome_units = []
    contig_to_genome = {}

    for input_file in input_files:
        if not os.path.isfile(input_file):
            logger.error(f"Input file not found: {input_file}")
            raise FileNotFoundError(f"Input file not found: {input_file}")
        records = [(record.header.name, record.seq) for record in Parser(input_file)]
        if not records:
            logger.error(f"No sequences parsed from: {input_file}. Is it empty or not in FASTA format?")
            raise ValueError(f"No sequences parsed from: {input_file}. Is it empty or not in FASTA format?")
        logger.debug(f"{input_file}: {len(records)} records")

        if input_file == snakemake.params.input_single_contig_genomes:
            # Split by contig ID
            contig_groups = {}
            for name, seq in records:
                contig_to_genome[name] = name
                contig_id = name.rsplit("_", 1)[0]
                contig_groups.setdefault(contig_id, []).append((name, seq))
            logger.debug(f"{input_file}: {len(contig_groups)} contig groups (treated as genomes)")
            genome_units.extend(contig_groups.values())
        else:
            # Treat each vMAG file as one genome unit
            for record in records:
                contig_to_genome[record[0]] = input_file
            genome_units.append(records)

    logger.info(f"Total sequences to check: {len(contig_to_genome):,} ({len(genome_units):,} genomes)")

    # Chunk into processing batches
    chunk_size = max(1, len(genome_units) // num_workers)
    genome_chunks = [
        genome_units[i:i + chunk_size] for i in range(0, len(genome_units), chunk_size)
    ]

    with Pool(processes=num_workers) as pool:
        async_results = [
            pool.apply_async(
                process_genomes_batch_wrapper,
                args=(chunk, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq)
            )
            for chunk in genome_chunks
        ]

        results = []
        for r in tqdm(async_results, desc="Checking circularity", unit="chunk", total=len(async_results)):
            results.append(r.get())

    # Combine all results
    combined_genomes = {}
    for result in results:
        combined_genomes.update(result)

    return combined_genomes, contig_to_genome

def main():
    input_fasta = snakemake.params.input_single_contig_genomes
    input_vmag_fastas = snakemake.params.input_vmag_fastas
    output_path = snakemake.params.circularity_tbl
    tr_min_len = snakemake.params.tr_min_len
    tr_max_len = snakemake.params.tr_max_len
    tr_max_count = snakemake.params.tr_max_count
    tr_max_ambig = snakemake.params.tr_max_ambig
    tr_max_basefreq = snakemake.params.tr_max_basefreq
    kmer_max_freq = snakemake.params.kmer_max_freq
    k = snakemake.params.k
    mem_limit = snakemake.resources.mem
    num_workers = snakemake.threads
    set_memory_limit(mem_limit)

    logger.info("Sequence circularity check starting...")
    
    if os.path.exists(input_vmag_fastas):
        vmag_fasta_files = [os.path.join(input_vmag_fastas, fna) for fna in os.listdir(input_vmag_fastas) if fna.endswith('.fna')]
        logger.debug(f"vMAG fna files found: {vmag_fasta_files}")
        if os.path.exists(input_fasta):
            logger.debug(f"Single contig genomes file found: {input_fasta}")
            input_files = [input_fasta] + vmag_fasta_files
        else:
            logger.debug("No single contig genomes file found.")
            input_files = vmag_fasta_files
    else:
        vmag_fasta_files = None
        logger.debug("No vMAG fna files found.")
        if os.path.exists(input_fasta):
            logger.debug(f"Single contig genomes file found: {input_fasta}")
            input_files = [input_fasta]
        else:
            logger.error("No single contig genomes file found.")
            raise FileNotFoundError("No input single-contig virus genome or vMAG files found.")

    # Process genomes in parallel
    sequences, contig_to_genome = parallel_processing(input_fasta, input_files, k, tr_min_len, tr_max_len, tr_max_count, tr_max_ambig, tr_max_basefreq, kmer_max_freq, num_workers)
    
    # Write results to output
    with open(output_path, "w") as out:
        header = [
            "contig", "contig_length", "kmer_freq", "prediction_type",
            "repeat_length", "repeat_count", "repeat_n_freq", 
            "repeat_mode_base_freq", "repeat_seq"
        ]
        out.write("\t".join(header) + "\n")
        for genome in sequences.values():
            if genome.tr.type is not None:
                row = [
                    genome.id, genome.length, genome.kmer_freq, genome.tr.type,
                    genome.tr.length, genome.tr.count, genome.tr.n_freq, 
                    genome.tr.mode_freq, genome.tr.seq
                ]
                out.write("\t".join(map(str, row)) + "\n")

    logger.info(f"Number of sequences checked: {len(sequences):,} ({len(set(contig_to_genome.values())):,} genomes)")
    logger.info(f"Number of circular sequences detected: {len([g for g in sequences.values() if g.tr.type is not None]):,}")
    logger.info("Circularity check completed.")

if __name__ == "__main__":
    main()