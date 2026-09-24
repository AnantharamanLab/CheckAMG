import os
from pathlib import Path

configfile: "./accessory_scripts/checkamg_annotate.config.yaml"

INPUT_BASE = Path(config["input_base"]).resolve()
OUTPUT_BASE = Path(config["output_base"]).resolve()

# Pattern: <input_base>/<subset>/<sample>.fna
# subset ∈ {"viromes","metagenomes","complete_virus_genomes"}
INP_GLOB = str(INPUT_BASE / "*" / "*" / "genes_reformatted.faa")

def outdir(subset, sample):
    return OUTPUT_BASE / subset / sample

# Discover all .faa inputs up front
INPUT_FILES = sorted(Path(p) for p in map(str, __import__("glob").glob(INP_GLOB)))
INPUT_FILES = [p for p in INPUT_FILES if p.parent.parent.name in {"viromes","metagenomes","complete_virus_genomes"}]
# Build expected output directories (as Snakemake "directory" targets) from discovered inputs
OUT_DIRS = [outdir(p.parent.parent.name, p.parent.name) for p in INPUT_FILES]

rule all:
    input:
        [directory(p) for p in OUT_DIRS]

rule checkamg_annotate:
    input:
        ptns=os.path.join(config["input_base"], "{subset}", "{sample}", "genes_reformatted.faa")
    output:
        directory(os.path.join(config["output_base"], "{subset}", "{sample}"))
    threads: config.get("threads", 16)
    # conda: "checkamg_pypi_test"
    conda: "CheckAMG-dev"
    log:
        "logs/{subset}/{sample}.log"
    params:
        tmpdir=directory(os.path.join(config["output_base"], "{subset}", "{sample}_tmp"))
    shell:
        r"""
        mkdir -p "$(dirname {log})" "{output}"
        checkamg annotate \
            -p "{input.ptns}" \
            --input-type prot \
            -o "{output}" \
            --keep-full-hmm-results \
            --save-to-parquet \
            -t {threads} \
            -d "{config[checkamg_db]}" \
            --mem 300 \
            > "{log}" 2>&1
        """
