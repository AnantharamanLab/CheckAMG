import os
from pathlib import Path

configfile: "./accessory_scripts/genomad.config.yaml"

INPUT_BASE = Path(config["input_base"]).resolve()
OUTPUT_BASE = Path(config["output_base"]).resolve()

# Pattern: <input_base>/<subset>/<sample>.fna
# subset ∈ {"viromes","metagenomes","complete_virus_genomes"}
INP_GLOB = str(INPUT_BASE / "*" / "*.fna")

def outdir(subset, sample):
    return OUTPUT_BASE / subset / sample

# Discover all .fna inputs up front
INPUT_FILES = sorted(Path(p) for p in map(str, __import__("glob").glob(INP_GLOB)))
INPUT_FILES = [p for p in INPUT_FILES if p.parent.name in {"viromes","metagenomes","complete_virus_genomes"}]
# Build expected output directories (as Snakemake "directory" targets) from discovered inputs
OUT_DIRS = [outdir(p.parent.name, p.stem) for p in INPUT_FILES]

rule all:
    input:
        [directory(p) for p in OUT_DIRS]

rule genomad:
    input:
        contigs=os.path.join(config["input_base"], "{subset}", "{sample}.fna")
    output:
        directory(os.path.join(config["output_base"], "{subset}", "{sample}"))
    threads: config.get("threads", 32)
    conda: "genomad"
    log:
        "logs/{subset}/{sample}.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})" "{output}"
        genomad end-to-end \
            --enable-score-calibration \
            --force-auto \
            --relaxed \
            --threads {threads} \
            "{input.contigs}" \
            "{output}" \
            "{config[genomad_db]}" \
            > "{log}" 2>&1
        """
