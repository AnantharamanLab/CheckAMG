import os
from pathlib import Path

configfile: "./accessory_scripts/dramv.config.yaml"

INPUT_BASE = Path(config["input_base"]).resolve()
OUTPUT_BASE = Path(config["output_base"]).resolve()

def discover_samples():
    for subset_dir in INPUT_BASE.iterdir():
        if not subset_dir.is_dir():
            continue
        for sample_dir in subset_dir.iterdir():
            if not sample_dir.is_dir():
                continue
            fd = sample_dir / "for-dramv"
            if fd.is_dir():
                fas = list(fd.glob("*.fa"))
                tabs = list(fd.glob("*.tab"))
                if len(fas) == 1 and len(tabs) == 1:
                    yield (subset_dir.name, sample_dir.name)

SAMPLES = sorted(set(discover_samples()))
ALL_DISTILLED = [
    str(OUTPUT_BASE / subset / sample / "distilled")
    for subset, sample in SAMPLES
]

# remove complete host genomes
for i, path in enumerate(ALL_DISTILLED):
    if "host_genomes" in path:
        ALL_DISTILLED.pop(i)

rule all:
    input:
        [directory(p) for p in ALL_DISTILLED]

def fa_input(wc):
    p = INPUT_BASE / wc.subset / wc.sample / "for-dramv"
    fas = sorted(p.glob("*.fa"))
    assert len(fas) == 1, f"Expected 1 .fa in {p}, found {len(fas)}"
    return str(fas[0])

def tab_input(wc):
    p = INPUT_BASE / wc.subset / wc.sample / "for-dramv"
    tabs = sorted(p.glob("*.tab"))
    assert len(tabs) == 1, f"Expected 1 .tab in {p}, found {len(tabs)}"
    return str(tabs[0])

rule dramv_annotate:
    input:
        fa=fa_input,
        tab=tab_input
    output:
        ann_tsv=os.path.join(config["output_base"], "{subset}", "{sample}", "annotations.tsv")
    threads:
        config.get("threads", 4)
    conda:
        "DRAM"
    log:
        "logs/{subset}/{sample}_annotate.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        outdir="$(dirname {output.ann_tsv})"
        if [ -d "$outdir" ]; then rm -rf "$outdir"; fi
        DRAM-v.py annotate \
            -i "{input.fa}" \
            -v "{input.tab}" \
            -o "$outdir" \
            --keep_tmp_dir \
            --threads {threads} \
            > "{log}" 2>&1
        """

rule dramv_distill:
    input:
        ann_tsv=os.path.join(config["output_base"], "{subset}", "{sample}", "annotations.tsv")
    output:
        distilled=directory(os.path.join(config["output_base"], "{subset}", "{sample}", "distilled"))
    threads:
        config.get("threads", 4)
    conda:
        "DRAM"
    log:
        "logs/{subset}/{sample}_distill.log"
    shell:
        r"""
        mkdir -p "$(dirname {log})"
        if [ -d "{output.distilled}" ]; then rm -rf "{output.distilled}"; fi
        DRAM-v.py distill \
            -i "{input.ann_tsv}" \
            -o "{output.distilled}" \
            --max_auxiliary_score 4 \
            --log_file_path "{output.distilled}/distill.log" \
            > "{log}" 2>&1
        """
