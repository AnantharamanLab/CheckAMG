# CheckAMG
[![PyPI](https://img.shields.io/pypi/v/checkamg)](https://pypi.org/project/checkamg/)
[![CheckAMG DB v1.1.1](https://img.shields.io/badge/CheckAMG%20DB-v1.1.1-blue)](https://zenodo.org/records/21776005)

**Automated discovery and curation of Auxiliary Metabolic Genes (AMGs), Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) encoded by viral genomes**

> ⚠️ **This tool is in active development and has not yet been peer-reviewed.**

CheckAMG identifies and curates high-confidence [auxiliary viral genes (AVGs: AMGs, AReGs, and APGs)](https://doi.org/10.1038/s41564-025-02095-4) in viral genomes. It combines two complementary approaches: a reference-based path that uses curated functional annotations and viral genome context, and a reference-independent path that uses a protein genome language model (Protein Set Transformer, PST). Its prediction approach reflects years of community-defined standards for identifying auxiliary genes, validating that they are virus-encoded, and filtering common misannotations.

CheckAMG supports:

* Nucleotide or protein input
* Single-contig viral genomes or vMAGs (multi-contig)
* Running on viral genomes or viromes/metagenomes directly
* CPU-only or GPU-accelerated execution

Full documentation, detailed installation instructions, parameter descriptions, an FAQ, and notebooks to reproduce the analyses in the manuscript are locaed in the **[CheckAMG Wiki](https://github.com/AnantharamanLab/CheckAMG/wiki)**.

## Installation

CheckAMG depends on PyTorch and the source-built `torch_scatter` extension (via PST), which must be installed before CheckAMG so the correct CPU/CUDA build is selected. The short version:

```bash
# Create a minimal conda environment
conda create -y -n CheckAMG python=3.11
conda activate CheckAMG
python -m ensurepip --upgrade
python -m pip install --upgrade pip setuptools wheel packaging uv

# Install torch + PyG extensions + FAISS (CPU example; see Wiki for GPU/CUDA)
uv pip install torch==2.8.0 --index-url https://download.pytorch.org/whl/cpu
uv pip install torch_geometric
uv cache clean torch_scatter 2>/dev/null || true
uv pip install torch_scatter --no-build-isolation --no-binary torch_scatter --refresh
uv pip install faiss-cpu

# Install CheckAMG (pulls PST and the remaining ML dependencies)
uv pip install checkamg
```

GPU builds, ABI/`torch_scatter` troubleshooting, the full dependency list, and `conda-pack` instructions are in the **[Installation guide](https://github.com/AnantharamanLab/CheckAMG/wiki/Installation)**.

Download the databases (~116 GB, reducible to ~21 GB with `--rm-hmm` and excluding the de novo module's database):

```bash
checkamg download -d /path/to/db/destination --rm-hmm
```

## Modules

CheckAMG predicts AVGs through two independent approaches that can be combined:

* **`annotate`** (annotation-based): predicts and curates AVGs from functional annotations (profile HMM searches) and viral genome context. It classifies predictions as metabolic (AMG), physiological (APG), or regulatory (AReG) and assigns a viral-origin confidence. This path is conservative and produces curated, trustworthy predictions on its own.
* **`de-novo`** (annotation-independent): embeds query proteins with the PST genome language model and predicts, by nearest-neighbor search against a labeled reference, whether each protein is virus-encoded and auxiliary. It does not assign specific functions, but it is **more sensitive** than `annotate` and can recover AVGs that have no informative annotation, which `annotate` cannot detect.
* **`aggregate`**: merges `annotate` and `de-novo` results into a single report.
* **`end-to-end`**: runs `annotate`, then `de-novo` (reusing the proteins from `annotate`), then `aggregate`, into one organized output.
* **`download`**: fetches the required databases.
* **`train`**: finetunes a custom CheckAMG-PST model for `de-novo`.

**We recommend running `end-to-end` by default when feasible.** Because `de-novo` is more sensitive than `annotate` alone, the combined run often identifies additional AVGs that annotation cannot. `de-novo` benefits substantially from a GPU, but GPU access is not required. If a GPU is unavailable or `de-novo` is impractical for your input size, running `annotate` alone is perfectly fine and will still produce curated, trustworthy AVG predictions.

## Quick start

Example data to test your installation are provided in [`examples/example_data`](https://github.com/AnantharamanLab/CheckAMG/tree/main/examples/example_data).

End-to-end (recommended):

```bash
checkamg end-to-end \
  -d /path/to/db/destination \
  -i examples/example_data/single_contig_viruses.fasta \
  -I examples/example_data/multi_contig_vMAGs \
  -o CheckAMG_example_out
```

Annotate only:

```bash
checkamg annotate \
  -d /path/to/db/destination \
  -i examples/example_data/single_contig_viruses.fasta \
  -I examples/example_data/multi_contig_vMAGs \
  -o CheckAMG_annotate_out
```

Run `checkamg -h` or `checkamg <module> -h` for all options, and see the **[Wiki](https://github.com/AnantharamanLab/CheckAMG/wiki)** for detailed descriptions of every parameter and output.

## Error reporting

To report bugs or request features, please use the [GitHub Issues page](https://github.com/AnantharamanLab/CheckAMG/issues). To help diagnose issues and debug, please re-run your command that led to the error with the `--debug` flag, and provide all log files (including module-specific logs within the `01_annotate`, `02_denovo`, and `03_aggregate` outputs if running `checkamg end-to-end`).

## Citation

A manuscript describing CheckAMG is in preparation. Until then, please cite the GitHub repository and the [CheckAMG database on Zenodo](https://zenodo.org/records/18407279).

Authors:

* James C. Kosmopoulos (**kosmopoulos [at] wisc [dot] edu**)
* Cody Martin
* Karthik Anantharaman (**karthik [at] bact [dot] wisc [dot] edu**)
