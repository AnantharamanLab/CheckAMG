# CheckAMG
[![PyPI](https://img.shields.io/pypi/v/checkamg)](https://pypi.org/project/checkamg/)

**Automated identification and curation of Auxiliary Metabolic Genes (AMGs), Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) in viral genomes and metagenomes**

> ⚠️ **This tool is in active development and has not yet been peer-reviewed.**

## Overview

CheckAMG is a pipeline for high-confidence identification and curation of auxiliary genes (AMGs, AReGs, APGs) in viral genomes. It leverages functional annotations, genomic context, and manually curated lists of AVG annotations. Its prediction approach reflects years of community-defined standards for identifying and auxiliary genes, validating that they are virus-encoded, and the filtering of common misannotations.

CheckAMG supports:

* Nucleotide or protein input
* Single-contig viral genomes or vMAGs (multi-contig)
* Sequences from viral genomes or metagenomes

## Dependencies

See `pyproject.toml` for all dependencies. Major packages:

* `python >=3.11, <3.13`
* [`metapyrodigal>=1.4.1`](https://github.com/cody-mar10/metapyrodigal)
* `polars-u64-idx>=1.30.0`
* [`pyfastatools==2.5.0`](https://github.com/cody-mar10/pyfastatools)
* [`pyhmmer==0.11.1`](https://github.com/lukas-schillinger/pyhmer)
* `sklearn==1.5.0`
* `snakemake==8.23.2`

## Installation
### Step 1 (reccomended): create a conda environment

```bash
conda create -n CheckAMG python=3.11
conda activate CheckAMG
```

### Step 2: Install from PyPI

```bash
pip install checkamg
```

### From bioconda

*Coming soon.*

## Quick start

Example data to test your installation of CheckAMG are provided in the `example_data` folder of this repository.

```
checkamg download -d /path/to/db_dir

checkamg annotate \
  -d /path/to/db_dir \
  -g example_data/single_contig_viruses.fasta \
  -vg example_data/multi_contig_vMAGs \
  -o CheckAMG_example_out
```

## Usage

CheckAMG has multiple modules. The main modules that will be used for AVG predictions are `annotate`, `de-novo`, and `end-to-end`. Currently, only the `annotate` module has been implemented, and its associated `download` module to download its required databases.

Run `checkamg -h` for full options and module descriptions:

```
usage: checkamg [-h] [-v] {download,annotate,de-novo,aggregate,end-to-end} ...

CheckAMG: automated identification and curation of Auxiliary Metabolic Genes (AMGs), Auxiliary
Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) in viral genomes.

positional arguments:
  {download,annotate,de-novo,aggregate,end-to-end}
                        CheckAMG modules
    download            Download the databases required by CheckAMG.
    annotate            Predict and curate auxiliary genes in viral genomes based on functional
                        annotations and genomic context.
    de-novo             (Not yet implemented) Predict auxiliary genes with an annotation-
                        independent method using a protein-based genome language model (Protein
                        Set Transformer).
    aggregate           (Not yet implemented) Aggregate the results of the CheckAMG annotate and
                        de-novo modules to produce a final report of auxiliary gene predictions.
    end-to-end          (Not yet implemented) Executes CheckAMG annotate, de-novo, and aggregate
                        in tandem.

options:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
```

### CheckAMG annotate

The `annotate` module is for the automated prediction and curation of auxiliary genes in viral genomes based on functional annotations and genomic context.

**Basic usage:**

```
checkamg annotate -g <genomes.fna> -d <db_dir> -o <output_dir>
```

Basic arguments:

* `-g`, `--genomes`: Path to single-contig nucleotide genomes in a single FASTA file
* `-vg`, `--vmags`: Path to a folder containing multi-contig nucleotide genomes
* `-p`, `--proteins`: Path to single-contig amino acid input in a single FASTA file
* `-vp`, `--vmag_proteins`: Path to a folder containing multi-contig amino acid input
* `-d`, `--db_dir`: Path to CheckAMG databases download with `checkamg download`
* `-o`, `--output`: Path to the output directory where CheckAMG intermediates and results will be written

Notes:

* At least one of `--genomes` or `--vmags`, or one of `--proteins` or `--vmag_proteins`, must be provided
* Both nucleotide and protein input types cannot be mixed
* Providing single-contig genomes or vMAGs only affects the labeling and organization of results, and does not affect AVG predictions
* Protein headers must be in [prodigal format](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output#protein-translations) (e.g. `>Contig1_1 # 144 # 635 # 1` or `>Contig1_2 # 1535 # 635 # -1`)

**Full usage:**

```
usage: checkamg annotate [-h] -d DB_DIR -o OUTPUT [-g GENOMES] [-vg VMAGS] [-p PROTEINS]
                         [-vp VMAG_PROTEINS] [--input_type INPUT_TYPE] [-l MIN_LEN] [-f MIN_ORF]
                         [-n MIN_ANNOT] [-c COV_FRACTION] [-Z WINDOW_SIZE] [-F MAX_FLANK]
                         [-V MIN_FLANK_VSCORE] [-H | --use_hallmark | --no-use_hallmark]
                         [-t THREADS] [-m MEM] [--debug | --no-debug]

Predict and curate auxiliary genes in viral genomes based on functional annotations and genomic
context.

options:
  -h, --help            show this help message and exit
  --input_type INPUT_TYPE
                        Specifies whether the input files are nucleotide genomes (nucl) or
                        translated amino-acid genomes (prot). Providing proteins as input will
                        skip the pyrodigal-gv step, but it will be unable to tell whether viral
                        genomes are circular, potentially losing additional evidence for
                        verifying the viral origin of putative auxiliary genes. (default: nucl).
  -l MIN_LEN, --min_len MIN_LEN
                        Minimum length in base pairs for input sequences (default: 5000).
  -f MIN_ORF, --min_orf MIN_ORF
                        Minimum number of open reading frames (proteins) inferred by pyrodigal-gv
                        for input sequences (default: 4).
  -n MIN_ANNOT, --min_annot MIN_ANNOT
                        Minimum percentage (0.0-1.0) of genes in a genome/contig required to have
                        been assigned a functional annotation using the CheckAMG database to be
                        considered for contextual analysis. (default: 0.2).
  -c COV_FRACTION, --cov_fraction COV_FRACTION
                        Minimum covered fraction for HMM alignments (default: 0.5).
  -Z WINDOW_SIZE, --window_size WINDOW_SIZE
                        Size in base pairs of the window used to calculate the average VL-score
                        of genes on a contig (default: 25000).
  -F MAX_FLANK, --max_flank MAX_FLANK
                        Maximum length in base pairs to check on the left/right flanks of
                        potentially auxiliary genes when checking for virus-like genes and non-
                        virus-like genes (default: 5000).
  -V MIN_FLANK_VSCORE, --min_flank_Vscore MIN_FLANK_VSCORE
                        Minimum V-score of genes in flanking regions required to verify a
                        potential auxiliary gene as viral and not host sequence contamination
                        (0.0-10.0) (default: 10.0).
  -H, --use_hallmark, --no-use_hallmark
                        Use viral hallmark gene annotations instead of V-scores when chekcing
                        flanking regions of potential auxiliary genes for viral verification
                        (default: False).
  -t THREADS, --threads THREADS
                        Number of threads to use for pyrodigal-gv and pyhmmer (default: 10).
  -m MEM, --mem MEM     Maximum amount of memory allowed to be allocated in GB (default: 80% of
                        available).
  --debug, --no-debug   Log CheckAMG genome with debug-level detail (default: False).

required arguments:
  -d DB_DIR, --db_dir DB_DIR
                        Path to CheckAMG database files (Required). (default: None)
  -o OUTPUT, --output OUTPUT
                        Output directory for all generated files and folders (Required).
                        (default: None)
  -g GENOMES, --genomes GENOMES
                        Input viral genome(s) in nucleotide fasta format (.fna or .fasta).
                        Expectation is that individual virus genomes are single contigs.
                        (default: None)
  -vg VMAGS, --vmags VMAGS
                        Path to folder containing vMAGs (multiple contigs) rather than single-
                        contig viral genomes. Expectation is that the folder contains one .fna or
                        .fasta file per virus genome and that each genome contains multiple
                        contigs. (default: None)
  -p PROTEINS, --proteins PROTEINS
                        Input viral genome(s) in amino-acid fasta format (.faa or .fasta).
                        Required if --input_type is prot. Expectations are that the amino-acid
                        sequence headers are in Prodigal format (>[CONTIG NAME]_[CDS NUMBER] #
                        START # END # FRAME # ...) and that each contig encoding proteins
                        represents a single virus genome. (default: None)
  -vp VMAG_PROTEINS, --vmag_proteins VMAG_PROTEINS
                        Path to folder containing vMAGs (multiple contigs) in amino-acid fasta
                        format (.faa or .fasta) rather than single-contig viral genomes.
                        Expectation is that the folder contains one .faa or .fasta file per virus
                        genome and that each genome file contains amino-acid sequences encoded on
                        multiple contigs. Required if --input_type is 'prot'. (default: None)
```

#### Outputs

The CheckAMG annotate output folder will have the following structure:

```
CheckAMG_annotate_output
├── CheckAMG_annotate.log
├── config_annotate.yaml
├── results/
│   ├── faa_metabolic/
│   │   ├── AMGs_all.faa
│   │   ├── AMGs_high_confidence.faa
│   │   ├── AMGs_low_confidence.faa
│   │   └── AMGs_medium_confidence.faa
│   ├── faa_physiology/
│   │   ├── APGs_all.faa
│   │   ├── APGs_high_confidence.faa
│   │   ├── APGs_low_confidence.faa
│   │   └── APGs_medium_confidence.faa
│   ├── faa_regulatory/
│   │   ├── AReGs_all.faa
│   │   ├── AReGs_high_confidence.faa
│   │   ├── AReGs_low_confidence.faa
│   │   └── AReGs_medium_confidence.faa
│   ├── final_results.tsv
│   ├── gene_annotations.tsv
│   ├── genes_genomic_context.tsv
│   ├── metabolic_genes_curated.tsv
│   ├── physiology_genes_curated.tsv
│   └── regulation_genes_curated.tsv
└── wdir/
```

* `CheckAMG_annotate.log`: Log file for the CheckAMG annotate run
* `config_annotate.yaml`: Snakemake pipeline configuration
* `results/`: Main results directory
    * `faa_metabolic/`, `faa_physiology/`, `faa_regulatory/`: Predicted AVGs by type and confidence
    * `final_results.tsv`: Summary table of AVG predictions
    * `gene_annotations.tsv`: All gene annotations
    * `genes_genomic_context.tsv`: Gene-level genomic context for confidence assignment
    * `*_genes_curated.tsv`: Curated lists of metabolic, physiological, and regulatory genes after filtering false positives
* `wdir/`: Intermediate files

Examples of these output files are provided in the `example_outputs` folder of this repository.

### CheckAMG de-novo
Coming soon.

### CheckAMG end-to-end
Coming soon.

## Important Notes / FAQs
### 1. What is an *AVG*?
An AVG is an **A**uxiliary **V**iral **G**ene, a virus-encoded gene that is non-essential for viral replication but augments host metabolism (AMGs), physiology(APGs), or regulation (AReGs). In the past, all AVGs have been referred to as AMGs, but recently the term AVG has been adopted to include broader host-modulating functions, not just metabolism.

Examples:
* A virus-encoded *psbA* or *soxY* would be an AMG because they encode proteins with functions in host photosynthesis and sulfide oxidation
* A virus-encoded *VasG* type VI secretion system protein or *HicA* toxin would be an APG because they are involved in host physiology
* A *LuxR* transcriptional regulator or an `AsiA` anti-sigma factor protein would be an AReG because they are likely involved in the regulation of host gene expression

Despite the name "CheckAMG", this tool also predicts APGs and AReGs using the same pipeline, differing only by functional annotation criteria.

### 2. How does CheckAMG classify and curate its predictions?

CheckAMG applies a two-stage filtering process:

1. Use a list of curated profile HMMs that represent [metabolic](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AMGs.tsv), [physiological](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/APGs.tsv), and [regulatory](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AReGs.tsv) genes to come up with initial AVG candidates
2. Use a second list of curated keywords/substrings that will be used to filter 'false' [AMGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_amgs.csv), [APGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_apgs.csv), and [AReGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_aregs.csv)
    * **Hard filters** exclude genes with highly suspicious functional annotations
    * **Soft filters** apply much stricter bitscore cutoffs to avoid ambiguous cases

*Unclassified* genes are those with annotations that don't meet thresholds for confident AVG classification, not necessarily unannotated.

### 3. What do the *viral origin confidence* assignments to predicted AVGs mean?

> **TL;DR** It reflects the likelihood that a gene is virus-encoded (vs host/MGE)

AVGs often resemble host genes and can result from contamination. CheckAMG uses local genome context to assign **high**, **medium**, or **low** viral origin confidence based on:

1. Proximity to virus-like or viral hallmark genes
2. Proximity to transposases or other non-viral mobilization genes
3. Local viral gene content (determined using [V- and VL-scores](https://github.com/AnantharamanLab/V-Score-Search) [Zhou et al., 2025](https://www.biorxiv.org/content/10.1101/2024.10.24.619987v1))
4. Contig circularity

A Random Forest model, trained on real and simulated viral/non-viral data, makes these assignments. Confidence levels refer to the viral origin, not the functional annotation.

### 4. Which confidence levels should I use?

> **TL;DR** When in doubt, use high, but medium can be included if your input is virus enriched.

The precision and recall of each confidence level for predicting true viral proteins depends on the input dataset. Whether you should use high, medium, and/or low-confidence AVGs will depend on your knowledge of your input data.

* **High-confidence**
    * CheckAMG assigns confidence levels such that high-confidence predictions can be almost always trusted (false-discovery rate < 0.1 in all tested cases)
    * To maintain the integrity of high-confidence predictions even in cases where viral proteins are relately rare in the input, high-confidence predictions are very conservative
    * **We reccomend using just high-confidence viral proteins when viral proteins are relatively rare in the input data (such as mixed-community metagenomes) or when the composition of the input data is unknown**
* **Medium-confidence**
    * Using medium-confidence predictions can dramatically increase the recovery of truly viral proteins, but they may not always be best to use
    * Medium-confidence predictions maintain false-discovery rates < 0.1 in datasets with at least 50% viral proteins, but as input sequences become increasingly non-viral in their protein composition, FDR begin to surpass 0.1 (see the table, below)
    * **We reccomend using both high- and medium-confidence viral proteins when you know that roughly half of your input sequences are viral, such as outputs from most virus prediction tools or viromes**
* **Low-confidence**
    * Low-confidence predictions are not filtered at all, so we only reccomend using them when you are certain that all of your input sequences are free of non-viral sequence contamination, or for testing

Below are preliminary results for benchmarking our viral origin confidence predictions against test datasets with varying sequence composition (% of proteins):

| Dataset              | % Viral Proteins | % MGE Proteins | % Host Proteins | Confidence Level  | Precision | Recall | F1 Score | FDR   | MCC   |
| -------------------- | ---------------- | -------------- | --------------- | ----------------- | --------- | ------ | -------- | ----- | ----- |
| near all MGE         | 5.00%            | 90.00%         | 5.00%           | High Confidence   | 0.94      | 0.2    | 0.329    | 0.06  | 0.423 |
| near all MGE         | 5.00%            | 90.00%         | 5.00%           | Medium Confidence | 0.57      | 0.848  | 0.682    | 0.43  | 0.676 |
| near all MGE         | 5.00%            | 90.00%         | 5.00%           | Low Confidence    | 0.05      | 1      | 0.095    | 0.95  | 0     |
| near all host        | 5.00%            | 5.00%          | 90.00%          | High Confidence   | 0.904     | 0.221  | 0.355    | 0.096 | 0.436 |
| near all host        | 5.00%            | 5.00%          | 90.00%          | Medium Confidence | 0.269     | 0.86   | 0.41     | 0.731 | 0.438 |
| near all host        | 5.00%            | 5.00%          | 90.00%          | Low Confidence    | 0.05      | 1      | 0.095    | 0.95  | 0     |
| MGE enriched         | 12.50%           | 75.00%         | 12.50%          | High Confidence   | 0.98      | 0.209  | 0.345    | 0.02  | 0.428 |
| MGE enriched         | 12.50%           | 75.00%         | 12.50%          | Medium Confidence | 0.736     | 0.857  | 0.792    | 0.264 | 0.762 |
| MGE enriched         | 12.50%           | 75.00%         | 12.50%          | Low Confidence    | 0.125     | 1      | 0.222    | 0.875 | 0     |
| host enriched        | 12.50%           | 12.50%         | 75.00%          | High Confidence   | 0.963     | 0.215  | 0.351    | 0.037 | 0.429 |
| host enriched        | 12.50%           | 12.50%         | 75.00%          | Medium Confidence | 0.518     | 0.863  | 0.648    | 0.482 | 0.61  |
| host enriched        | 12.50%           | 12.50%         | 75.00%          | Low Confidence    | 0.125     | 1      | 0.222    | 0.875 | 0     |
| equal source         | 33.30%           | 33.30%         | 33.30%          | High Confidence   | 0.991     | 0.204  | 0.338    | 0.009 | 0.378 |
| equal source         | 33.30%           | 33.30%         | 33.30%          | Medium Confidence | 0.845     | 0.863  | 0.854    | 0.155 | 0.78  |
| equal source         | 33.30%           | 33.30%         | 33.30%          | Low Confidence    | 0.333     | 1      | 0.5      | 0.667 | 0     |
| equal viral/nonviral | 50.00%           | 25.00%         | 25.00%          | High Confidence   | 0.996     | 0.21   | 0.347    | 0.004 | 0.341 |
| equal viral/nonviral | 50.00%           | 25.00%         | 25.00%          | Medium Confidence | 0.916     | 0.864  | 0.889    | 0.084 | 0.786 |
| equal viral/nonviral | 50.00%           | 25.00%         | 25.00%          | Low Confidence    | 0.5       | 1      | 0.667    | 0.5   | 0     |
| virus enriched       | 75.00%           | 12.50%         | 12.50%          | High Confidence   | 0.998     | 0.209  | 0.346    | 0.002 | 0.248 |
| virus enriched       | 75.00%           | 12.50%         | 12.50%          | Medium Confidence | 0.971     | 0.861  | 0.913    | 0.029 | 0.72  |
| virus enriched       | 75.00%           | 12.50%         | 12.50%          | Low Confidence    | 0.75      | 1      | 0.857    | 0.25  | 0     |
| near all virus       | 90.00%           | 5.00%          | 5.00%           | High Confidence   | 1         | 0.209  | 0.346    | 0     | 0.16  |
| near all virus       | 90.00%           | 5.00%          | 5.00%           | Medium Confidence | 0.989     | 0.861  | 0.921    | 0.011 | 0.567 |
| near all virus       | 90.00%           | 5.00%          | 5.00%           | Low Confidence    | 0.9       | 1      | 0.947    | 0.1   | 0     |

### 5. Snakemake

CheckAMG modules are executes as [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines. If a run is interrupted, it can resume from the last complete step as long as intermediate files exist.

## Error reporting

To report bugs or request features, please use the [GitHub Issues](https://github.com/AnantharamanLab/CheckAMG/issues).
## Citation

*Coming soon.*

Authors:

* James C. Kosmopoulos (**kosmopoulos [at] wisc [dot] edu**)
* Cody Martin
* Karthik Anantharaman (**karthik [at] bact [dot] wisc [dot] edu**)<sup>*</sup>

<sup>*</sup>To whom correspondence should be addressed.