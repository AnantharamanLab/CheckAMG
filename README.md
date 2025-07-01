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

```bash
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

```bash
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

The CheckAMG output folder will have the following structure:

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

**Example output tables:**

`final_results.tsv`:

| Protein                                        | Contig                                       | Genome                                       | Protein Classification | Protein Viral Origin Confidence | KEGG KO | KEGG KO Name                                                           | FOAM ID       | FOAM Annotation                                                                          | Pfam Accession                  | Pfam Name                                        | CAZy Family                               | CAZy Activities                                     | METABOLIC db ID                           | METABOLIC Annotation                                        | PHROG Number                                                           | PHROG Annotation                      | Best Scoring HMM           | Best Scoring HMM Annotation | Best Scoring HMM Origin |
| ---------------------------------------------- | -------------------------------------------- | -------------------------------------------- | ---------------------- | ------------------------------- | ------- | ---------------------------------------------------------------------- | ------------- | ---------------------------------------------------------------------------------------- | ------------------------------- | ------------------------------------------------ | ----------------------------------------- | --------------------------------------------------- | ----------------------------------------- | ----------------------------------------------------------- | ---------------------------------------------------------------------- | ------------------------------------- | -------------------------- | --------------------------- | ----------------------- |
| vRhyme_bin_8__LASCr2D2E2F_000000298731_22      | vRhyme_bin_8__LASCr2D2E2F_000000298731       | LASCr2D2E2F_vRhyme_bin_8                     | metabolic              | high                            |         |                                                                        | HMMsoil68357  | alaA; alanine-synthesizing transaminase EC:2.6.1.66 2.6.1.2                              |                                 |                                                  |                                           |                                                     | HMMsoil68357                              | alaA; alanine-synthesizing transaminase EC:2.6.1.66 2.6.1.2 | FOAM                                                                   |
| vRhyme_bin_10__MGr2D2E2F_000000095924_3        | vRhyme_bin_10__MGr2D2E2F_000000095924        | MGr2D2E2F_vRhyme_bin_10                      | metabolic              | medium                          |         |                                                                        | HMMsoil57593  | cysH; phosphoadenosine phosphosulfate reductase EC:1.8.4.8 1.8.4.10                      |                                 | phrog_424                                        | phosphoadenosine phosphosulfate reductase | phrog_424                                           | phosphoadenosine phosphosulfate reductase | PHROG                                                       |
| vRhyme_bin_11__CRr2D2E2F_000000051715_10       | vRhyme_bin_11__CRr2D2E2F_000000051715        | CRr2D2E2F_vRhyme_bin_11                      | metabolic              | medium                          | K02078  | acyl carrier protein                                                   | HMMsoil119881 | NDUFAB1; NADH dehydrogenase (ubiquinone) 1 alpha/beta subcomplex 1, acyl-carrier protein | PF00550.29                      | Phosphopantetheine attachment site               |                                           |                                                     |                                           |                                                             | K02078                                                                 | acyl carrier protein                  | KEGG                       |
| vRhyme_bin_11__LABRr1A1B1C_000000203956_1      | vRhyme_bin_11__LABRr1A1B1C_000000203956      | LABRr1A1B1C_vRhyme_bin_11                    | metabolic              | medium                          |         |                                                                        | HMMsoil44705  | rpe, RPE; ribulose-phosphate 3-epimerase EC:5.1.3.1                                      | PF00834.23                      | Ribulose-phosphate 3 epimerase family            |                                           |                                                     |                                           |                                                             | PF00834.23                                                             | Ribulose-phosphate 3 epimerase family | Pfam                       |
| vRhyme_unbinned_109__SEr3G3H3I_000000174344_18 | vRhyme_unbinned_109__SEr3G3H3I_000000174344  | vRhyme_unbinned_109__SEr3G3H3I_000000174344  | metabolic              | low                             |         |                                                                        | HMMsoil31643  | cysH; phosphoadenosine phosphosulfate reductase EC:1.8.4.8 1.8.4.10                      | PF01507.23                      | Phosphoadenosine phosphosulfate reductase family |                                           |                                                     |                                           | HMMsoil31643                                                | cysH; phosphoadenosine phosphosulfate reductase EC:1.8.4.8 1.8.4.10    | FOAM                                  |
| vRhyme_bin_11__BAr1A1B1C_000000245585_10       | vRhyme_bin_11__BAr1A1B1C_000000245585        | BAr1A1B1C_vRhyme_bin_11                      | physiological          | medium                          | K03704  | cold shock protein                                                     |               | PF00313.26                                                                               | 'Cold-shock' DNA-binding domain |                                                  |                                           | phrog_729                                           | cold shock protein                        | K03704                                                      | cold shock protein                                                     | KEGG                                  |
| vRhyme_bin_12__LASCr2D2E2F_000000190496_2      | vRhyme_bin_12__LASCr2D2E2F_000000190496      | LASCr2D2E2F_vRhyme_bin_12                    | physiological          | medium                          | K06348  | sporulation inhibitor KapD                                             | HMMsoil65112  | dnaQ; DNA polymerase III subunit epsilon EC:2.7.7.7                                      | PF00929.28                      | Exonuclease                                      |                                           |                                                     |                                           | phrog_192                                                   | DNA polymerase exonuclease subunit                                     | K06348                                | sporulation inhibitor KapD | KEGG                        |
| vRhyme_unbinned_10__LASCr2D2E2F_000000024005_2 | vRhyme_unbinned_10__LASCr2D2E2F_000000024005 | vRhyme_unbinned_10__LASCr2D2E2F_000000024005 | physiological          | low                             | K02456  | general secretion pathway protein G                                    | PF07596.15    | Protein of unknown function (DUF1559)                                                    |                                 |                                                  |                                           |                                                     | PF07596.15                                | Protein of unknown function (DUF1559)                       | Pfam                                                                   |
| vRhyme_unbinned_153__SEr3G3H3I_000000269915_1  | vRhyme_unbinned_153__SEr3G3H3I_000000269915  | vRhyme_unbinned_153__SEr3G3H3I_000000269915  | physiological          | low                             | K27087  | ESX secretion system protein EccD                                      | PF08817.14    | WXG100 protein secretion system (Wss), protein YukD                                      |                                 |                                                  | PF08817.14                                | WXG100 protein secretion system (Wss), protein YukD | Pfam                                      |
| vRhyme_unbinned_17__BAr2D2E2F_000000115566_6   | vRhyme_unbinned_17__BAr2D2E2F_000000115566   | vRhyme_unbinned_17__BAr2D2E2F_000000115566   | physiological          | low                             | K10859  | DNA oxidative demethylase [EC:1.14.11.33]                              | PF13532.10    | 2OG-Fe(II) oxygenase superfamily                                                         |                                 |                                                  | phrog_21963                               | alkylated DNA repair                                | K10859                                    | DNA oxidative demethylase [EC:1.14.11.33]                   | KEGG                                                                   |
| vRhyme_bin_11__MGr1A1B1C_000000108133_4        | vRhyme_bin_11__MGr1A1B1C_000000108133        | MGr1A1B1C_vRhyme_bin_11                      | regulatory             | high                            | K03497  | ParB family transcriptional regulator, chromosome partitioning protein | PF02195.22    | ParB/Sulfiredoxin domain                                                                 |                                 |                                                  |                                           | phrog_1033                                          | ParB-like partition protein               | phrog_1033                                                  | ParB-like partition protein                                            | PHROG                                 |
| vRhyme_bin_10__LAWAr2D2E2F_000000097578_27     | vRhyme_bin_10__LAWAr2D2E2F_000000097578      | LAWAr2D2E2F_vRhyme_bin_10                    | regulatory             | medium                          | K04078  | chaperonin GroES                                                       |               |                                                                                          |                                 |                                                  |                                           |                                                     |                                           |                                                             |                                                                        | K04078                                | chaperonin GroES           | KEGG                        |
| vRhyme_bin_10__LAWAr2D2E2F_000000171335_10     | vRhyme_bin_10__LAWAr2D2E2F_000000171335      | LAWAr2D2E2F_vRhyme_bin_10                    | regulatory             | medium                          | K05788  | integration host factor subunit beta                                   | PF00216.25    | Bacterial DNA-binding protein                                                            |                                 |                                                  | phrog_379                                 | DNA binding protein                                 | K05788                                    | integration host factor subunit beta                        | KEGG                                                                   |
| vRhyme_bin_10__MHr2D2E2F_000000121054_7        | vRhyme_bin_10__MHr2D2E2F_000000121054        | MHr2D2E2F_vRhyme_bin_10                      | regulatory             | medium                          | K03497  | ParB family transcriptional regulator, chromosome partitioning protein | PF02195.22    | ParB/Sulfiredoxin domain                                                                 |                                 |                                                  |                                           |                                                     |                                           | K03497                                                      | ParB family transcriptional regulator, chromosome partitioning protein | KEGG                                  |
| vRhyme_unbinned_82__LABRr1A1B1C_000000211601_4 | vRhyme_unbinned_82__LABRr1A1B1C_000000211601 | vRhyme_unbinned_82__LABRr1A1B1C_000000211601 | regulatory             | low                             | K07733  | prophage regulatory protein                                            | PF12728.11    | Helix-turn-helix domain                                                                  |                                 |                                                  |                                           |                                                     |                                           | PF12728.11                                                  | Helix-turn-helix domain                                                | Pfam                                  |
| vRhyme_bin_10__BAr1A1B1C_000000114367_27       | vRhyme_bin_10__BAr1A1B1C_000000114367        | BAr1A1B1C_vRhyme_bin_10                      | unclassified           | high                            |         |                                                                        | HMMsoil3470   | UGDH, ugd; UDPglucose 6-dehydrogenase EC:1.1.1.22                                        | PF01381.26                      | Helix-turn-helix                                 |                                           |                                                     |                                           |                                                             |                                                                        | PF01381.26                            | Helix-turn-helix           | Pfam                        |
| vRhyme_bin_10__BAr1A1B1C_000000114367_32       | vRhyme_bin_10__BAr1A1B1C_000000114367        | BAr1A1B1C_vRhyme_bin_10                      | unclassified           | high                            |         |                                                                        |               |                                                                                          | PF13479.10                      | AAA domain                                       |                                           |                                                     |                                           |                                                             | phrog_1366                                                             |                                       | phrog_1366                 |                             | PHROG                   |
| vRhyme_bin_10__BAr1A1B1C_000000114367_38       | vRhyme_bin_10__BAr1A1B1C_000000114367        | BAr1A1B1C_vRhyme_bin_10                      | unclassified           | high                            | K02314  | replicative DNA helicase [EC:5.6.2.3]                                  | PF03796.19    | DnaB-like helicase C terminal domain                                                     |                                 |                                                  | phrog_19                                  | DnaB-like replicative helicase                      | K02314                                    | replicative DNA helicase [EC:5.6.2.3]                       | KEGG                                                                   |
| vRhyme_bin_10__CRr2D2E2F_000000254623_113      | vRhyme_bin_10__CRr2D2E2F_000000254623        | CRr2D2E2F_vRhyme_bin_10                      | unclassified           | medium                          | K01159  | crossover junction endodeoxyribonuclease RuvC [EC:3.1.21.10]           |               |                                                                                          |                                 |                                                  | phrog_159                                 | RuvC-like Holliday junction resolvase               | phrog_159                                 | RuvC-like Holliday junction resolvase                       | PHROG                                                                  |
| vRhyme_bin_10__CRr2D2E2F_000000254623_114      | vRhyme_bin_10__CRr2D2E2F_000000254623        | CRr2D2E2F_vRhyme_bin_10                      | unclassified           | medium                          | K03553  | recombination protein RecA                                             | PF00154.25    | recA bacterial DNA recombination protein                                                 |                                 | phrog_97                                         | UvsX-like recombinase                     | K03553                                              | recombination protein RecA                | KEGG                                                        |

`gene_annotations.tsv`:

| Protein                                  | Contig                                | Genome                  | KEGG_V-score | Pfam_V-score | PHROG_V-score | KEGG_hmm_id | KEGG_Description                | KEGG_score | KEGG_coverage | FOAM_hmm_id                           | FOAM_Description | FOAM_score | FOAM_coverage | Pfam_hmm_id                                | Pfam_Description              | Pfam_score | Pfam_coverage | dbCAN_hmm_id | dbCAN_Description | dbCAN_score | dbCAN_coverage | METABOLIC_hmm_id | METABOLIC_Description | METABOLIC_score | METABOLIC_coverage | PHROG_hmm_id | PHROG_Description                                     | PHROG_score | PHROG_coverage | top_hit_hmm_id                             | top_hit_description                                   | top_hit_db | Circular_Contig | Viral_Origin_Confidence | Viral_Flanking_Genes_Upstream | Viral_Flanking_Genes_Downstream | MGE_Flanking_Genes |
| ---------------------------------------- | ------------------------------------- | ----------------------- | ------------ | ------------ | ------------- | ----------- | ------------------------------- | ---------- | ------------- | ------------------------------------- | ---------------- | ---------- | ------------- | ------------------------------------------ | ----------------------------- | ---------- | ------------- | ------------ | ----------------- | ----------- | -------------- | ---------------- | --------------------- | --------------- | ------------------ | ------------ | ----------------------------------------------------- | ----------- | -------------- | ------------------------------------------ | ----------------------------------------------------- | ---------- | --------------- | ----------------------- | ----------------------------- | ------------------------------- | ------------------ |
| vRhyme_bin_14__BAr1A1B1C_000000019024_1  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 | 10           | 10           | 10            | K03546      | DNA repair protein SbcC/Rad50   | 95.094658  | 0.956         |                                       |                  |            |               | PF13476.10                                 | AAA domain                    | 38.264549  | 0.546         |              |                   |             |                |                  |                       |                 |                    | phrog_77     | SbcC-like subunit of palindrome specific endonuclease | 114.687996  | 0.962          | phrog_77                                   | SbcC-like subunit of palindrome specific endonuclease | PHROG      | FALSE           | medium                  |                               |                                 |                    |
| vRhyme_bin_14__BAr1A1B1C_000000019024_2  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 |              |              |               |             |                                 |            |               |                                       |                  |            |               |                                            |                               |            |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             |                | KEGG                                       | FALSE                                                 | medium     |                 |                         |                               |
| vRhyme_bin_14__BAr1A1B1C_000000019024_3  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 |              |              |               |             |                                 |            |               |                                       |                  |            |               |                                            |                               |            |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             |                | KEGG                                       | FALSE                                                 | medium     |                 |                         |                               |
| vRhyme_bin_14__BAr1A1B1C_000000019024_4  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 | 10           | 10           | 10            | K04077      | chaperonin GroEL [EC:5.6.1.7]   | 234.647522 | 0.913         |                                       |                  |            |               | PF00118.28                                 | TCP-1/cpn60 chaperonin family | 80.925293  | 0.328         |              |                   |             |                |                  |                       |                 |                    | phrog_5725   | chaperonin groEL                                      | 270.102173  | 0.894          | phrog_5725                                 | chaperonin groEL                                      | PHROG      | FALSE           | medium                  |                               |                                 |                    |
| vRhyme_bin_14__BAr1A1B1C_000000019024_5  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 |              |              |               |             |                                 |            |               |                                       |                  |            |               |                                            |                               |            |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             |                | KEGG                                       | FALSE                                                 | medium     |                 |                         |                               |
| vRhyme_bin_14__BAr1A1B1C_000000019024_6  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 |              |              |               |             |                                 |            |               |                                       |                  |            |               |                                            |                               |            |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             |                | KEGG                                       | FALSE                                                 | medium     |                 |                         |                               |
| vRhyme_bin_14__BAr1A1B1C_000000019024_7  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 | 10           | 10           | 10            | K09935      | N-glycosidase YbiA [EC:3.2.2.-] | 209.760986 | 0.979         |                                       |                  |            |               | PF08719.15                                 | NADAR domain                  | 122.976585 | 0.952         |              |                   |             |                |                  |                       |                 |                    | phrog_1017   |                                                       | 125.373451  | 0.911          | K09935                                     | N-glycosidase YbiA [EC:3.2.2.-]                       | KEGG       | FALSE           | medium                  |                               |                                 |                    |
| vRhyme_bin_14__BAr1A1B1C_000000019024_8  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 | 10           | 10           |               |             |                                 |            | HMMsoil40289  | FANCM; fanconi anemia group M protein | 19.662205        | 0.235      | PF04851.19    | Type III restriction enzyme, res subunit   | 60.936756                     | 0.294      |               |              |                   |             |                |                  |                       |                 | phrog_16           | DNA helicase | 204.102753                                            | 0.82        | phrog_16       | DNA helicase                               | PHROG                                                 | FALSE      | medium          |                         |                               |                                 |
| vRhyme_bin_14__BAr1A1B1C_000000019024_9  | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 | 10           |              |               |             |                                 |            |               |                                       |                  |            | PF05272.15    | Virulence-associated protein E-like domain | 48.510647                     | 0.203      |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             | PF05272.15     | Virulence-associated protein E-like domain | Pfam                                                  | FALSE      | medium          |                         |                               |                                 |
| vRhyme_bin_14__BAr1A1B1C_000000019024_10 | vRhyme_bin_14__BAr1A1B1C_000000019024 | BAr1A1B1C_vRhyme_bin_14 |              |              |               |             |                                 |            |               |                                       |                  |            |               |                                            |                               |            |               |              |                   |             |                |                  |                       |                 |                    |              |                                                       |             |                | KEGG                                       | FALSE                                                 | medium     |                 |                         |                               |

`genes_genomic_context.tsv`:

| genome                  | contig                                | protein                                  | gene_number | contig_pos_start | contig_pos_end | frame | is_vMAG | METABOLIC_score | PHROG_score | KEGG_score | Pfam_score | dbCAN_score | FOAM_score | METABOLIC_hmm_id | PHROG_hmm_id | KEGG_hmm_id | Pfam_hmm_id | dbCAN_hmm_id | FOAM_hmm_id  | METABOLIC_coverage | PHROG_coverage | KEGG_coverage | Pfam_coverage | dbCAN_coverage | FOAM_coverage | Pfam_hmm_coverage | Pfam_V-score | Pfam_VL-score | KEGG_hmm_coverage | KEGG_V-score | KEGG_VL-score | PHROG_hmm_coverage | PHROG_V-score | PHROG_VL-score | gene_length_bases | prot_length_AAs | contig_avg_KEGG_V-score | contig_avg_Pfam_V-score | contig_avg_PHROG_V-score | contig_avg_KEGG_VL-score | contig_avg_Pfam_VL-score | contig_avg_PHROG_VL-score | circular_contig | window_avg_KEGG_VL-score | window_avg_Pfam_VL-score | window_avg_PHROG_VL-score | window_avg_KEGG_V-score | window_avg_Pfam_V-score | window_avg_PHROG_V-score | KEGG_verified_flank_up | KEGG_verified_flank_down | Pfam_verified_flank_up | Pfam_verified_flank_down | PHROG_verified_flank_up | PHROG_verified_flank_down | KEGG_MGE_flank | Pfam_MGE_flank | PHROG_MGE_flank | RF_viral_prob | Viral_Origin_Confidence |
| ----------------------- | ------------------------------------- | ---------------------------------------- | ----------- | ---------------- | -------------- | ----- | ------- | --------------- | ----------- | ---------- | ---------- | ----------- | ---------- | ---------------- | ------------ | ----------- | ----------- | ------------ | ------------ | ------------------ | -------------- | ------------- | ------------- | -------------- | ------------- | ----------------- | ------------ | ------------- | ----------------- | ------------ | ------------- | ------------------ | ------------- | -------------- | ----------------- | --------------- | ----------------------- | ----------------------- | ------------------------ | ------------------------ | ------------------------ | ------------------------- | --------------- | ------------------------ | ------------------------ | ------------------------- | ----------------------- | ----------------------- | ------------------------ | ---------------------- | ------------------------ | ---------------------- | ------------------------ | ----------------------- | ------------------------- | -------------- | -------------- | --------------- | ------------- | ----------------------- |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_1  | 1           | 2                | 1018           | \-1   | TRUE    |                 | 114.687996  | 95.094658  | 38.264549  |             |            |                  | phrog_77     | K03546      | PF13476.10  |              |              |                    | 0.962          | 0.956         | 0.546         |                |               | 0.546             | 10           | 3.99973935    | 0.956             | 10           | 4.33543784    | 0.962              | 10            | 3.85009459     | 1017              | 339             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | FALSE                  | TRUE                     | FALSE                  | TRUE                     | FALSE                   | TRUE                      | FALSE          | FALSE          | FALSE           | 0.78827778    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_2  | 2           | 1119             | 1385           | \-1   | TRUE    |                 |             |            |            |             |            |                  |              |             |             |              |              |                    |                |               |               |                |               |                   |              |               |                   |              |               |                    |               |                | 267               | 89              | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_3  | 3           | 1385             | 1801           | \-1   | TRUE    |                 |             |            |            |             |            |                  |              |             |             |              |              |                    |                |               |               |                |               |                   |              |               |                   |              |               |                    |               |                | 417               | 139             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_4  | 4           | 1798             | 3582           | \-1   | TRUE    |                 | 270.102173  | 234.647522 | 80.925293  |             |            |                  | phrog_5725   | K04077      | PF00118.28  |              |              |                    | 0.894          | 0.913         | 0.328         |                |               | 0.328             | 10           | 3.65753389    | 0.913             | 10           | 3.66651798    | 0.894              | 10            | 3.65195607     | 1785              | 595             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80577778    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_5  | 5           | 3605             | 3802           | \-1   | TRUE    |                 |             |            |            |             |            |                  |              |             |             |              |              |                    |                |               |               |                |               |                   |              |               |                   |              |               |                    |               |                | 198               | 66              | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_6  | 6           | 3799             | 4185           | \-1   | TRUE    |                 |             |            |            |             |            |                  |              |             |             |              |              |                    |                |               |               |                |               |                   |              |               |                   |              |               |                    |               |                | 387               | 129             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_7  | 7           | 4220             | 4657           | \-1   | TRUE    |                 | 125.373451  | 209.760986 | 122.976585 |             |            |                  | phrog_1017   | K09935      | PF08719.15  |              |              |                    | 0.911          | 0.979         | 0.952         |                |               | 0.952             | 10           | 3.40036527    | 0.979             | 10           | 3.42667389    | 0.911              | 10            | 3.33304403     | 438               | 146             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80605556    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_8  | 8           | 4667             | 6133           | \-1   | TRUE    |                 | 204.102753  |            | 60.936756  |             | 19.662205  |                  | phrog_16     |             | PF04851.19  |              | HMMsoil40289 | 0.82               |                | 0.294         |               | 0.235          | 0.294         | 10                | 4.53420385   |               |                   |              | 0.82          | 10                 | 4.16988007    | 1467           | 489               | 10              | 8.72142857              | 10                      | 3.8395377                | 3.66184335               | 3.74776275               | FALSE                     | 3.8395377       | 3.66184335               | 3.74776275               | 10                        | 8.72142857              | 10                      | TRUE                     | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | FALSE                     | FALSE          | FALSE          | 0.80966667      | medium        |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_9  | 9           | 6214             | 8361           | 1     | TRUE    |                 |             |            | 48.510647  |             |            |                  |              |             | PF05272.15  |              |              |                    |                |               | 0.203         |                |               | 0.203             | 10           | 3.88671627    |                   |              |               |                    |               |                | 2148              | 716             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |
| BAr1A1B1C_vRhyme_bin_14 | vRhyme_bin_14__BAr1A1B1C_000000019024 | vRhyme_bin_14__BAr1A1B1C_000000019024_10 | 10          | 8942             | 9988           | 1     | TRUE    |                 |             |            |            |             |            |                  |              |             |             |              |              |                    |                |               |               |                |               |                   |              |               |                   |              |               |                    |               |                | 1047              | 349             | 10                      | 8.72142857              | 10                       | 3.8395377                | 3.66184335               | 3.74776275                | FALSE           | 3.8395377                | 3.66184335               | 3.74776275                | 10                      | 8.72142857              | 10                       | TRUE                   | TRUE                     | TRUE                   | TRUE                     | TRUE                    | TRUE                      | FALSE          | FALSE          | FALSE           | 0.80911111    | medium                  |

`metabolic_genes_curated.tsv`:

| Protein                                  | Contig                                | Genome                  | KEGG_V-score | Pfam_V-score | PHROG_V-score | KEGG_hmm_id | KEGG_Description                    | KEGG_score    | KEGG_coverage                                                                 | FOAM_hmm_id                                                                 | FOAM_Description                                                                         | FOAM_score | FOAM_coverage | Pfam_hmm_id        | Pfam_Description                   | Pfam_score | Pfam_coverage | dbCAN_hmm_id | dbCAN_Description | dbCAN_score | dbCAN_coverage | METABOLIC_hmm_id | METABOLIC_Description | METABOLIC_score | METABOLIC_coverage | PHROG_hmm_id | PHROG_Description    | PHROG_score   | PHROG_coverage                                                                | top_hit_hmm_id            | top_hit_description                 | top_hit_db | Circular_Contig | Viral_Origin_Confidence | Viral_Flanking_Genes_Upstream | Viral_Flanking_Genes_Downstream | MGE_Flanking_Genes |
| ---------------------------------------- | ------------------------------------- | ----------------------- | ------------ | ------------ | ------------- | ----------- | ----------------------------------- | ------------- | ----------------------------------------------------------------------------- | --------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------- | ---------- | ------------- | ------------------ | ---------------------------------- | ---------- | ------------- | ------------ | ----------------- | ----------- | -------------- | ---------------- | --------------------- | --------------- | ------------------ | ------------ | -------------------- | ------------- | ----------------------------------------------------------------------------- | ------------------------- | ----------------------------------- | ---------- | --------------- | ----------------------- | ----------------------------- | ------------------------------- | ------------------ |
| vRhyme_bin_14__BAr1A1B1C_000000077132_21 | vRhyme_bin_14__BAr1A1B1C_000000077132 | BAr1A1B1C_vRhyme_bin_14 | 0.8          |              |               |             |                                     |               | HMMsoil123373                                                                 | AADAT, KAT2; kynurenine/2-aminoadipate aminotransferase EC:2.6.1.7 2.6.1.39 | 12.580647                                                                                | 0.431      | PF04864.17    | Allinase           | 47.909111                          | 0.584      |               |              |                   |             |                |                  |                       |                 |                    |              |                      |               | PF04864.17                                                                    | Allinase                  | Pfam                                | FALSE      | medium          | TRUE                    | TRUE                          | FALSE                           |
| vRhyme_bin_17__BAr1A1B1C_000000024206_9  | vRhyme_bin_17__BAr1A1B1C_000000024206 | BAr1A1B1C_vRhyme_bin_17 | 10           | 7.52         |               | K02078      | acyl carrier protein                | 63.143044     | 0.967                                                                         | HMMsoil119881                                                               | NDUFAB1; NADH dehydrogenase (ubiquinone) 1 alpha/beta subcomplex 1, acyl-carrier protein | 45.002964  | 0.703         | PF00550.29         | Phosphopantetheine attachment site | 37.689743  | 0.714         |              |                   |             |                |                  |                       |                 |                    |              |                      |               |                                                                               | K02078                    | acyl carrier protein                | KEGG       | FALSE           | medium                  | TRUE                          | TRUE                            | FALSE              |
| vRhyme_bin_23__BAr1A1B1C_000000238114_13 | vRhyme_bin_23__BAr1A1B1C_000000238114 | BAr1A1B1C_vRhyme_bin_23 | 2.61         | 1.5          |               | K00859      | dephospho-CoA kinase [EC:2.7.1.24]  | 113.974876    | 0.952                                                                         |                                                                             |                                                                                          |            |               | PF01121.24         | Dephospho-CoA kinase               | 68.866798  | 0.878         |              |                   |             |                |                  |                       |                 |                    |              |                      |               |                                                                               | K00859                    | dephospho-CoA kinase [EC:2.7.1.24]  | KEGG       | FALSE           | low                     | TRUE                          | TRUE                            | FALSE              |
| vRhyme_bin_23__BAr1A1B1C_000000238114_29 | vRhyme_bin_23__BAr1A1B1C_000000238114 | BAr1A1B1C_vRhyme_bin_23 | 0.84         | 0.61         | 0.43          | K00858      | NAD+ kinase [EC:2.7.1.23]           | 221.885864    | 0.976                                                                         |                                                                             |                                                                                          |            |               | PF20143.3          | ATP-NAD kinase C-terminal domain   | 111.877022 | 0.432         |              |                   |             |                |                  |                       |                 |                    | phrog_13155  | 175.80899            | 0.986         | K00858                                                                        | NAD+ kinase [EC:2.7.1.23] | KEGG                                | FALSE      | low             | TRUE                    | TRUE                          | FALSE                           |
| vRhyme_bin_3__BAr1A1B1C_000000145806_5   | vRhyme_bin_3__BAr1A1B1C_000000145806  | BAr1A1B1C_vRhyme_bin_3  |              |              |               |             |                                     | HMMsoil91     | AANAT; arylalkylamine N-acetyltransferase EC:2.3.1.87                         | 19.609152                                                                   | 0.074                                                                                    |            |               |                    |                                    |            |               |              |                   |             |                |                  |                       |                 |                    |              |                      | HMMsoil91     | AANAT; arylalkylamine N-acetyltransferase EC:2.3.1.87                         | FOAM                      | FALSE                               | medium     | TRUE            | FALSE                   | FALSE                         |
| vRhyme_bin_7__BAr1A1B1C_000000112348_1   | vRhyme_bin_7__BAr1A1B1C_000000112348  | BAr1A1B1C_vRhyme_bin_7  | 10           |              |               |             |                                     |               | HMMsoil79621                                                                  | algL; poly(beta-D-mannuronate) lyase EC:4.2.2.3                             | 37.23978                                                                                 | 0.207      | PF00754.29    | F5/8 type C domain | 42.205345                          | 0.306      |               |              |                   |             |                |                  |                       |                 |                    |              |                      |               | PF00754.29                                                                    | F5/8 type C domain        | Pfam                                | FALSE      | medium          | FALSE                   | TRUE                          | FALSE                           |
| vRhyme_bin_15__BAr2D2E2F_000000273613_66 | vRhyme_bin_15__BAr2D2E2F_000000273613 | BAr2D2E2F_vRhyme_bin_15 | 10           | 7.52         | 4             | K02078      | acyl carrier protein                | 86.53978      | 0.915                                                                         | HMMsoil119881                                                               | NDUFAB1; NADH dehydrogenase (ubiquinone) 1 alpha/beta subcomplex 1, acyl-carrier protein | 68.364075  | 0.817         | PF00550.29         | Phosphopantetheine attachment site | 48.093494  | 0.817         |              |                   |             |                |                  |                       |                 |                    | phrog_5934   | acyl carrier protein | 73.046074     | 0.927                                                                         | K02078                    | acyl carrier protein                | KEGG       | FALSE           | medium                  | FALSE                         | FALSE                           | FALSE              |
| vRhyme_bin_15__BAr2D2E2F_000000362517_5  | vRhyme_bin_15__BAr2D2E2F_000000362517 | BAr2D2E2F_vRhyme_bin_15 | 10           | 10           | 10            | K00472      | prolyl 4-hydroxylase [EC:1.14.11.2] | 227.06604     | 0.706                                                                         | HMMsoil135564                                                               | P4HA; prolyl 4-hydroxylase EC:1.14.11.2                                                  | 53.87035   | 0.251         | PF13640.10         | 2OG-Fe(II) oxygenase superfamily   | 57.089725  | 0.376         |              |                   |             |                |                  |                       |                 |                    | phrog_61     | 2OG-Fe(II) oxygenase | 79.976585     | 0.717                                                                         | K00472                    | prolyl 4-hydroxylase [EC:1.14.11.2] | KEGG       | FALSE           | medium                  | FALSE                         | TRUE                            | FALSE              |
| vRhyme_bin_17__BAr2D2E2F_000000399195_14 | vRhyme_bin_17__BAr2D2E2F_000000399195 | BAr2D2E2F_vRhyme_bin_17 | 2.61         | 1.5          |               | K00859      | dephospho-CoA kinase [EC:2.7.1.24]  | 155.239059    | 0.936                                                                         |                                                                             |                                                                                          |            |               | PF01121.24         | Dephospho-CoA kinase               | 116.545601 | 0.871         |              |                   |             |                |                  |                       |                 |                    |              |                      |               |                                                                               | K00859                    | dephospho-CoA kinase [EC:2.7.1.24]  | KEGG       | FALSE           | medium                  | TRUE                          | FALSE                           | TRUE               |
| vRhyme_bin_2__BAr2D2E2F_000000004288_12  | vRhyme_bin_2__BAr2D2E2F_000000004288  | BAr2D2E2F_vRhyme_bin_2  |              |              |               |             |                                     | HMMsoil130274 | korD, oorD; 2-oxoglutarate ferredoxin oxidoreductase subunit delta EC:1.2.7.3 | 20.312967                                                                   | 0.851                                                                                    |            |               |                    |                                    |            |               |              |                   |             |                |                  |                       |                 |                    |              |                      | HMMsoil130274 | korD, oorD; 2-oxoglutarate ferredoxin oxidoreductase subunit delta EC:1.2.7.3 | FOAM                      | FALSE                               | medium     | TRUE            | TRUE                    | FALSE                         |

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