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
### Step 1 (recommended): create a conda environment

```bash
conda create -n CheckAMG python=3.11
conda activate CheckAMG
```

### Step 2: Install from PyPI

```bash
pip install checkamg
```

### Or from bioconda

*Coming soon.*

### Step 3: Download the databases required by CheckAMG

About 40 GB of free disk space will be required to download the databases. This can be reduced to about ~21 GB after downloading finishes if the human-readable HMM files are removed by providing the `--rm_hmm` argument.

```
checkamg download -d /path/to/db/destination --rm_hmm
```

## Quick start

Example data to test your installation of CheckAMG are provided in the [`examples/example_data`](https://github.com/AnantharamanLab/CheckAMG/tree/main/examples/example_data) folder of this repository.

```
checkamg download -d /path/to/db/destination

checkamg annotate \
  -d /path/to/db/destination \
  -g examples/example_data/single_contig_viruses.fasta \
  -vg examples/example_data/multi_contig_vMAGs \
  -o CheckAMG_example_out
```

## Usage

CheckAMG has multiple modules. The main modules that will be used for AVG predictions are `annotate`, `de-novo`, and `end-to-end`. Currently, only the `annotate` module has been implemented, and its associated `download` module to download its required databases.

Run `checkamg -h` for full options and module descriptions:

```
usage: checkamg [-h] [-v] {download,annotate,de-novo,aggregate,end-to-end} ...

CheckAMG: automated identification and curation of Auxiliary Metabolic Genes (AMGs),
Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) in viral
genomes.

positional arguments:
  {download,annotate,de-novo,aggregate,end-to-end}
                        CheckAMG modules
    download            Download the databases required by CheckAMG.
    annotate            Predict and curate auxiliary genes in viral genomes based on
                        functional annotations and genomic context.
    de-novo             (Not yet implemented) Predict auxiliary genes with an annotation-
                        independent method using Protein Set Transformer.
    aggregate           (Not yet implemented) Aggregate the results of the CheckAMG
                        annotate and de-novo modules to produce a final report of
                        auxiliary gene predictions.
    end-to-end          (Not yet implemented) Executes CheckAMG annotate, de-novo, and
                        aggregate in tandem.

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
* Providing single-contig genomes or vMAGs only affects the labelling and organization of results, and does not affect AVG predictions
* Protein headers must be in [prodigal format](https://github.com/hyattpd/prodigal/wiki/understanding-the-prodigal-output#protein-translations) (e.g. `>Contig1_1 # 144 # 635 # 1` or `>Contig1_2 # 1535 # 635 # -1`)

**Full usage:**

```
usage: checkamg annotate [-h] -d DB_DIR -o OUTPUT [-g GENOMES] [-vg VMAGS]
                         [-p PROTEINS] [-vp VMAG_PROTEINS] [--input_type INPUT_TYPE]
                         [-l MIN_LEN] [-f MIN_ORF] [-n MIN_ANNOT] [-c COV_FRACTION]
                         [-e EVALUE] [-b BIT_SCORE] [-bh BITSCORE_FRACTION_HEURISTIC]
                         [-Z WINDOW_SIZE] [-F MAX_FLANK] [-V MIN_FLANK_VSCORE]
                         [-H | --use_hallmark | --no-use_hallmark] [-t THREADS]
                         [-m MEM] [--debug | --no-debug]

Predict and curate auxiliary genes in viral genomes based on functional annotations
and genomic context.

options:
  -h, --help            show this help message and exit
  --input_type INPUT_TYPE
                        Specifies whether the input files are nucleotide genomes
                        (nucl) or translated amino-acid genomes (prot). Providing
                        proteins as input will skip the pyrodigal-gv step, but it will
                        be unable to tell whether viral genomes are circular,
                        potentially losing additional evidence for verifying the viral
                        origin of putative auxiliary genes. (default: nucl).
  -l MIN_LEN, --min_len MIN_LEN
                        Minimum length in base pairs for input sequences (default:
                        5000).
  -f MIN_ORF, --min_orf MIN_ORF
                        Minimum number of open reading frames (proteins) inferred by
                        pyrodigal-gv for input sequences (default: 4).
  -n MIN_ANNOT, --min_annot MIN_ANNOT
                        Minimum percentage (0.0-1.0) of genes in a genome/contig
                        required to have been assigned a functional annotation using
                        the CheckAMG database to be considered for contextual
                        analysis. (default: 0.2).
  -c COV_FRACTION, --cov_fraction COV_FRACTION
                        Minimum covered fraction (of the user viral protein) for HMM
                        searches (default: 0.5).
  -e EVALUE, --evalue EVALUE
                        Maximum fallback E-value for HMM searches when database-
                        provided cutoffs are not available (default: 1e-05).
  -b BIT_SCORE, --bit_score BIT_SCORE
                        Minimum fallback bit score for HMM searches when database-
                        provided cutoffs are not available (default: 50).
  -bh BITSCORE_FRACTION_HEURISTIC, --bitscore_fraction_heuristic BITSCORE_FRACTION_HEURISTIC
                        Retain HMM hits scoring at least this fraction of the
                        database-provided threshold under heuristic filtering
                        (default: 0.5).
  -Z WINDOW_SIZE, --window_size WINDOW_SIZE
                        Size in base pairs of the window used to calculate the average
                        VL-score of genes on a contig (default: 25000).
  -F MAX_FLANK, --max_flank MAX_FLANK
                        Maximum length in base pairs to check on the left/right flanks
                        of potentially auxiliary genes when checking for virus-like
                        genes and non-virus-like genes (default: 5000).
  -V MIN_FLANK_VSCORE, --min_flank_Vscore MIN_FLANK_VSCORE
                        Minimum V-score of genes in flanking regions required to
                        verify a potential auxiliary gene as viral and not host
                        sequence contamination (0.0-10.0) (default: 10.0).
  -H, --use_hallmark, --no-use_hallmark
                        Use viral hallmark gene annotations instead of V-scores when
                        checking flanking regions of potential auxiliary genes for
                        viral verification (default: False).
  -t THREADS, --threads THREADS
                        Number of threads to use for pyrodigal-gv and pyhmmer
                        (default: 10).
  -m MEM, --mem MEM     Maximum amount of memory allowed to be allocated in GB
                        (default: 80% of available).
  --debug, --no-debug   Log CheckAMG with debug-level detail (default: False).

required arguments:
  -d DB_DIR, --db_dir DB_DIR
                        Path to CheckAMG database files (Required). (default: None)
  -o OUTPUT, --output OUTPUT
                        Output directory for all generated files and folders
                        (Required). (default: None)
  -g GENOMES, --genomes GENOMES
                        Input viral genome(s) in nucleotide fasta format (.fna or
                        .fasta). Expectation is that individual virus genomes are
                        single contigs. (default: None)
  -vg VMAGS, --vmags VMAGS
                        Path to folder containing vMAGs (multiple contigs) rather than
                        single-contig viral genomes. Expectation is that the folder
                        contains one .fna or .fasta file per virus genome and that
                        each genome contains multiple contigs. (default: None)
  -p PROTEINS, --proteins PROTEINS
                        Input viral genome(s) in amino-acid fasta format (.faa or
                        .fasta). Required if --input_type is prot. Expectations are
                        that the amino-acid sequence headers are in Prodigal format
                        (>[CONTIG NAME]_[CDS NUMBER] # START # END # FRAME # ...) and
                        that each contig encoding proteins represents a single virus
                        genome. (default: None)
  -vp VMAG_PROTEINS, --vmag_proteins VMAG_PROTEINS
                        Path to folder containing vMAGs (multiple contigs) in amino-
                        acid fasta format (.faa or .fasta) rather than single-contig
                        viral genomes. Expectation is that the folder contains one
                        .faa or .fasta file per virus genome and that each genome file
                        contains amino-acid sequences encoded on multiple contigs.
                        Required if --input_type is 'prot'. (default: None)
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

Examples of these output files are provided in the [`examples/example_outputs`](https://github.com/AnantharamanLab/CheckAMG/tree/main/examples/example_outputs) folder of this repository.

### CheckAMG de-novo
*Coming soon.*

### CheckAMG end-to-end
*Coming soon.*

## Important Notes / FAQs
### 1. What is an *AVG*?
An AVG is an **A**uxiliary **V**iral **G**ene, a virus-encoded gene that is non-essential for viral replication but augments host metabolism (AMGs), physiology (APGs), or regulation (AReGs). In the past, all AVGs have been referred to as AMGs, but recently the term AVG has been adopted to include broader host-modulating functions, not just metabolism.

Examples:
* A virus-encoded *psbA* or *soxY* would be an AMG because they encode proteins with functions in host photosynthesis and sulfide oxidation
* A virus-encoded VasG type VI secretion system protein or HicA toxin would be an APG because they are involved in host physiology
* A LuxR transcriptional regulator or an AsiA anti-sigma factor protein would be an AReG because they are likely involved in the regulation of host gene expression

Despite the name "CheckAMG", this tool also predicts APGs and AReGs using the same pipeline, differing only by functional annotation criteria.

### 2. How does CheckAMG classify and curate its predictions?

CheckAMG applies a two-stage filtering process:

1. Use a list of curated profile HMMs that represent [metabolic](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AMGs.tsv), [physiological](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/APGs.tsv), and [regulatory](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AReGs.tsv) genes to come up with initial AVG candidates
2. Use a second list of curated keywords/substrings that will be used to filter 'false' [AMGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_amgs.csv), [APGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_apgs.csv), and [AReGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_aregs.csv)
    * **Hard filters** exclude genes with highly suspicious functional annotations
    * **Soft filters** apply much stricter bitscore and coverage cutoffs to avoid ambiguous cases (see [*additional HMMsearch filtering for curating auxiliary gene functions*](#additional-hmmsearch-filtering-for-curating-auxiliary-gene-functions))

*Unclassified* genes are those with annotations that don't meet thresholds for confident AVG classification, not necessarily unannotated.

### 3. What do the *viral origin confidence* assignments to predicted AVGs mean?

> **TL;DR** It reflects the likelihood that a gene is virus-encoded (vs host/MGE)

AVGs often resemble host genes and can result from contamination. CheckAMG uses local genome context to assign **high**, **medium**, or **low** viral origin confidence based on:

1. Proximity to virus-like or viral hallmark genes
2. Proximity to transposases or other [non-viral mobilization genes](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/mobile_genes.csv)
3. Local viral gene content, determined using [V- and VL-scores](https://github.com/AnantharamanLab/V-Score-Search) ([Zhou et al., 2025](https://www.biorxiv.org/content/10.1101/2024.10.24.619987v1))
4. Contig circularity

A LightGBM model, trained on real and simulated viral/non-viral data, makes these assignments. Confidence levels refer to the viral origin, not the functional annotation.

### 4. Which confidence levels should I use?

> **TL;DR** When in doubt, use high, but medium should be included if your input is virus enriched.

The precision and recall of each confidence level for predicting true viral proteins depends on the input dataset. Whether you should use high, medium, and/or low-confidence AVGs will depend on your knowledge of your input data.

* **High-confidence**
    * CheckAMG assigns confidence levels such that high-confidence predictions can be almost always trusted (false-discovery rate < 0.1 in all tested cases)
    * To maintain the integrity of high-confidence predictions even in cases where viral proteins are relatively rare in the input, high-confidence predictions are conservative
    * **We recommend using just high-confidence AVGs when viral proteins are relatively rare in the input data (such as mixed-community metagenomes) or when the composition of the input data is unknown**
* **Medium-confidence**
    * Using medium-confidence predictions can significantly increase the recovery of truly viral proteins, but they may not always be best to use
    * Medium-confidence predictions maintain false-discovery rates < 0.1 in datasets with at least 50% viral proteins, but as input sequences become increasingly non-viral in their protein composition, FDRs begin to surpass 0.1 (see the figure and table, below)
    * **We recommend using both high- and medium-confidence AVGs when you know that roughly half of your input sequences are viral, such as outputs from most virus prediction tools or viromes**
* **Low-confidence**
    * Low-confidence predictions are not filtered at all, so we only recommend using them when you are certain that all of your input sequences are free of non-viral sequence contamination, or for testing

Below are preliminary results for benchmarking our viral origin confidence predictions against test datasets with varying sequence composition (% of proteins, see the table below for composition):

<img src="precision_recall_plot.png" alt="Precision-Recall Plot" style="max-width: 100%; height: auto;">

<div style="overflow-x:auto;">
<table style="width:100%; max-width:100%; table-layout:auto;">
  <thead>
    <tr style="background-color:#111827; color:#fff;">
      <th align="center">Dataset</th>
      <th align="center">% Viral Proteins</th>
      <th align="center">% MGE Proteins</th>
      <th align="center">% Host Proteins</th>
    </tr>
  </thead>
  <tbody>
    <tr style="background-color:#212C40; color:#fff;">
      <td align="center">Near all virus</td>
      <td align="center">90</td>
      <td align="center">5</td>
      <td align="center">5</td>
    </tr>
    <tr style="background-color:#3c643c; color:#fff;">
      <td align="center">Virus enriched</td>
      <td align="center">75</td>
      <td align="center">12.5</td>
      <td align="center">12.5</td>
    </tr>
    <tr style="background-color:#f3f4f6; color:#222;">
      <td align="center">Equal viral/nonviral</td>
      <td align="center">50</td>
      <td align="center">25</td>
      <td align="center">25</td>
    </tr>
    <tr style="background-color:#ffe4b5; color:#222;">
      <td align="center">Half viral/host</td>
      <td align="center">50</td>
      <td align="center">0</td>
      <td align="center">50</td>
    </tr>
    <tr style="background-color:#e3e1f7; color:#222;">
      <td align="center">Equal viral/MGE/host</td>
      <td align="center">33.3</td>
      <td align="center">33.3</td>
      <td align="center">33.3</td>
    </tr>
    <tr style="background-color:#15513a; color:#fff;">
      <td align="center">MGE enriched</td>
      <td align="center">12.5</td>
      <td align="center">75</td>
      <td align="center">12.5</td>
    </tr>
    <tr style="background-color:#fcd34d; color:#222;">
      <td align="center">Host enriched</td>
      <td align="center">12.5</td>
      <td align="center">12.5</td>
      <td align="center">75</td>
    </tr>
    <tr style="background-color:#212C40; color:#fff;">
      <td align="center">Near all MGE</td>
      <td align="center">5</td>
      <td align="center">90</td>
      <td align="center">5</td>
    </tr>
    <tr style="background-color:#f3f4f6; color:#222;">
      <td align="center">Near all host</td>
      <td align="center">5</td>
      <td align="center">5</td>
      <td align="center">90</td>
    </tr>
  </tbody>
</table>
</div>

### 5. How does CheckAMG assign functions to proteins?

If you're curious about the internal mechanics of how CheckAMG annotates proteins for their function, here's a breakdown of the behavior. These settings are designed to balance sensitivity (not missing true hits) and specificity (excluding weak/ambiguous matches), with additional database-specific optimizations for functional reliability.

1. **Homology Searching Method**

   * CheckAMG uses `pyhmmer` for fast and reproducible HMM searches of user proteins against profile HMMs

2. **Profile HMM Databases**

   * CheckAMG relies on the following databases:
     * [KEGG Orthology (KO)](https://www.genome.jp/kegg/ko.html) ([Kanehisa et al., 2016](https://doi.org/10.1093/nar/gkv1070))
     * [Functional Ontology Assignments for Metagenomes (FOAM) database](https://osf.io/5ba2v/?view_only=) ([Prestat et al., 2014](https://doi.org/10.1093/nar/gku702))
     * [Pfam-A](http://pfam.xfam.org/) ([Mistry et al., 2021](https://doi.org/10.1093/nar/gkaa913))
     * [Prokaryotic Virus Remote Homologous Groups database (PHROGs)](https://phrogs.lmge.uca.fr/) ([Terzian et al., 2021](https://doi.org/10.1093/nargab/lqab067))
     * [dbCAN CAZyme domain HMM database](https://bcb.unl.edu/dbCAN2/) ([Zheng et al., 2023](https://doi.org/10.1093/nar/gkad328))
     * [The METABOLIC HMM database](https://github.com/AnantharamanLab/METABOLIC/tree/master) ([Zhou et al., 2022](https://doi.org/10.1186/s40168-021-01213-8))
   * These databases can be downloaded and processed using the `checkamg download` module

3. **E-value Threshold**

   * An *initial*, permissive E-value cutoff of `0.1` is applied during `hmmsearch` to minimize missed hits due to chunking or memory differences when parallelizing, which can affect search reproducibility

4. **Database-Specific Thresholds**

   * CheckAMG applies specialized rules depending on the HMM source:
     * **Pfam:** Applies sequence-level *gathering threshold* (GA); hits below GA are excluded
     * **FOAM & KEGG:** Use database-defined bit score thresholds, but apply a relaxed fallback heuristic (see below)
     * **METABOLIC:** Uses GA cutoffs derived from its underlying Pfam/TIGRFAM sources, where available

5. **Fallback Heuristic (KEGG & FOAM)**

   * KEGG (and consequently, FOAM) thresholds can sometimes be overly strict, especially for environmental viruses, filtering out hits that are biologically valid
   * To recover these valid hits, CheckAMG applies a relaxed fallback heuristic inspired by the Anvi'o [`anvi-run-kegg-kofams`](https://anvio.org/help/7.1/programs/anvi-run-kegg-kofams/#how-does-it-work) strategy:
     * If a hit falls below the database-provided trusted threshold (e.g., KEGG TC), it is still retained **if the bit score is at least 50% of the threshold value** *and* **E-value is below 1e-5**
   * This heuristic improves annotation recovery without compromising too much on precision ([Kananen et al., 2025](https://doi.org/10.1093/bioadv/vbaf039))

6. **Fallback Filtering for Other Databases**

   * If the HMM source doesn't have defined cutoffs, such as dbCAN, PHROGs, and some profiles in the METABOLIC database, CheckAMG enforces:
     * A minimum sequence coverage of the user's protein of `0.5` (this is to ensure functional inferences aren't drawn solely from small, individual domains)
     * A minimum bit score of `50`
     * Both are configurable by the user if desired

7. **Result Consolidation and Best-Hit Filtering**

   * Each input protein is searched against **each HMM source database** (KEGG, FOAM, Pfam, PHROG, dbCAN, and METABOLIC)
   * All domain hits are first filtered using the criteria above
   * Then, for **each database**, only the **single best hit per protein** is retained:
     * Preference is given to the hit with the **lowest E-value**
     * If E-values are equal, the hit with the **higher bit score** is selected

These defaults provide a balance between accuracy and recall, and are based on benchmarking and community best practices. Users may modify thresholds using the `--bit_score`, `--bitscore_fraction_heuristic`, and `--evalue` arguments.

#### Additional HMMsearch filtering for curating auxiliary gene functions
As mentioned in the section [*How does CheckAMG classify and curate its predictions?*](#how-does-checkamg-classify-and-curate-its-predictions), CheckAMG applies *hard* and *soft* filters to functional annotations to reduce the chances of incorrectly assigning a metabolic/physiological/regulatory function to often mis-annotated viral genes. Annotations with hard filter keywords are excluded entirely from auxiliary gene predictions, but those containing soft filter keywords are retained if they meet stricter HMMsearch thresholds. Instead of hard-coding these thresholds, CheckAMG has the `--scaling_factor` argument:

* The value provided by `--scaling_factor` will be used to multiply the minimum bit score and minimum covered fraction provided by the `--bit_score` and `--cov_fraction` arguments to come up with the heuristic, stricter thresholds

* The default `--scaling_factor` is `1.5`, users can increase this value but we do not recommend decreasing it
  * If the default `--bit_score` and `--cov_fraction` values are also used (`50` and `0.5`), this means that a suspicious annotation containing a soft filter keyword must have had a bit score of at least `75` *and* and a sequence coverage of at least `0.75` from the HMMsearch to make it into the final auxiliary gene predictions
  * If non-default `--bit_score` and `--cov_fraction` values are used, `--scaling_factor` will apply to the user-provided values

* If database-provided, trusted bit score cutoffs are available for the matching HMMs (KEGG and FOAM only), those are used instead of the calculated heuristic threshold for minimum bit score, but the scaling factor is still applied to the minimum coverage

* This heuristic is NOT used for assigning overall gene functions as detailed above, only during the annotation curation step of CheckAMG
  * This means that you may end up with, for example, many *glycoside hydrolase* functional annotations in the `gene_annotations.tsv` output, but much fewer glycoside hydrolases classified as "metabolic" in the `final_results.tsv`, with the rest being marked as "unclassified"
  * This is because the genes with glycoside hydrolase annotations that were classified as AMGs met the heuristic thresholds defiend by `--scaling_factor`, `--bit_score`, and `--cov_fraction`, and can be considered as functioning in host carbohydrate metabolism
  * On the other hand, genes with glycoside hydrolase annotations that were NOT classified as AMGs and ended up as "unclassified" did not met the stricter thresholds, and they are more likely to be involved in functions other than host carbohydrate metabolism (like breaking down the surface of the host's cell wall, for example)

### 6. Snakemake

CheckAMG modules are executed as [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines. If a run is interrupted, it can resume from the last complete step as long as intermediate files exist.

## Error reporting

To report bugs or request features, please use the [GitHub Issues page](https://github.com/AnantharamanLab/CheckAMG/issues).

## Citation

*Coming soon.*

Authors:

* James C. Kosmopoulos (**kosmopoulos [at] wisc [dot] edu**)
* Cody Martin
* Karthik Anantharaman (**karthik [at] bact [dot] wisc [dot] edu**)
