# CheckAMG

**Automated curation of Auxiliary Metabolic Genes (AMGs), Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) in viral genomes.**

> ⚠️ **This tool is in active development and has not yet been peer-reviewed.**

## Quick Usage

```bash
checkamg download -d /path/to/db/destination

checkamg annotate \
  -d /path/to/db/destination \
  -g examples/example_data/single_contig_viruses.fasta \
  -vg examples/example_data/multi_contig_vMAGs \
  -o CheckAMG_example_out
```

## Features

* Input: nucleotide or protein sequences
* Handles single-contig viral genomes and multi-contig vMAGs
* Functional annotation + viral genome context-based curation
* Outputs curated lists and amino-acid sequences of AMGs, AReGs, and APGs

## Command-line Modules

```bash
checkamg -h
```

* `download`: Get required databases
* `annotate`: Predict and curate AVGs
* `de-novo`, `aggregate`, `end-to-end`: Coming soon

## Example Output

* FASTA files of predicted AVGs (by confidence and function class)
* Tabular summary of predictions (`final_results.tsv`, `gene_annotations.tsv`)

## License

GPL-3.0-or-later

**Example data and full documentation:**
[https://github.com/AnantharamanLab/CheckAMG](https://github.com/AnantharamanLab/CheckAMG)
