# CheckAMG

**Automated discovery and curation of Auxiliary Metabolic Genes (AMGs), Auxiliary Regulatory Genes (AReGs), and Auxiliary Physiology Genes (APGs) encoded by viral genomes**

> This tool is in active development and has not yet been peer-reviewed.

CheckAMG identifies and curates high-confidence auxiliary viral genes (AVGs) using two complementary approaches: an annotation-based path (curated functional annotations plus viral genome context) and a annotation-independent path (the Protein Set Transformer genome language model).

## Quick usage

```bash
checkamg download -d /path/to/db/destination

checkamg end-to-end \
  -d /path/to/db/destination \
  -i examples/example_data/single_contig_viruses.fasta \
  -I examples/example_data/multi_contig_vMAGs \
  -o CheckAMG_example_out
```

## Features

* Input: nucleotide or protein sequences
* Handles single-contig viral genomes and multi-contig vMAGs/bins
* Reference-based annotation and genome-context curation (`annotate`)
* Reference-independent prediction with a protein genome language model (`de-novo`)
* Combined reporting (`aggregate`) and one-command runs (`end-to-end`)
* Outputs curated lists and amino-acid sequences of AMGs, AReGs, and APGs

## Modules

```bash
checkamg -h
```

* `download`: get the required databases
* `annotate`: reference-based AVG prediction and curation
* `de-novo`: reference-independent AVG prediction with PST
* `aggregate`: merge `annotate` and `de-novo` results
* `end-to-end`: run annotate, de-novo, and aggregate in tandem
* `train`: finetune a custom CheckAMG-PST model

## License

GPL-3.0-or-later

**Installation, example data, and full documentation:**
[https://github.com/AnantharamanLab/CheckAMG](https://github.com/AnantharamanLab/CheckAMG) and the [CheckAMG Wiki](https://github.com/AnantharamanLab/CheckAMG/wiki).
