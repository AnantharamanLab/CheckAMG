CheckAMG end-to-end results
===========================

CheckAMG version: 1.2

Start with 03_aggregate/avg_predictions.tsv. It lists every auxiliary viral gene (AVG) supported by
annotate, de-novo, or both, and 03_aggregate/README.txt explains its columns.

01_annotate/    Output of CheckAMG annotate: annotation-based AVG predictions at every viral-origin
                confidence level (see 01_annotate/results/README.txt).
02_denovo/      Output of CheckAMG de-novo: protein language model predictions for every protein
                (predictions.tsv).
03_aggregate/   Combined results of annotate and de-novo. Start here.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
