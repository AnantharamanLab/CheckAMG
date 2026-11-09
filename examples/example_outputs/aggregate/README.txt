CheckAMG aggregate results
==========================

CheckAMG version: 1.2

This run
--------
Proteins:  3,840
AVGs supported by annotate, de-novo, or both: 83
  annotate and de-novo: 33
  annotate only: 11
  de-novo only: 39
  by final AVG confidence: very high 18, high 25, medium 29, low 11

Start here
----------
1. avg_predictions.tsv    One row per auxiliary viral gene (AVG) supported by annotate, de-novo, or
                          both. This is the main result.
2. genome_summary.tsv     One row per genome, with AVG counts by final AVG confidence, AVG type,
                          and evidence class.
3. detailed/              Supporting tables:
                            all_proteins.tsv                 One row per protein, including proteins
                                                             not supported as AVGs.
                            aggregated_results_detailed.tsv  Every probability and confidence level
                                                             from both modules, for every protein.
                            run_summary.tsv                  Inputs, support definitions, and counts.

AVG types
---------
AMG   auxiliary metabolic gene
APG   auxiliary physiology gene
AReG  auxiliary regulatory gene

Which confidence to use
-----------------------
Final AVG Confidence Level  From de-novo. Combines the probability that the protein is auxiliary
                            (AVG-like) with the probability that it is viral: very high, high,
                            medium, or low. When it is available, use it to judge AVG predictions.
Viral Origin Confidence (annotate)
                            From annotate's genome-context classifier: high, medium, or low. This
                            is supporting evidence for viral origin.

Evidence Class
--------------
Every row of avg_predictions.tsv is supported by at least one module:
  annotate and de-novo  Annotate classifies the protein as an AVG with high or medium viral-origin
                        confidence, and de-novo assigns a final AVG confidence of very high, high,
                        or medium. These have the most support.
  annotate only         Supported by annotate but not by de-novo (Final AVG Confidence Level is
                        below those levels, or empty if de-novo made no prediction for the protein).
  de-novo only          Supported by de-novo but not by annotate. This includes proteins that
                        annotate classified as AVGs at a lower viral-origin confidence, and proteins
                        with no annotate AVG type (AVG Type is empty), which annotation-based
                        prediction alone would miss.

Columns in avg_predictions.tsv
------------------------------
Protein, Contig, Genome     Identifiers. For multi-contig bins/vMAGs, Genome is the bin name.
Evidence Class              See above.
AVG Type                    AMG, APG, or AReG, from annotate. Empty if annotate did not classify the
                            protein as an AVG.
Function                    From annotate. For AVGs, the description of the reference family that
                            supports the AVG call. Otherwise, the best-scoring HMM annotation, or empty
                            if the protein has no annotation.
Final AVG Confidence Level  See above.
Final AVG Probability       de-novo probability (0-1) underlying the final AVG confidence level.
Viral Origin Confidence (annotate)
                            See above.
AMG Weight, Functional Categories (Level 1), Functional Category Paths (L1 > L2 > L3),
Supporting HMM Database, Supporting HMM ID, Supporting HMM Bit Score, Protein in Strict Viral Region
                            From annotate, as described in the annotate README.txt. Empty for proteins
                            without an annotate AVG type. Functional categories are "Unknown" for
                            proteins with no annotation and "Other" for annotated proteins with no
                            assigned category.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
FAQ: https://github.com/AnantharamanLab/CheckAMG/wiki/FAQ
