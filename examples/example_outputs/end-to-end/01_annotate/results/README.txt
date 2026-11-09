CheckAMG annotate results
=========================

CheckAMG version: 1.2

This run
--------
Genomes analyzed:   116 (of 118 input genomes)
Proteins analyzed:  3,840
AVGs predicted:     70 (25 AMGs, 15 APGs, 30 AReGs)
  by viral-origin confidence: high 24, medium 20, low 26
AVG candidates removed by curation: 44 (see detailed/filtered_candidates.tsv)

Start here
----------
1. avg_predictions.tsv    One row per predicted auxiliary viral gene (AVG). This is the main result.
2. genome_summary.tsv     One row per genome, with AVG counts by type and viral-origin confidence.
3. faa/                   Amino-acid sequences of predicted AVGs: AMGs_*, APGs_*, AReGs_*, split by
                          viral-origin confidence (high, medium, low) plus *_all.faa.
4. detailed/              Supporting tables:
                            all_proteins.tsv         One row per analyzed protein, including proteins
                                                     that are not AVGs.
                            filtered_candidates.tsv  Proteins that matched an AVG reference but were
                                                     removed during curation, with the reason (one row
                                                     per protein and candidate AVG type).
                            run_summary.tsv          What was run (version, database, parameters), how
                                                     many inputs were analyzed or removed by the
                                                     length/ORF filters, and how many AVGs were found.
                          plus the full tables used to build the files above (annotations, genome
                          context features, curated tables, long-format categories, filter audits).

If this is part of a CheckAMG end-to-end run, the combined results of annotate and de-novo are in
03_aggregate/ (see 03_aggregate/README.txt).

AVG types
---------
AMG   auxiliary metabolic gene
APG   auxiliary physiology gene
AReG  auxiliary regulatory gene

Choosing a viral-origin confidence level
----------------------------------------
Every AVG is reported with one of three levels. Choose which to keep based on what your input is:
  high     Recommended for mixed-community metagenomes, or when input composition is unknown.
  medium   Can be added to high for viromes, confident predictions from virus-prediction tools,
           and isolate or complete viral genomes.
  low      Not recommended for any analysis. These proteins lack sufficient support for viral
           origin and may be encoded by host or other non-viral sequence.
Counting every row of avg_predictions.tsv without filtering by confidence overstates the number of
virus-encoded AVGs.

Columns in avg_predictions.tsv
------------------------------
Protein, Contig, Genome     Identifiers. For multi-contig bins/vMAGs, Genome is the bin name.
AVG Type                    AMG, APG, or AReG.
Function                    Description of the reference family that supports the AVG call. Among
                            the protein's HMM hits to reference families of its AVG type (for AMGs,
                            only families at or above the minimum AMG weight), the hit with the
                            highest bit score is used. A protein that qualifies as more than one
                            AVG type is reported once in this table and in detailed/all_proteins.tsv,
                            with priority AMG, then APG, then AReG (it still appears in each
                            matching FASTA file in faa/).
Viral Origin Confidence     high, medium, or low likelihood that the protein is virus-encoded,
                            from the genome-context LightGBM classifier. Refers to viral origin,
                            not to the functional annotation.
Viral Origin Probability    LightGBM probability (0-1) underlying the confidence level.
AMG Weight                  AMGs only. AMG weight (0-1) of the reference family that supplied
                            Function. It combines how specifically the family is associated with
                            metabolism (its metabolic ratio, from KEGG BRITE, FOAM, GO, and EC
                            classifications) with how rare it is in viral genomes (its VL-score).
                            Higher values mean more specifically metabolic and less common in
                            viruses. AMGs require at least the minimum AMG weight
                            (0.6 in this run, default 0.6).
Functional Categories (Level 1)
                            Top-level functional categories of the reference family that
                            supplied Function, separated by "; ". Useful for filtering.
Functional Category Paths (L1 > L2 > L3)
                            Every category assignment of that family as a full path, separated
                            by "; ". A family can belong to several pathways. This lists
                            potential pathways and does not mean the gene acts in all (or any)
                            of them. For counts per category, use the long-format
                            detailed/*_gene_categories.tsv files (one row per assignment).
Supporting HMM Database     Database of the HMM hit that supports the AVG call.
Supporting HMM ID           ID of that HMM.
Supporting HMM Bit Score    Bit score of that hit.
Protein in Strict Viral Region
                            True if the protein lies in a strict viral region, as identified by
                            CheckAMG annotate's V- and VL-score-based viral region algorithm (see
                            "What are strict viral regions?" in the wiki FAQ). This is an optional
                            stringency flag. Viral Origin Confidence is the primary viral-origin
                            call. Use this flag only for extra-stringent filtering (the
                            --filter-ambig-regions option applies it automatically).

Columns in detailed/all_proteins.tsv
------------------------------------
Protein Classification      metabolic (AMG), physiological (APG), regulatory (AReG), or
                            unclassified.
Annotation Status           AVG; AVG candidate, filtered (see detailed/filtered_candidates.tsv);
                            Annotated, not an AVG; or Unannotated.
Protein Viral Origin Confidence / Probability
                            As above, for every protein.
Protein in Strict Viral Region
                            As above.
Function                    For AVGs, chosen as described for avg_predictions.tsv. For all other
                            proteins, the best-scoring HMM annotation across all databases.
<database> columns          Best hit per database (KEGG, FOAM, Pfam, dbCAN/CAZy, METABOLIC,
                            CAMPER, PHROG).
Best Scoring HMM, Best Scoring HMM Annotation, Best Scoring HMM Origin
                            Highest-scoring hit across all databases. This can differ from the
                            hit supporting an AVG call.

Notes on other detailed/ tables
-------------------------------
*_score columns are HMM bit scores. *_coverage columns are the fraction of the HMM profile
covered by the alignment. *_V-score and *_VL-score columns are viral-likeness scores of the
annotation (higher = more virus-like). Distances are in base pairs.

Documentation: https://github.com/AnantharamanLab/CheckAMG/wiki
FAQ: https://github.com/AnantharamanLab/CheckAMG/wiki/FAQ
