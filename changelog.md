# 1.1.1

* Updated the database download script as was required for the 1.1 release
* Minor edits to code comments and log messages

# 1.1

* Modified the AMG weight calculation to slightly down-scale VL-scores in the formula (multiplied by 0.75)
* Added additional filtering keywords for AReGs and APGs
* Minor edits to `checkamg end-to-end` logging to point to module-specific log files when errors are encountered
* Updated the `de-novo` database to v1.1, containing newly-trained CheckAMG-PST model checkpoints and reference embeddings

# 1.0

This is a major release that adds annotation-independent AVG prediction with a protein genome language model, alongside the existing annotation-based approach, and the modules to combine and train them.

## New modules

* `de-novo`: predict auxiliary viral genes with a protein-based genome language model, the Protein Set Transformer (PST).
  * Reference-independent, deep-learning approach that takes query sequences (contigs, genomes, or proteins) and predicts whether each encoded protein is (A) viral and (B) auxiliary.
  * Runs on CPU or GPU (GPU recommended).
  * Accepts viral genomes/metagenomic contigs (.fna/.fasta) or translated proteins (.faa/.fasta, in Prodigal format), gzipped allowed, as single files or folders of bins. Also accepts pre-computed ESM2 or PST embeddings from previous runs on the same input to save time.
  * Embeds query proteins with PST (via ESM2), then compares them against the pre-trained reference embeddings with a nearest-neighbor search to produce, per protein, an `AVG-like prob` and a `Viral prob`, plus a combined `Final AVG prob`. Each is summarized into very-high-, high-, medium-, or low-confidence calls using thresholds calibrated on benchmark datasets.
  * Does not assign specific functions or classify predictions as metabolic, physiological, or regulatory on its own.

* `aggregate`: merge the results of `annotate` and `de-novo` into a single report combining functional classification, viral-origin confidence, and the annotation-independent AVG/viral predictions.

* `end-to-end`: run `annotate`, then `de-novo` (reusing the proteins predicted by `annotate`), then `aggregate`, writing each module to its own sub-directory (`01_annotate`, `02_denovo`, `03_aggregate`) in one output. `de-novo` can be directed to a GPU independently while `annotate` runs on CPU.

* `train`: finetune a CheckAMG-PST model to predict whether proteins are viral and auxiliary, for use with `de-novo`. Allows updating the released model or finetuning on user data with labeled viral/non-viral and AVG/non-AVG proteins. Trained models are supplied to `de-novo` via `--model-ckpt` with one of `--train-data-file`, `--train-embed-file`, or (`--train-index-file` and `--train-labels-file`).

## Databases

* Updated the `annotate` database to v1.1:
  * Updated Pfam-A from v38.0 to v38.2.
  * Updated KEGG KOfam from 2025-12-01 to 2026-02-01.
  * Added database-derived bitscore cutoffs for METABOLIC profile HMMs.
* Added a separate `de-novo` database (v1) containing the pre-trained CheckAMG-PST model checkpoint, reference protein embeddings, FAISS index, and labels for running `de-novo` with `--db-dir`.
* `checkamg download` now retrieves both databases by default and exposes `--download-annotate-db`/`--download-denovo-db` to fetch them selectively.

## Annotation filtering: AMG weight and reworked presets

* Added the **AMG weight** (`--min-amg-weight`, default 0.6) as the primary continuous filter for candidate AMGs, replacing the older binary hard/soft keyword scheme. It is a capped geometric mean of (i) the ratio of auxiliary-like to non-auxiliary-like metabolic annotations for a function and (ii) the function's rarity in viral genomes (VL-score). Candidates are kept only if their annotation's AMG weight meets the threshold.
* Reworked `--filter-presets`: renamed `allow_glycosyl` to `allow_glucan`, added `all_filters` (filter all `allow_*` categories, strict), and kept `no_filter`. The default now enables all four `allow_*` categories (`allow_glucan,allow_nucleotide,allow_methyl,allow_lipid`), so CheckAMG relies on the AMG weight plus essential-function filters by default while keeping glycosyl, nucleotide, methyl, and lipid annotations.
* Updated and expanded the curated reference AVG lists (`AMGs.tsv`, `APGs.tsv`, `AReGs.tsv`) and the AVG filter tables (`AMG_filters.tsv`, `APG_filters.tsv`, `AReG_filters.tsv`).

## Genome context and curation

* Strict viral regions are now identified from window-average VL-scores and refined using per-gene V-scores (or viral hallmark genes if `--use-hallmark`). Added the `--filter-ambig-regions` option to exclude predictions outside strict viral regions (off by default).
* Added the "Protein in Strict Viral Region" column to `annotate`'s `final_results.tsv`, indicating whether a protein falls within a strict viral region (one that would *not* be filtered under `--filter-ambig-regions`).
* AVG array filtering is now controlled by `--filter-avg-arrays` (on by default) and `--avg-array-limit` (default 3), replacing `--max-avg-array-length`.
* Retrained and updated the LightGBM viral-origin confidence model using the new relaxed V/VL-score cutoffs.
* Assigned reference AMGs, APGs, and AReGs to multi-level categories and report all applicable categories in `metabolic_gene_categories.tsv`, `physiology_gene_categories.tsv`, and `regulatory_gene_categories.tsv` under `results/`.

## HMM annotation changes in `annotate`

* Now enforces a minimum HMM-profile coverage **and** a minimum query-sequence coverage to keep hits, with an updated default of `0.4` for both.
* Raised the across-the-board minimum fallback bit score to `60` (`--bitscore`) for databases without defined cutoffs.
* Uses the newly obtained METABOLIC cutoffs, and enforces HMMsearch cutoffs separately for V/VL-score assignment (relaxed) and functional annotation (strict).
* Revised chunk-size calculations to produce more chunks and reduce memory usage on large inputs, plus other memory-usage improvements.

## CLI and workflow infrastructure

* Reorganized scripts into module-specific directories (`annotate`, `denovo`, `aggregate`, `end_to_end`, `download`, `train_pst`, `common`).
* Standardized logging across all scripts, including consistent step banners written cleanly to both console and log files.
* Consolidated shared helpers (Snakemake setup, memory limits, parquet paths, pyrodigal-gv, logger initialization) into `common` to reduce duplication.
* Corrected command logging so the full set of arguments is logged with their correct names, and made re-runnable commands accurate.
* Non-CheckAMG messages (errors from other libraries or base Python, and progress bars) are now appended to the log file alongside CheckAMG log messages.

## Input validation

* Added checks to the nucleotide and protein filtering steps to reject duplicate sequence names and avoid downstream mismatching.
* Added FASTA checks for duplicate sequence names, spaces in sequence names, and more robust verification that proteins are in Prodigal format, refactored into shared helper libraries under `common` for use by both `annotate` and `de-novo`.

## Bug fixes

* Fixed a crash in `annotate` during HMMsearch result filtering (`polars.exceptions.NoDataError`) that occurred when a search yielded no hits to a reference database. This happened with small inputs (< 100 proteins) lacking detectable homology to the smaller specialized databases (dbCAN, METABOLIC, CAMPER).
* Fixed a genome-context bug that reported V/VL-scores for nonexistent flanking MGE genes and sometimes matched flanking MGE genes to the wrong genes.
* Updated `hmm_annotate.py` HMM-profile-name retrieval for compatibility with pyhmmer > 0.11.


# 0.10.0

* Added genome-context curation to `curate_annots.py` module:
  * Flag proteins located outside strict viral regions or directly adjacent to their boundaries, based on window average VL-scores and V-scores of genes (or hallmark genes if enabled). These proteins are not removed from final AVG predictions by default, since strict viral region calls can be too conservative when viral origin confidence remains high.
  * Flag proteins in contiguous runs of 3 or more AVGs in a row, excluded by default as this indicates non-auxiliary function.
  * Parameter values are configurable with `--min-flank-vscore`, `--min-window-avg-vlscore`, and `--max-avg-array-length`.

* Updated the viral origin confidence LGBM
  * Previous versions trained and evaluated the model using train/test splits that could inadvertently exclude some proteins from a contig.
  * As a result, the model was sometimes trained with incomplete genome context information, though this did not introduce data leakage between training, validation, and test sets.
  * This has been corrected. The model was retrained and re-evaluated using datasets that retain all proteins encoded by each contig (including a new test dataset comprised of only host chromosomes with integrated proviruses and some integrated MGEs).
  * Inference is now parallelized by batching multiple contigs per prediction call to improve throughput.

* Additional changes:
  * Added support for gzipped FASTA file inputs (.fasta.gz, .fna.gz and .faa.gz).
  * Now logs how many AVGs of each type were filtered during the curation module.
  * Changed the order of some logging messages in `organize_proteins.py`.
  * Now writes full HMMsearch results to parquet instead of tsv when `--keep_full_hmm_results`is enabled.
  * Added additional methyltransferase annotations from Pfam to the AMG and AMG filter lists.
  * Added additional defense/anti-defense annotations to the APG and AReG lists
  * Added a feature to optionally save all intermediate and final tables to parquet instead of TSV to reduce filesize on large input datasets
  * Changed some argument names (e.g., `--genomes` to `--input-contigs`, `--vmags` to `--input-bins`) and renamed variables and log messages to be less specific to genomes/vMAGs and more generalized to contigs/bins

# 0.9.0

* Updated `checkamg download` to retrieve and extract a pre-built, standardized CheckAMG database containing all required profile HMMs and cutoff files, rather than downloading individual databases from their original sources.
  * This ensures reproducibility across CheckAMG and database versions and avoids failures caused by upstream download links changing or disappearing.
  * Added the notebook `build_checkamg_db.ipynb`, which documents how the standardized CheckAMG database is assembled, including data sources and formatting steps.

# 0.8.1

* Modified **annotate_hmm.py** so it resumes HMM searches from the last completed database instead of restarting all HMM searches across all databases.
  * This allows long runs with very large inputs that crash due to memory issues to resume where they left off, saving time when rerunning with more memory.
  * HMM search parameters, filtering strategy, and other aspects of the annotation pipeline are unchanged.

* Added the `--keep_full_hmm_results` option to CheckAMG annotate to control whether full HMM search results are written.
  * Previously, full results were always written by default, which can use substantial disk space for large inputs. This option now defaults to `False`.

# 0.8.0

* Removed the split between "hard" and "soft" keyword filters for AVG annotation filtering.
  * All keywords are now treated as a single filter set, including those previously classified as "soft".
  * These filter hits are no longer bypassed based on exceptional profile HMM matches.
  * As a result, the `--scaling_factor` argument has been removed.
  
* Genome context now reports the distance to contig ends for each gene

# 0.7.0

This release expands AVGs annotations, standardizes HMM annotation and filtering across databases, improves HMMsearch reporting and filtering, and adds reproducibility assets for rebuilding reference tables used by CheckAMG.

Major changes include:

* Expansion of the curated annotations used by CheckAMG (AMGs, APGs, AReGs), plus a large expansion of FOAM and KEGG reference annotations.   
* Added CAMPER profile HMMs ([McGivern et al., 2024](https://doi.org/10.1101/2023.09.24.559193)) to the CheckAMG database.
* Added reproducibility assets for rebuilding the required tables/files used by CheckAMG in the `notebooks` folder.
* KEGG AMGs were expanded using BRITE KO classifications (beyond the previous KOs sourced from VIBRANT). 
* False-positive filtering is now driven by explicit, standardized and inspected (see `make_checkamg_required_tables.ipynb`), pre-flagged HMM ID tables (hard/soft and exception categories) instead of only keyword lists. 
* Refined terms used to filter false-positives that were either too strict or lenient.
* HMMsearch reporting is more explicit: the pipeline now carries per-hit "kept vs removed" information (and rationale) and writes a best-hit-per-sequence filtered output in addition to the full hit table. 
* Default annotation thresholds were updated: `--scaling_factor` -> 3.0, `--bit_score` -> 30, `--cov_fraction` -> 0.30 (and `cov_fraction` is now  HMM profile coverage, not sequence coverage).  

# 0.6.2

* **AMGs.tsv**, **AReGs.tsv**, **hmm_id_to_name.csv**:
  * Fixed FOAM profile HMM annotations that had valid KO labels but did not get their names/descriptions properly mapped

# 0.6.1

* **filter_by_cds.py**:
  * Fixed a bug where vMAGs containing too few ORFs to pass the filter set by `--min_orf` were still written, but as empty files under `filtered_faa_by_cds`, causing the annotation step to crash

# 0.6.0

* **lgbm_model.joblib**, **lgbm_feature_names.joblib**, **lgbm_thresholds.joblib**:
  * Improvements to the viral origin confidence LGBM
  * Added additional features that consider the V/VL-scores of the 3 nearest genes on the left and right flanks of each gene, and the V/VL-scores of the nearest mobile genes

* **AMGs.tsv**, **APGs.tsv**, **AReGs.tsv**, **FOAM.tsv**, **hmm_id_to_name.csv**, **mobile_genes.csv**, **viral_hallmark_genes.csv**:
  * Slight modifications to the AMG, APG, AReG, viral hallmark, and mobile genes lists
  * Updated the all-HMM list with updated dbCAN and missing FOAM annotations

* **download\_db.py**:
  * Added functionality to download and prepare the database-provided bitscore thresholds for KEGG and FOAM using the same versions as the downloaded HMMs
  * Fixed incorrect version label for the dbCAN HMM files

* **CheckAMG_annotate.smk**:
  * Updated the KEGG and FOAM threshold file locations from the 'files' directory (with the AMG, APG, AReG, etc. tables) to the 'db' directory (with the HMM profiles)
  * Now puts the snakemake `*.done` files in their own folder

* **CheckAMG_annotate.py**:
  * Set up a subfolder for snakemake files

* **__main__.py**:
  * Changed the default `--scaling_factor` from `1.6` to `1.8` due to KEGG threshold updates

* **pyproject.toml**:
  * Added a missing `scikit-learn` dependency

# 0.5.3

* **download\_db.py**:
  * Updated URLs to download dbCAN v14

# 0.5.2

* Minor fix to the formatting in the ['false' AMG keywords](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/false_amgs.csv)

# v0.5.1

* Added annotation keyword filtering presets to allow certain functions to be retained when curating AVG predictions, if desired

# v0.5.0

* Improvements to the viral origin confidence LGBM
  * Added additional engineered features (min, max, n, log1p, and delta for V/VL-score and distance features)

* Updated the list of [MGE-like genes](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/mobile_genes.tsv)

* Removed overlapping HMMs across possible [AMGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AMGs.tsv), [APGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/APGs.tsv), and [AReGs](https://github.com/AnantharamanLab/CheckAMG/blob/main/CheckAMG/files/AReGs.tsv) and moved HMMs related to quorum sensing from APGs to AReGs

# v0.4.2

* Improvements to the viral origin confidence LGBM

# v0.4.1

* **__main__.py**
  * Changed the default `--window_size` from `25000` to `5000`
  
* Improvements to the viral origin confidence LGBM

# v0.4.0

* **annotate\_hmm.py**:
  * Changed parallelization from multiple concurrent processes (1 CPU each) to a single process using all available CPUs via pyhmmer's native threading.
  * Now splits input sequences into chunks only for very large inputs.
  * Reduces memory usage and improves speed for large datasets; similar runtime for smaller datasets.

* **filter\_by\_length.py**, **filter\_by\_cds.py**:
  * Revised parallelization to prevent bottlenecks from a few large sequences slowing all jobs.

* **check\_circular.py**:
  * Improved batching and parallelization, drastically reducing runtime when there are some very long sequences.

* **genome\_context.py**, **curate\_annots.py**:
  * Removed `max_flank_length` enforcement; no longer report binary flanking gene flags.
  * Instead, now report exact distances to the nearest V-score=10, viral hallmark, or MGE gene as features.
  * Stopped imputing missing features since LGBM  can handle `np.nan` directly.

* **__main__.py**, **CheckAMG\_annotate.py**, **CheckAMG\_annotate.smk**, **organize\_proteins.py**, **make\_final\_table.py**:
  * Minor compatibility updates for the above changes.

* Re-trained the viral origin confidence LGBM on a more robust dataset
  * Updated the model `.joblib` files for the new model
  * Updated the precision-recall curve plot and table
