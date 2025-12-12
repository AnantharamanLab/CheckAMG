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
