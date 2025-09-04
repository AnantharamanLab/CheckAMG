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
