# Manuscript analysis notebooks

The Jupyter notebooks and R Markdown files in this directory produce the analyses, figures, and tables in the CheckAMG manuscript. Notebooks compute results and write tables to `tables/<analysis>/`. R Markdown files read those tables and write figures to `plots/<analysis>/` and figure source data to `tables/<analysis>/source_data/`. Scripts in `accessory_scripts/` are called from the notebooks.

Most notebooks read and write absolute paths on the analysis server. Expensive intermediate results are cached outside the repository to save runtime, and the notebooks skip a stage whose cached output exists. functional_propagation.ipynb builds the neighbor searches, sequence identities, and sequence clusterings that it and novel_avgs.ipynb read, and novel_avgs.ipynb runs the Foldseek searches. Large inputs and intermediate datasets are available as the Supplemental Data on Zenodo (DOI 10.5281/zenodo.22904215).

## Run order

Each file depends only on files above it. The database builds, protein embedding, model training, the benchmark and application runs, and the cached searches and clusterings in functional_propagation.ipynb and novel_avgs.ipynb are computationally expensive.

| step | file | produces |
|---|---|---|
| 1 | make_checkamg_required_tables.ipynb | CheckAMG reference tables in `CheckAMG/files/` (AMG, APG, and AReG tables, AMG weights, filters, categories) |
| 2 | build_annotate_db.ipynb | CheckAMG annotate database |
| 3 | train_test_datasets.ipynb | ground-truth datasets from geNomad and proGenomes |
| 4 | train_test_split.ipynb | training and test splits, `tables/lgbm/train_test_dataset_composition.tsv` |
| 5 | train_lgbm.ipynb | LightGBM viral-origin classifier and its evaluation tables in `tables/lgbm/` |
| 6 | lgbm_figures.Rmd | Extended Data Figures 2 and 3 |
| 7 | esm_embed_data.ipynb | ESM-2 protein embeddings |
| 8 | pst_embed_data.ipynb | PST graph-formatted embeddings |
| 9 | train_pst.ipynb | CheckAMG-PST fine-tuning |
| 10 | eval_pst.ipynb | CheckAMG-PST evaluation and thresholds in `tables/pst/` |
| 11 | pst_figures.Rmd | Figure 4 |
| 12 | build_denovo_db.ipynb | CheckAMG de-novo database |
| 13 | compare_pst_lgbm.ipynb | CheckAMG-PST versus LightGBM tables in `tables/compare_pst_lgbm/` |
| 14 | compare_pst_lgbm_figures.Rmd | Extended Data Figure 8 |
| 15 | strict_viral_regions.ipynb | strict viral region examples in `tables/genome_viz/` |
| 16 | strict_viral_regions_figures.Rmd | Extended Data Figure 1 |
| 17 | amg_benchmark_predictions.ipynb | CheckAMG, DRAM-V, VIBRANT, and geNomad runs on the benchmark datasets |
| 18 | amg_benchmark_data.ipynb | compiled benchmark predictions |
| 19 | amg_annotations.ipynb | compiled per-gene benchmark annotations |
| 20 | dramv_pfam_analysis.ipynb | DRAM-V Pfam analysis tables in `tables/DRAMV_pfam/` |
| 21 | dramv_pfam_analysis_figures.Rmd | Extended Data Figure 7 |
| 22 | amg_benchmark_comparisons.ipynb | benchmark comparison tables in `tables/amg_benchmark/` |
| 23 | amg_benchmark_figures.Rmd | Figure 3 and Extended Data Figures 4, 5, and 6 |
| 24 | amg_reference_weights_figures.Rmd | Figure 2 |
| 25 | run_pipeline.ipynb | CheckAMG annotate and de-novo runs on the MetaVR soil and human-gut genomes |
| 26 | functional_propagation.ipynb | functional label propagation, its evaluation, and the sequence clustering comparison in `tables/propagation/` |
| 27 | novel_avgs.ipynb | structural, DefenseFinder, and family analysis of the sequence-similarity invisible AVGs in `tables/novel_avgs/` |
| 28 | soil_gut_avgs.ipynb | soil and human-gut AVG tables in `tables/soil_gut_avgs/` |
| 29 | functional_propagation_figures.Rmd | Extended Data Figure 9 |
| 30 | soil_gut_avgs_figures.Rmd | Figures 5 and 6, Extended Data Figure 10, and the Supplemental Table 14 source |
