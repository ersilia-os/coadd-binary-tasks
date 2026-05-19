# Scripts

Run in order. All scripts expect to be executed from the repo root with the `lazyqsar` environment.

| Script | Description |
|---|---|
| `00_prepare_configs.py` | Extract available organisms, strains, and assay types from raw CO-ADD and SPARK files. Outputs skeleton CSVs to `config/` for manual review. |
| `01_explore_coadd.py` | Summarise CO-ADD screening coverage (compound counts per pathogen × strain × assay). Produces coverage and strain-overlap plots. |
| `02_parse_coadd_inhibition.py` | Preprocess raw CO-ADD inhibition data: standardise SMILES, convert concentrations to µM, extract operators. |
| `03_binarise_coadd_inhibition.py` | Binarise inhibition data at 50/75/90% cutoffs per strain and pathogen. Produces per-strain and merged CSV files. |
| `04_parse_coadd_mic.py` | Preprocess raw CO-ADD MIC (dose-response) data: standardise SMILES, convert to µM, extract censored operators. |
| `05_binarise_coadd_mic.py` | Binarise MIC data at 10/25/50 µM cutoffs per strain and pathogen, handling censored measurements. |
| `06_coadd_analysis.py` | Extract cytotoxicity datasets (CC50, HC10); generate bioactivity and cytotoxicity scatter plots; produce the `selection_table.csv` that drives downstream modelling. |
| `07_lq_coadd.py` | Train LazyQSAR models for all selected CO-ADD datasets (inhibition + MIC, three cutoff tiers). Runs 5-fold CV and saves final models. |
| `08_lq_cytotox.py` | Train LazyQSAR models for CC50 and HC10 cytotoxicity datasets. Predicts on a ChEMBL sample and cross-correlates with bioactivity model predictions. |
| `09_explore_spark.py` | Overview of SPARK data: lab-contribution vs summary-file coverage, organism/assay compound counts, overlap with processed CO-ADD data. |
| `10_parse_spark_mic.py` | Parse and binarise SPARK curated MIC data (all organisms/strains) in the same format as CO-ADD scripts 04–05. |
| `11_predict_euopenscreen.py` | Predict the EU-OPENSCREEN library (~106k SMILES) against all CoADD models for overlapping organisms. Saves per-model JSON reports (y_hat, y_score, y_rank). |
