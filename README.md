# Cross-ancestry portability of canonical Mendelian randomization estimates

Analysis code repository for the pre-registered audit of ten canonical Mendelian randomization exposure–outcome pairs across European, East Asian, and African-ancestry populations.

**Author:** Hayden Farquhar MBBS MPHTM (Independent Researcher, Finley, NSW, Australia). ORCID [0009-0002-6226-440X](https://orcid.org/0009-0002-6226-440X).

**Pre-registration:** OSF DOI [10.17605/OSF.IO/U8TX4](https://doi.org/10.17605/OSF.IO/U8TX4) (frozen 2026-05-21).

**Preprint:** to be posted; this README will be updated with the preprint DOI on deposit.

**Zenodo archive:** concept DOI [10.5281/zenodo.20340928](https://doi.org/10.5281/zenodo.20340928) (cite-all-versions); v0.1.0 specific DOI [10.5281/zenodo.20340929](https://doi.org/10.5281/zenodo.20340929).

## Overview

The study is a pre-registered, multi-pair, multi-ancestry audit of canonical MR-derived causal claims (LDL → CAD, urate → gout, BMI → T2D, BP → stroke, Lp(a) → CAVS, HbA1c → T2D, lipid → stroke-subtype) using MR-APSS as the primary estimator alongside a seven-method sensitivity panel, colocalisation (`coloc.abf`), multivariable lipid MR, and pre-specified cross-ancestry heterogeneity testing.

This repository contains the analysis code that produces every quantitative finding, table, and figure reported in the accompanying manuscript. Manuscript-preparation code, submission scripts, and editorial-workflow tools are intentionally excluded — this is a reproducibility repository, not a publishing workflow.

## Data sources

All input data are publicly accessible GWAS summary statistics. The repository does **not** redistribute the raw files (combined ~57 GB). Acquisition instructions are in [`data/raw/README.md`](data/raw/README.md).

| Source | Phenotype(s) | Ancestries | Access |
|--------|--------------|------------|--------|
| GLGC 2021 (Graham et al., *Nature*) | LDL, HDL, TG, TC | EUR / EAS / AFR | Free — U Michigan + GWAS Catalog |
| CKDGen 2019 (Tin et al., *Nat Genet*) | Urate | EUR + trans-ethnic | Free |
| GIANT 2018 (Yengo et al., *Hum Mol Genet*) | BMI | EUR | Free |
| ICBP + UKB 2018 (Evangelou et al., *Nat Genet*) | SBP, DBP | EUR | Free — GWAS Catalog |
| GSCAN 2022 (Saunders et al., *Nature*) | Smoking | EUR / EAS / AFR | Free |
| Sinnott-Armstrong 2021 (*Nat Genet*) | Lp(a) and other UKB biomarkers | EUR (UKB) | Free |
| MAGIC 2021 (Chen et al., *Nat Genet*) | HbA1c | EUR / EAS / AFR | Free |
| FinnGen R12 | CAD, T2D, gout, CAVS, fracture | EUR (Finnish) | Free — GCS bucket |
| MEGASTROKE 2018 (Malik et al., *Nat Genet*) | Ischemic stroke + subtypes | EUR | Free — GWAS Catalog mirror |
| AGEN-T2D 2020 (Spracklen et al., *Nature*) | T2D | EAS | Free |
| Pan-UKB AFR | Lp(a), urate, HbA1c, BMI, SBP | AFR (UKB) | Free — S3 manifest |
| 1000 Genomes phase 3 PLINK panels | LD reference | EUR / EAS / AFR | Free — Zenodo 6614170 |
| Bulik-Sullivan LDSC LD score files | LD scores for MR-APSS | EUR HM3 baseline | Free |

Several additional data sources (BioBank Japan, Million Veteran Program CAD, AWI-Gen blood pressure, DIAMANTE/T2DGGI EUR T2D) were pre-registered as **DEFERRED** at pre-registration and did not enter the analysis. Reproducing their absence is part of the audit's substantive finding.

## Requirements

The pipeline is written in R with one Python figure-generation step. Tested on macOS 25.5 with R 4.4.2 and Python 3.11.

**R packages** (pinned via `renv.lock`; install via `Rscript scripts/00_setup_renv.R` once):

```
TwoSampleMR (0.6.4), MendelianRandomization (0.10.0), MRAPSS (0.1.2),
MRPRESSO (1.0), RadialMR (1.1), coloc (5.2.3), metafor (4.6-0),
ieugwasr (1.0.1), data.table (1.16.0), ggplot2, forestplot, here,
digest, openssl, jsonlite, dplyr, tidyr, readr.
```

**Python packages**:

```
matplotlib >= 3.7
pandas >= 2.0
streamlit >= 1.30   # for the optional dashboard
```

Install via:

```bash
pip install -r requirements.txt
```

**External binaries:** PLINK 1.9 (LD clumping) on `$PATH`.

## Reproduction

Steps below assume the data acquisition in [`data/raw/README.md`](data/raw/README.md) has been completed and `data/raw/`, `data/ld_panels/`, and `data/ld_panels/ldscores/` are populated. Expected total runtime end-to-end is **~6–8 hours** on a 16-core desktop (M-series Apple silicon or comparable); the LDSC-based MR-APSS step is the slowest single component (~30–60 min total).

```bash
# 0. Set up R environment (one-time; ~10 minutes)
Rscript scripts/00_setup_renv.R

# 1. Verify data presence (~30 seconds; non-zero exit if any required file missing)
bash scripts/audit_data_presence.sh

# 2. Acquire data fingerprints (~1 minute)
Rscript scripts/01_fetch_summary_stats.R

# 3. Verify ancestry-matched LD panels (~1 minute)
Rscript scripts/02_build_ld_panels.R

# 4. Extract instruments per pair × ancestry (~10 minutes)
Rscript scripts/03_extract_instruments.R

# 5. Harmonise instruments to outcomes (~15 minutes)
Rscript scripts/04_harmonize.R

# 6. Primary MR via MR-APSS (~30–60 minutes; reads 47M–8M-row sumstats per pair)
Rscript scripts/05_mr_primary_apss.R

# 7. Seven-method sensitivity panel (~15 minutes)
Rscript scripts/06_mr_sensitivity_panel.R

# 8. Colocalisation per pair × locus × prior (~20 minutes)
Rscript scripts/07_coloc.R

# 9. Power diagnostics + heterogeneity tests (~5 minutes)
Rscript scripts/08_power_calc.R
Rscript scripts/09_heterogeneity_meta.R

# 10. Per-pair forest plots (~2 minutes)
Rscript scripts/10_make_figures.R

# 11. H2 descriptive: MR-APSS vs IVW divergence by overlap stratum (~1 minute)
Rscript scripts/11_mr_apss_vs_ivw.R

# 12. Multivariable lipid MR for HDL → CAD (~15 minutes; reads 3 lipid GWAS)
Rscript scripts/12_mvmr_lipid_cad.R

# 13. Manuscript figures 1–3 (PRISMA flow, headline forest, H2 strip) (~30 seconds)
python3 scripts/13_produce_manuscript_figures.py
```

Optional: launch the interactive results dashboard.

```bash
streamlit run app/streamlit_app.py
```

## Script descriptions

| Script | Purpose | Outputs |
|--------|---------|---------|
| `scripts/00_setup_renv.R` | One-time R environment setup via renv | `renv.lock` |
| `scripts/01_fetch_summary_stats.R` | Verify all required exposure/outcome files present + checksummed | `outputs/tables/canonical_pairs_prespec.csv` |
| `scripts/02_build_ld_panels.R` | Verify 1000G phase 3 PLINK panels + LDSC files | — |
| `scripts/03_extract_instruments.R` | Per-ancestry GWS clumping → instrument list per slot | `outputs/tables/instrument_strength.csv` |
| `scripts/04_harmonize.R` | TwoSampleMR `harmonise_data(action=2)` + APOE exclusion + Steiger filter | 16 `harmonized_*.rds` + `outputs/tables/harmonisation_log.csv` |
| `scripts/05_mr_primary_apss.R` | Pre-registered primary estimator: MR-APSS per pair | `outputs/tables/mr_apss_estimates.csv` |
| `scripts/06_mr_sensitivity_panel.R` | 7-method sensitivity panel: IVW, Egger, median, mode, PRESSO, Radial | `outputs/tables/mr_estimates_per_method.csv` + `radial_outliers_*.csv` |
| `scripts/07_coloc.R` | `coloc.abf` per pair × top-3 locus × 3 priors | `outputs/tables/coloc_results.csv` |
| `scripts/08_power_calc.R` | Burgess 2014 detectable-OR floors per slot | `outputs/tables/power_diagnostics.csv` |
| `scripts/09_heterogeneity_meta.R` | Cross-ancestry Cochran's Q via metafor::rma | `outputs/tables/heterogeneity_meta.csv` |
| `scripts/10_make_figures.R` | Per-pair forest plots showing MR-APSS + sensitivity panel | 15 `forest_*.png` |
| `scripts/11_mr_apss_vs_ivw.R` | H2 descriptive: per-stratum MR-APSS vs IVW divergence | `outputs/tables/mr_apss_vs_ivw.csv` |
| `scripts/12_mvmr_lipid_cad.R` | Pre-registered multivariable MR (HDL\|TG,LDL → CAD) | `outputs/tables/mvmr_lipid_cad.csv` |
| `scripts/13_produce_manuscript_figures.py` | Figures 1 (PRISMA), 2 (headline forest), 3 (H2 strip) | 3 PNGs at `outputs/figures/figure_{1,2,3}_*.png` |
| `scripts/utils.R` | Shared helpers: paths, schema-aware sumstats readers, LD canonicalisation | — |
| `scripts/audit_data_presence.sh` | Verify every required input file exists + parses | exit 0 = clean |
| `tests/test_parsers.R` | Smoke tests for the 10 schema-aware sumstats parsers | console output |
| `app/streamlit_app.py` | Interactive results dashboard (optional) | served on localhost |

## Outputs and paper references

| Repository output | Manuscript element |
|------|----|
| `outputs/figures/figure_1_prisma_flow.png` | Figure 1 — PRISMA 2020 pair-selection flow |
| `outputs/figures/figure_2_headline_forest.png` | Figure 2 — Headline forest plots (LDL → LAS, HDL → CAD, HbA1c → T2D EAS) |
| `outputs/figures/figure_3_h2_strip.png` | Figure 3 — H2 sample-overlap-stratified IVW–MR-APSS divergence |
| `outputs/figures/forest_*.png` | Supplementary Figures S1–S15 (per-pair sensitivity-panel forests) |
| `outputs/tables/canonical_pairs_prespec.csv` | Supplementary Table S1 (pre-specified pair × ancestry) |
| `outputs/tables/instrument_strength.csv` | Supplementary Table S3 (per-slot instrument strength) |
| `outputs/tables/harmonisation_log.csv` | Supplementary Table S4 (post-harmonisation counts) |
| `outputs/tables/power_diagnostics.csv` | Supplementary Table S5 (per-slot detectable-OR) |
| `outputs/tables/coloc_results.csv` | Supplementary Table S6 (144-row coloc.abf posteriors) |
| `outputs/tables/mr_estimates_per_method.csv` | Supplementary Table S8 (per-method MR estimates) |
| `outputs/tables/mr_apss_estimates.csv` | Manuscript Table 2 (MR-APSS primary, 16 rows) |
| `outputs/tables/heterogeneity_meta.csv` | §3.4 (IVW-based H1 sensitivity) |
| `outputs/tables/mr_apss_vs_ivw.csv` | Manuscript Table 4 + §3.7 (H2 descriptive) |
| `outputs/tables/mvmr_lipid_cad.csv` | §3.3 multivariable lipid MR |

## Citation

If you use this code, please cite the accompanying manuscript and this repository:

```
Farquhar H. Cross-ancestry portability of canonical Mendelian randomization
estimates: a pre-registered audit reveals a structural evidence-base gap.
Preprint: to be posted. Code: https://doi.org/10.5281/zenodo.20340928
(concept DOI; cite-all-versions). Pre-registration: OSF DOI
10.17605/OSF.IO/U8TX4.
```

## License

Code is released under the **MIT License** (see [`LICENSE`](LICENSE)). Documentation, tables, and figures in this repository are released under **CC-BY 4.0**.

## Contact

For questions or to report issues with the analysis pipeline, open a GitHub issue or contact the author at `hayden.farquhar@icloud.com`.
