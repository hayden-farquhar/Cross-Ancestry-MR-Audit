# Raw data acquisition

This repository does not redistribute the ~57 GB of raw GWAS summary statistics used by the analysis pipeline. All inputs are publicly accessible. Below are the acquisition instructions, ordered by source.

After acquiring, the expected layout is:

```
data/
├── raw/
│   ├── exposures/
│   │   ├── GLGC_Graham2021/{EUR,EAS,AFR}/...
│   │   ├── GIANT_BMI/...
│   │   ├── CKDGen_Tin2019/...
│   │   ├── ICBP_Evangelou2018/...
│   │   ├── GSCAN_Saunders2022/...
│   │   ├── SinnottArmstrong_UKB/...
│   │   ├── MAGIC_Chen2021/...
│   │   └── PanUKB_AFR/...
│   └── outcomes/
│       ├── FinnGen_R12/...
│       ├── MEGASTROKE_Malik2018/...
│       └── AGEN_Spracklen2020/...
├── ld_panels/
│   ├── 1000G_EUR/{1000G_EUR.bed,bim,fam,afreq}
│   ├── 1000G_EAS/...
│   ├── 1000G_AFR/...
│   └── ldscores/{EUR,EAS}/...
└── processed/   (created by scripts/04_harmonize.R)
```

Verify presence and parseable-headers via `bash scripts/audit_data_presence.sh` (exit 0 = clean).

## LD reference panels — 1000 Genomes phase 3 PLINK (EUR / EAS / AFR)

- Zenodo deposit [6614170](https://doi.org/10.5281/zenodo.6614170)
- Download `1000G_EUR.tar.gz`, `1000G_EAS.tar.gz`, `1000G_AFR.tar.gz` and extract under `data/ld_panels/` so each ancestry has its own `1000G_*/` folder with `.bed`, `.bim`, `.fam`, and `.afreq` files.

## LDSC LD-score files for MR-APSS

- EUR HM3 baseline: Bulik-Sullivan 2015 HM3 EUR LD scores. Zenodo deposit [8182036](https://doi.org/10.5281/zenodo.8182036) (`eur_w_ld_chr.tar.gz`). Extract under `data/ld_panels/ldscores/EUR/`.
- EAS partitioned baseline: Zenodo deposit [10515792](https://doi.org/10.5281/zenodo.10515792) (`1000G_Phase3_EAS_baseline_v1.2_ldscores.tgz`). Extract the single-column LDSC files under `data/ld_panels/ldscores/EAS/`.
- AFR LDSC files: not required for the present 16-slot analytical state (no AFR slots are POWERED in the current analysis).

## Exposures

### Lipids — GLGC 2021 (Graham et al., *Nature*)

Per-ancestry inverse-rank-normalised summary statistics for LDL, HDL, TG (logTG), TC.

- EBI GWAS Catalog mirror (recommended; uniform schema across ancestries):
  - LDL EUR: GCST90239654, EAS: GCST90239655, AFR: GCST90239656
  - HDL EUR: GCST90239652, EAS: GCST90239651, AFR: GCST90239650
  - logTG EUR: GCST90239664, EAS: GCST90239663, AFR: GCST90239662
  - TC EUR: GCST90239676, EAS: GCST90239675, AFR: GCST90239674
- Direct download FTP base: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90239001-GCST90240000/`
- File naming convention used in this pipeline: `{LDL,HDL,logTG,TC}_INV_{EUR,EAS,AFR}_HRC_1KGP3_others_ALL.meta.singlevar.results.gz`. The pipeline auto-detects both the U-Michigan and EBI schemas (`rsID/CHROM/EFFECT_SIZE/pvalue_GC` vs `variant_id/chromosome/beta/p_value`).

### Urate — CKDGen 2019 (Tin et al., *Nat Genet*)

- Source: CKDGen consortium downloads page.
- Files used: `urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz` (trans-ethnic), and the EUR-only equivalent.

### BMI — GIANT 2018 (Yengo et al., *Hum Mol Genet*)

- Source: [GIANT consortium downloads](http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files).
- File: `GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz`.

### Blood pressure — ICBP + UKB 2018 (Evangelou et al., *Nat Genet*)

- EBI GWAS Catalog FTP for accessions [GCST006624](https://www.ebi.ac.uk/gwas/studies/GCST006624) (SBP), [GCST006629](https://www.ebi.ac.uk/gwas/studies/GCST006629) (PP), [GCST006630](https://www.ebi.ac.uk/gwas/studies/GCST006630) (DBP).
- File naming used: `Evangelou_30224653_SBP.txt.gz`, `Evangelou_30224653_DBP.txt.gz`.

### Smoking — GSCAN Phase 2 2022 (Saunders et al., *Nature*)

- Source: GSCAN dbGaP / consortium downloads. Per-ancestry stratified deposits (EUR / EAS / AFR; each ZIP contains the five GSCAN traits — only SmkInit and CigDay are used here).

### UKB biomarkers — Sinnott-Armstrong 2021 (*Nat Genet*)

- Source: [biobankengine.stanford.edu downloads](https://biobankengine.stanford.edu/).
- Files used (`{trait}.imp.gz`): `Lipoprotein_A.imp.gz`, `Urate.imp.gz`, `HbA1c.imp.gz`, `LDL_direct_adjstatins.imp.gz`, `HDL.imp.gz`, `Triglycerides.imp.gz`, `Cholesterol_adjstatins.imp.gz`. Note: Sinnott files do not report per-variant N; the pipeline supplies a per-source fallback N (Sinnott UKB Lp(a) N ≈ 337,000 per registration §7.1).

### HbA1c — MAGIC 2021 (Chen et al., *Nat Genet*)

- Source: [MAGIC consortium downloads](https://magicinvestigators.org/downloads/).
- Files used: `MAGIC1000G_HbA1c_{TA,EUR,EAS,AA,SAS,HISP}.tsv.gz`.

### Pan-UKB AFR per-phenotype

- Source: [Pan-UKB S3 bucket](https://pan.ukbb.broadinstitute.org/), per-phenotype manifest.
- Files used: Lp(a), urate, HbA1c, BMI, SBP. Each file contains all six ancestries' columns; AFR effects are extracted in script 03.

## Outcomes

### FinnGen R12

- GCS bucket: `gs://finngen-public-data-r12/summary_stats/release/`
- Manifest: `https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv`
- Per-endpoint URL pattern: `gs://finngen-public-data-r12/summary_stats/release/finngen_R12_{ENDPOINT}.gz`
- Endpoints used: `I9_CHD`, `I9_AF`, `N14_ACUTERENFAIL`, `T2D`, `M13_GOUT`, `I9_CAVS_OPERATED`, `ST19_FRACT_FOREA`.

### MEGASTROKE 2018 (Malik et al., *Nat Genet*)

- EBI GWAS Catalog mirror (recommended; bypasses the MEGASTROKE.org form that has been intermittently 500'ing): GCST005838 — GCST006910 series.
- Files used: AS, AIS, LAS, CES, SVS EUR + AIS trans-ancestry.

### AGEN-T2D 2020 (Spracklen et al., *Nature*)

- Source: AGEN consortium downloads via GWAS Catalog accession [GCST010118](https://www.ebi.ac.uk/gwas/studies/GCST010118).
- File: `SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz`.

## DEFERRED outcome cohorts (not in current analysis)

The following cohorts were pre-registered as DEFERRED pending Data Access Committee applications and could not enter the H1 cross-ancestry heterogeneity denominator post-hoc per the pre-registration §17.5. They are listed here for completeness; they are NOT required to reproduce the analyses in this repository.

- **BioBank Japan (BBJ)** — EAS CAD and ischemic stroke outcome cohorts. DAC application required via the [BBJ data-sharing process](https://biobankjp.org/en/).
- **Million Veteran Program (MVP)** — AFR CAD and gout outcome cohorts. Sponsoring U.S. institution required.
- **AWI-Gen** — AFR blood pressure outcome cohort. Sponsoring institution required.
- **DIAMANTE / T2DGGI** — EUR T2D additional cohort. Click-through agreement.

## Provenance log

Once downloads complete, every input file's SHA-256 hash, source URL, and download date is recorded by `scripts/01_fetch_summary_stats.R` to `data/processed/data_provenance.tsv`. The OSF deposit (DOI 10.17605/OSF.IO/U8TX4) also archives a snapshot log of the input fingerprints at pre-registration freeze.
