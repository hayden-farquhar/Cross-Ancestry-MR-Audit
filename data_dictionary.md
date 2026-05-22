# Data dictionary

Variable definitions for every column in every CSV under `outputs/tables/`. Variables marked `[std]` use the project's standardised schema across all source GWAS, applied by the `read_sumstats()` parser in `scripts/utils.R`.

## Standardised GWAS-row schema (used internally)

| Column | Type | Description |
|--------|------|-------------|
| `SNP` `[std]` | character | rsID (rs[0-9]+) where available; chr:pos:ref:alt for non-rsID-native sources before canonicalisation |
| `CHR` `[std]` | integer | Chromosome number (1–22) |
| `BP` `[std]` | integer | Base-pair position (build varies by source — see §2.9 of accompanying manuscript) |
| `EA` `[std]` | character | Effect allele (A/C/G/T) |
| `NEA` `[std]` | character | Non-effect allele (A/C/G/T) |
| `EAF` `[std]` | numeric | Effect allele frequency in source population |
| `BETA` `[std]` | numeric | Per-allele effect estimate on natural-log scale (log-OR for binary outcomes) |
| `SE` `[std]` | numeric | Standard error of `BETA` |
| `P` `[std]` | numeric | Two-sided p-value (genomic-control-corrected where source uses GC) |
| `N` `[std]` | integer | Per-variant sample size (or fallback total-N where source omits per-variant N) |

## `canonical_pairs_prespec.csv`

The pre-specified pair × ancestry matrix, frozen at OSF pre-registration. One row per slot.

| Column | Description |
|--------|-------------|
| `pair_id` | Short identifier (e.g., `LDL_CAD`, `BMI_T2D`, `HBA1C_T2D`) |
| `exposure` | Exposure phenotype |
| `outcome` | Outcome phenotype |
| `ancestry` | EUR / EAS / AFR |
| `pair_class` | `canonical` (1–10) or `negative_control` (NC-1, NC-2, NC-3) |
| `pre_reg_status` | `POWERED` / `LIMITED` / `DEFERRED` / `POWER_FLOOR_FAIL` at pre-registration |
| `exposure_source` | Source GWAS for the exposure |
| `outcome_source` | Source GWAS for the outcome |
| `n_cases` / `n_controls` | Outcome case/control counts |
| `pre_reg_rationale` | One-line rationale from pre-registration §5.1 |

## `instrument_strength.csv`

Per-slot post-clumping instrument inventory before harmonisation.

| Column | Description |
|--------|-------------|
| `pair_slot` | Pair × ancestry identifier (e.g., `GLGC_LDL_EUR`) |
| `source` | Exposure GWAS source label |
| `ancestry` | EUR / EAS / AFR |
| `n_snps` | Count of GWS-clumped instruments retained |
| `median_F` | Median per-SNP F-statistic (Bowden 2016: F = (β/SE)²) |
| `min_F` | Minimum per-SNP F-statistic |
| `pct_variance_explained` | Cumulative R² explained by the instrument set in the exposure GWAS |

## `harmonisation_log.csv`

Per-pair-ancestry counts after `TwoSampleMR::harmonise_data(action=2)`, APOE region exclusion, and the Steiger direction-of-causation filter.

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `n_pre_apoe` | Instrument count before APOE region exclusion |
| `apoe_dropped` | Count dropped because they fall in chr19:44.4M–46.5M (GRCh37) |
| `n_outcome_match` | Count with a matching SNP in the outcome GWAS |
| `n_post_harm` | Count after harmonisation (palindrome handling + allele alignment) |
| `n_steiger_pass` | Count passing the per-SNP Steiger filter (Hemani 2017) |
| `direction_uncertain` | `TRUE` if < 80% of SNPs pass Steiger (registration §13.1 directional-confidence threshold) |

## `mr_estimates_per_method.csv`

One row per pair × ancestry × method × estimator type.

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `ancestry` | EUR / EAS / AFR |
| `method` | MR estimator name (`Inverse variance weighted`, `MR Egger`, `Weighted median`, `Weighted mode`, `MR-PRESSO (raw)`, `MR-PRESSO (outlier-corrected)`, `MR-PRESSO Global` (heterogeneity p), `MR-Egger intercept`, `Radial-IVW`, `Q heterogeneity`, `MR-APSS`) |
| `b` | Effect estimate on natural-log scale (log-OR for binary outcomes) |
| `se` | Standard error |
| `pval` | Two-sided p-value (stored as character for sub-double precision values like p < 1e-308) |
| `nsnp` | Instrument count used by this method (may differ from harmonised count for outlier-corrected estimators) |
| `note` | Optional note (e.g., MR-PRESSO outlier flags, Egger intercept value, etc.) |

## `mr_apss_estimates.csv`

The pre-registered primary estimator's per-slot results.

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `ancestry` | EUR / EAS / AFR |
| `method` | Always `MR-APSS` |
| `b` / `se` / `pval` | Effect estimate, SE, p (numeric) |
| `nsnp` | Instrument count after HM3 + L2 retention |
| `C_off_diag` | LDSC-derived C-matrix off-diagonal (empirical sample-overlap diagnostic; |C| > 0.05 flags non-trivial residual sample structure) |
| `C_diag1` / `C_diag2` | LDSC C-matrix diagonals (exposure / outcome variance-inflation factors) |

## `mr_apss_vs_ivw.csv`

H2 descriptive: per-slot MR-APSS vs IVW divergence, classified by pre-registered sample-overlap stratum.

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `ancestry` | EUR / EAS / AFR |
| `overlap` | `no` (< 5% shared sample), `moderate` (5–30%), or `high` (> 30%) — pre-classified at registration §10.2 |
| `b_ivw` / `se_ivw` / `pval_ivw` | IVW estimate |
| `b_apss` / `se_apss` / `pval_apss` | MR-APSS estimate |
| `abs_diff_b` | |β_MR-APSS − β_IVW|, the H2 divergence statistic |

## `coloc_results.csv`

`coloc.abf` posteriors per pair × top-3 instrument locus × prior parameterisation. 144 rows (16 pairs × 3 loci × 3 priors).

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `locus_index` | 1, 2, or 3 (top exposure loci by p-value) |
| `locus_lead_snp` | Lead SNP at the locus |
| `prior_class` | `primary` (p₁₂ = 1e-5), `stringent` (1e-6), or `liberal` (1e-4) |
| `H0` – `H4` | Posterior probabilities of hypotheses H0 (no variant), H1 (only exposure), H2 (only outcome), H3 (independent signals), H4 (shared causal variant) |
| `n_common_snps` | Number of SNPs in the 500-kb window common to both exposure and outcome GWAS |

## `power_diagnostics.csv`

Burgess 2014 two-sample MR power per slot.

| Column | Description |
|--------|-------------|
| `pair` | Pair × ancestry identifier |
| `n_cases` / `n_controls` | Outcome case / control counts |
| `r2_exposure` | Cumulative R² of the instrument set in the exposure GWAS |
| `detectable_OR_80pct` | Two-sided detectable OR at α = 0.05 / K_universe, 80% power |
| `power_at_OR_1.5` | Power to detect OR = 1.5 at the same α (sanity check) |
| `power_floor_pass` | `TRUE` if detectable_OR ∈ [1/1.5, 1.5] (registration §19.2) AND case count ≥ 3,000 |

## `heterogeneity_meta.csv`

Cross-ancestry Cochran's Q (H1 primary test) via `metafor::rma()` with ancestry as a categorical moderator, on IVW point estimates (sensitivity to the MR-APSS-based primary test in §3.4 of the manuscript).

| Column | Description |
|--------|-------------|
| `pair` | Pair identifier (cross-ancestry; no ancestry suffix) |
| `n_ancestries` | Number of ancestry arms contributing |
| `test_kind` | `Pairwise Q (k=2)` (the present analysis has only K_universe=1) |
| `Cochran_Q` | Q statistic |
| `Cochran_Q_pval` | p-value of Q under χ²(df = k-1) |
| `I2_pct` | I² percent of variation due to heterogeneity (clipped to 0 when Q < df) |
| `method_used` | The MR method whose estimates fed into Q (default `Inverse variance weighted`) |
| `pval_Bonferroni` / `pval_BH` | Multiple-comparisons corrections (trivial at K=1) |
| `sig_at_Bonferroni_05` / `sig_at_BH_05` | Decision flags |

## `mvmr_lipid_cad.csv`

Multivariable lipid MR (LDL + HDL + TG → CAD) using the Burgess 2015 weighted regression formulation.

| Column | Description |
|--------|-------------|
| `exposure` | `LDL` / `HDL` / `TG` |
| `estimator` | `MVMR-IVW` |
| `b_cond` | Conditional effect estimate (log-OR per 1-SD exposure, conditional on the other two lipids) |
| `se_cond` | SE of `b_cond` |
| `pval_cond` | Two-sided p-value |
| `n_snps` | Number of independent SNPs in the multivariable instrument set |

## `susie_coloc_lpa.csv`

SuSiE-based multi-variant colocalisation for the Lp(a) → CAVS LPA-region pair (§4.4 of accompanying manuscript). One row per (exposure credible set × outcome credible set) pair.

| Column | Description |
|---|---|
| `pair` | Always `LPA_CAVS_EUR` |
| `method` | Always `SuSiE-coloc` |
| `hit1` | Lead variant (rsID) of the exposure credible set |
| `hit2` | Lead variant (rsID) of the outcome credible set |
| `idx1` / `idx2` | Credible-set indices in the exposure / outcome SuSiE fits |
| `PP.H0.abf` – `PP.H4.abf` | Posterior probabilities of the five colocalisation hypotheses (no causal variant / only exposure / only outcome / two independent causal variants / one shared causal variant) |

## `mvmr_lipid_cad_harmonised.csv`

Per-instrument harmonised multi-exposure data underlying the MVMR fit (one row per instrument retained after r²<0.001 clumping + cross-lipid + outcome matching + LDL-effect-allele frame alignment).

| Column | Description |
|--------|-------------|
| `SNP` | rsID |
| `EA` / `NEA` | Effect allele / non-effect allele (aligned to LDL effect-allele frame) |
| `b_ldl` / `se_ldl` | LDL effect and SE |
| `b_hdl` / `se_hdl` | HDL effect and SE (aligned to LDL EA) |
| `b_tg` / `se_tg` | TG effect and SE (aligned to LDL EA) |
| `b_out` / `se_out` | CAD outcome effect and SE (aligned to LDL EA) |

## `radial_outliers_*.csv`

Per-pair Radial-MR outlier-flagged SNPs (one file per pair). Columns vary by pair but typically include `SNP`, `Q_i`, `Q_i_pval`, `outlier_status`.
