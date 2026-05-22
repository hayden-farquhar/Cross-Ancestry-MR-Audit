# Phase 6c (H2 secondary) — MR-APSS vs IVW divergence by sample-overlap class.
#
# Per registration §4.1 H2: mean |β_MR-APSS − β_IVW| in the high-overlap stratum
# should exceed the no-overlap stratum, with cluster-bootstrap CI on the ratio.
# H2 is interpreted DESCRIPTIVELY when <3 slots in either stratum (§4.1 closing
# sentence). At current data state we have 2 high-overlap slots — descriptive
# reporting only.
#
# Run from project root (after R/05 completes):
#   Rscript R/11_mr_apss_vs_ivw.R
#
# Output: outputs/tables/mr_apss_vs_ivw.csv + console summary.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(data.table)
})

# --------------------------------------------------------------------------
# Overlap classification per pair (registration §10.1 + §10.2 logic).
# FinnGen R12 is 0 % UKB → any EUR exposure paired against FinnGen is
# effectively no-overlap regardless of the EUR exposure's UKB fraction.
# MEGASTROKE has < 10 % UKB → EUR exposures with moderate UKB share land in
# the "moderate" stratum. ICBP + UKB BP → MEGASTROKE is the canonical
# high-overlap pair (both UKB-heavy).
# --------------------------------------------------------------------------
overlap_class <- c(
  LDL_CAD_EUR       = "no",       # GLGC EUR (~28% UKB) → FinnGen (0% UKB)
  LDL_AIS_EUR       = "moderate", # GLGC EUR → MEGASTROKE
  LDL_LAS_EUR       = "moderate",
  LDL_CES_EUR       = "moderate",
  LDL_SVS_EUR       = "moderate",
  HDL_CAD_EUR       = "no",
  TG_CAD_EUR        = "no",
  TC_CAD_EUR        = "no",
  BMI_T2D_EUR       = "no",       # GIANT EUR → FinnGen
  URATE_GOUT_EUR    = "no",       # CKDGen → FinnGen
  SBP_AIS_EUR       = "high",     # ICBP+UKB BP → MEGASTROKE (both UKB-heavy)
  DBP_AIS_EUR       = "high",
  LPA_CAVS_EUR      = "no",       # Sinnott UKB → FinnGen
  HBA1C_T2D_EUR     = "no",       # MAGIC → FinnGen
  HBA1C_T2D_EAS     = "no",       # MAGIC EAS → AGEN-T2D EAS
  NC3_LDL_FRACT_EUR = "no"        # GLGC EUR → FinnGen
)

# --------------------------------------------------------------------------
# Load the canonical per-method estimate table (must contain both IVW and
# MR-APSS rows after R/05 has merged its output in).
# --------------------------------------------------------------------------
combined_path <- file.path(paths$outputs_tables, "mr_estimates_per_method.csv")
stopifnot(file.exists(combined_path))
mr <- fread(combined_path)

stopifnot("pair" %in% names(mr))
stopifnot("method" %in% names(mr))
stopifnot("b" %in% names(mr))

apss <- mr[method == "MR-APSS", .(pair, ancestry, b_apss = b, se_apss = se, pval_apss = pval)]
ivw  <- mr[method == "Inverse variance weighted",
           .(pair, ancestry, b_ivw = b, se_ivw = se, pval_ivw = pval)]

if (nrow(apss) == 0) stop("No MR-APSS rows in ", combined_path,
                          " — run R/05 first.")
if (nrow(ivw) == 0) stop("No IVW rows in ", combined_path,
                         " — check method labels.")

div <- merge(apss, ivw, by = c("pair", "ancestry"), all = FALSE)
div[, abs_diff_b := abs(b_apss - b_ivw)]
div[, overlap := overlap_class[pair]]
div[, overlap := factor(overlap, levels = c("no", "moderate", "high"))]

# --------------------------------------------------------------------------
# Stratum summary
# --------------------------------------------------------------------------
strat <- div[, .(
  n_slots      = .N,
  mean_abs_diff = mean(abs_diff_b),
  median_abs_diff = median(abs_diff_b),
  max_abs_diff = max(abs_diff_b),
  pairs = paste(pair, collapse = "; ")
), by = overlap][order(overlap)]

cat("\n=== H2 descriptive: |β_MR-APSS − β_IVW| by overlap stratum ===\n")
print(strat)

# Per registration §4.1, H2 inference requires ≥3 slots per stratum.
# Report descriptively only if any stratum is below the floor.
high_n <- strat[overlap == "high", n_slots]
no_n   <- strat[overlap == "no", n_slots]
infer_ok <- length(high_n) > 0 && length(no_n) > 0 &&
            high_n >= 3 && no_n >= 3

if (!infer_ok) {
  cat(sprintf("\nDescriptive reporting only — high-overlap n = %d, no-overlap n = %d.\n",
              if (length(high_n)) high_n else 0L,
              if (length(no_n))   no_n   else 0L))
  cat("(Pre-registration §4.1: H2 inference suspended below 3 slots per stratum.)\n")
} else {
  # Cluster bootstrap on the ratio of mean abs diff (high : no-overlap)
  set.seed(2026)
  B <- 1000
  ratio_boot <- replicate(B, {
    high_boot <- div[overlap == "high"][sample(.N, .N, replace = TRUE), mean(abs_diff_b)]
    no_boot   <- div[overlap == "no"]  [sample(.N, .N, replace = TRUE), mean(abs_diff_b)]
    high_boot / no_boot
  })
  ratio_obs <- strat[overlap == "high", mean_abs_diff] / strat[overlap == "no", mean_abs_diff]
  ratio_ci  <- quantile(ratio_boot, c(0.025, 0.975))
  cat(sprintf("\nRatio mean|diff| high / no-overlap = %.3f (95%% CI %.3f–%.3f, %d boots)\n",
              ratio_obs, ratio_ci[1], ratio_ci[2], B))
  cat(sprintf("H2 supported = %s (criterion: ratio ≥ 1.5 AND 95%% CI excludes 1.0)\n",
              ifelse(ratio_obs >= 1.5 && ratio_ci[1] > 1.0, "TRUE", "FALSE")))
}

# --------------------------------------------------------------------------
# Persist per-slot divergence table
# --------------------------------------------------------------------------
out_path <- file.path(paths$outputs_tables, "mr_apss_vs_ivw.csv")
fwrite(div[, .(pair, ancestry, overlap, b_ivw, se_ivw, pval_ivw,
               b_apss, se_apss, pval_apss, abs_diff_b)], out_path)
cat(sprintf("\nWrote %d slot-level divergence rows to %s\n", nrow(div), out_path))

cat("\n=== Per-slot divergence ===\n")
print(div[order(overlap, -abs_diff_b),
          .(pair, ancestry, overlap, b_ivw = round(b_ivw, 3),
            b_apss = round(b_apss, 3), abs_diff_b = round(abs_diff_b, 3))])
