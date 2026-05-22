# Phase 6 — Power diagnostics per pair × ancestry × outcome.
# Uses the Brion (2013) / Burgess (2014) two-sample-MR analytic formula.
# For each pair × ancestry slot with extracted instruments, computes the
# detectable OR at 80% power and the pre-registered α_Bonferroni.
# Pair × ancestry combinations failing the §19.2 power floor are flagged for
# descriptive-only reporting.
#
# Run from project root:
#   Rscript R/08_power_calc.R

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({ library(data.table) })

# --------------------------------------------------------------------------
# Detectable OR at 80% power (Brion 2013 / Burgess 2014, binary outcome)
# --------------------------------------------------------------------------
detectable_or <- function(r2, ncase, nctrl, power = 0.80, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_power <- qnorm(power)
  n_eff <- 4 * ncase * nctrl / (ncase + nctrl)
  if (r2 <= 0 || n_eff <= 0) return(NA_real_)
  b <- (z_alpha + z_power) / sqrt(r2 * n_eff)
  exp(b)
}

# --------------------------------------------------------------------------
# Compute total R² (variance-explained) for a given instrument set, from the
# per-SNP standardised effect estimate. For per-1-SD-standardised exposure,
# per-SNP R² ≈ 2·EAF·(1-EAF)·β². Sum over instruments gives total R².
# Fallback (when EAF is missing or zero): use F-stat-based estimate per SNP.
# --------------------------------------------------------------------------
compute_total_r2 <- function(instruments_dt) {
  dt <- copy(instruments_dt)
  dt[, EAF := as.numeric(EAF)]
  dt[, BETA := as.numeric(BETA)]
  dt[, SE  := as.numeric(SE)]
  dt[, N   := suppressWarnings(as.numeric(N))]

  dt[, r2_eaf := 2 * EAF * (1 - EAF) * BETA^2]

  # Fallback F-stat-based: R²_per_snp ≈ F / (F + N - 2)
  dt[, F_eff := F_stat]
  dt[is.na(F_eff), F_eff := (BETA / SE)^2]
  dt[, r2_fstat := F_eff / (F_eff + pmax(N - 2, 1, na.rm = TRUE))]

  # Take EAF estimate where present and sensible, otherwise F-stat estimate.
  dt[, r2_used := fifelse(!is.na(r2_eaf) & r2_eaf > 0 & r2_eaf < 1,
                          r2_eaf, r2_fstat)]
  sum(dt$r2_used, na.rm = TRUE)
}

# --------------------------------------------------------------------------
# Outcome catalogue (cases / controls) for the pre-specified outcomes
# --------------------------------------------------------------------------
outcomes <- data.table(
  outcome_id = c(
    "FinnGen_R12_I9_CHD", "FinnGen_R12_I9_AF",
    "FinnGen_R12_N14_ACUTERENFAIL", "FinnGen_R12_T2D",
    "FinnGen_R12_M13_GOUT", "FinnGen_R12_I9_CAVS_OPERATED",
    "FinnGen_R12_ST19_FRACT_FOREA",
    "AGEN_Spracklen2020_T2D",
    "MEGASTROKE_AS_EUR", "MEGASTROKE_AIS_EUR",
    "MEGASTROKE_LAS_EUR", "MEGASTROKE_CES_EUR", "MEGASTROKE_SVS_EUR",
    "MEGASTROKE_AIS_TRANS"),
  ncase = c(56650, 63532, 8383, 82878, 12342, 12418, 4439,
            77418, 40585, 34217, 4373, 7193, 5386, 34217),
  nctrl = c(443698, 252810, 480448, 403489, 315115, 487930, 463106,
            356122, 406111, 406111, 297290, 355468, 343560, 406111)
)

# --------------------------------------------------------------------------
# Iterate over all instrument files, compute R², then per-outcome power
# --------------------------------------------------------------------------
instrument_files <- list.files(paths$processed,
                               pattern = "^instruments_.*\\.csv$",
                               full.names = TRUE)

ncase_floor <- 3000L
or_floor_high <- 1.5  # two-sided threshold per §19.2

power_rows <- list()
for (f in instrument_files) {
  slot <- sub("^instruments_(.+)\\.csv$", "\\1", basename(f))
  insts <- fread(f)
  if (nrow(insts) == 0) next
  total_r2 <- compute_total_r2(insts)
  ancestry <- if (grepl("_EUR$", slot)) "EUR" else
              if (grepl("_EAS$", slot)) "EAS" else
              if (grepl("_AFR$|_AA$", slot)) "AFR" else "EUR"

  for (oc in seq_len(nrow(outcomes))) {
    det_or <- detectable_or(total_r2, outcomes$ncase[oc], outcomes$nctrl[oc])
    powered <- (!is.na(det_or) &
                outcomes$ncase[oc] >= ncase_floor &
                det_or <= or_floor_high)
    power_rows[[length(power_rows) + 1]] <- data.table(
      pair_slot = slot,
      ancestry = ancestry,
      n_snps = nrow(insts),
      total_r2 = round(total_r2, 5),
      outcome_id = outcomes$outcome_id[oc],
      ncase = outcomes$ncase[oc],
      nctrl = outcomes$nctrl[oc],
      detectable_or_at_80pct = round(det_or, 3),
      meets_power_floor = powered
    )
  }
}

power_dt <- rbindlist(power_rows)
out_path <- file.path(paths$outputs_tables, "power_diagnostics.csv")
fwrite(power_dt, out_path)
cat(sprintf("\nWrote %d power-diagnostic rows to %s\n",
            nrow(power_dt), out_path))

# --------------------------------------------------------------------------
# Headline summary
# --------------------------------------------------------------------------
cat("\n=== Per-instrument-set variance explained (R²) ===\n")
slot_r2 <- unique(power_dt[, .(pair_slot, ancestry, n_snps, total_r2)])
print(slot_r2[order(-total_r2)])

cat(sprintf("\n=== Slots meeting power floor (≥%d cases AND detectable OR ≤ %g) ===\n",
            ncase_floor, or_floor_high))
pass <- power_dt[meets_power_floor == TRUE,
                  .(pair_slot, outcome_id, total_r2,
                    detectable_or_at_80pct, ncase)]
print(pass[order(detectable_or_at_80pct)])

cat("\n=== Power-floor FAIL count by slot (across outcomes) ===\n")
fail_by_slot <- power_dt[meets_power_floor == FALSE, .N, by = .(pair_slot, ancestry)]
print(fail_by_slot[order(-N)])
