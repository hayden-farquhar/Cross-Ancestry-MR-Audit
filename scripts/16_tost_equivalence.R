# Phase 8 — Equivalence (TOST) test for the six "transport" pairs.
# Addresses reviewer major point #5: a non-significant cross-ancestry Cochran's Q
# is ABSENCE OF POWER, not evidence of equivalence. This recasts each non-urate
# pair as a two-one-sided-test (TOST) equivalence problem against a pre-specified
# relative margin, so each pair returns one of {equivalent, different, inconclusive}
# instead of a binary "transports / does not".
#
# Margin = m * |beta_EUR|, on the pre-registered relative-shift scale (registration
# §2.10): m = 0.25 is the "transport-discrepant" tier (the primary margin reported);
# m = 0.10 is the strict "transport-stable" tier (reported alongside).
#
# Reads:  outputs/tables/mr_apss_estimates.csv  (primary estimator = MR-APSS)
# Writes: outputs/tables/tost_equivalence.csv

source(here::here("R", "utils.R"))
suppressPackageStartupMessages(library(data.table))

apss <- fread(here::here("outputs", "tables", "mr_apss_estimates.csv"))
apss[, pair_outcome := sub("_(EUR|EAS|AFR)$", "", pair)]

# Cross-ancestry universe = stripped pairs with BOTH an EUR and an EAS arm.
# (URATE_GOUT_TRANS_EUR -> "URATE_GOUT_TRANS"; URATE_GOUT_EAS2 keeps its "2" and so
#  neither strips to a pairing partner — the sensitivity variants self-exclude.)
wide <- dcast(apss, pair_outcome ~ ancestry, value.var = c("b", "se"))
wide <- wide[!is.na(b_EUR) & !is.na(b_EAS)]

ALPHA <- 0.05
Z <- qnorm(1 - ALPHA)                       # 1.645 -> 90% CI for a 5% TOST

tost <- function(dB, seD, margin) {
  ci_lo <- dB - Z * seD; ci_hi <- dB + Z * seD
  # two one-sided tests against H0: |true diff| >= margin
  p_upper <- pnorm((dB - margin) / seD)             # evidence diff <  +margin
  p_lower <- pnorm((dB + margin) / seD, lower.tail = FALSE)  # evidence diff > -margin
  p_tost  <- max(p_upper, p_lower)
  equivalent <- (ci_lo > -margin) & (ci_hi < margin)
  different  <- (ci_lo > 0) | (ci_hi < 0)           # 90% CI excludes 0
  verdict <- if (equivalent) "equivalent"
             else if (different & (abs(dB) >= margin)) "different"
             else "inconclusive"
  list(ci_lo = ci_lo, ci_hi = ci_hi, margin = margin, p_tost = p_tost, verdict = verdict)
}

rows <- list()
for (i in seq_len(nrow(wide))) {
  w <- wide[i]
  dB  <- w$b_EAS - w$b_EUR
  seD <- sqrt(w$se_EAS^2 + w$se_EUR^2)
  for (m in c(0.25, 0.10)) {
    margin <- m * abs(w$b_EUR)
    r <- tost(dB, seD, margin)
    rows[[length(rows) + 1]] <- data.table(
      pair = w$pair_outcome, b_EUR = w$b_EUR, b_EAS = w$b_EAS,
      delta_beta = dB, se_delta = seD,
      ci90_lo = r$ci_lo, ci90_hi = r$ci_hi,
      margin_m = m, margin_abs = margin, p_tost = r$p_tost, verdict = r$verdict)
  }
}
res <- rbindlist(rows)
fwrite(res, here::here("outputs", "tables", "tost_equivalence.csv"))

cat("\n=== Cross-ancestry equivalence (TOST), margin = 0.25*|beta_EUR| ===\n")
cat("   (90% CI of EAS-EUR difference vs +/- margin; primary estimator MR-APSS)\n\n")
p <- res[margin_m == 0.25][order(verdict, pair)]
for (i in seq_len(nrow(p))) with(p[i], cat(sprintf(
  "  %-12s dB=%+.3f  90%%CI[%+.3f,%+.3f]  margin +/-%.3f  -> %s\n",
  pair, delta_beta, ci90_lo, ci90_hi, margin_abs, verdict)))
cat(sprintf("\n  m=0.25: equivalent=%d  inconclusive=%d  different=%d  (of %d pairs)\n",
    sum(p$verdict=="equivalent"), sum(p$verdict=="inconclusive"),
    sum(p$verdict=="different"), nrow(p)))
ps <- res[margin_m == 0.10]
cat(sprintf("  m=0.10: equivalent=%d  inconclusive=%d  different=%d\n",
    sum(ps$verdict=="equivalent"), sum(ps$verdict=="inconclusive"), sum(ps$verdict=="different")))
cat("\nHonest reading: pairs returning 'inconclusive' are UNDERPOWERED to demonstrate\n",
    "transport, not evidence of equivalence. Only 'equivalent' supports portability.\n")
