# Phase 6 — Cross-ancestry heterogeneity meta-analysis.
# Per pair: meta-regression via metafor::rma() with ancestry as moderator.
# Cochran's Q across ancestries. I² per pair. Bonferroni correction at pair
# level (10 tests). FDR-BH sensitivity (q < 0.05).
#
# Reads outputs/tables/mr_estimates_per_method.csv produced by R/05 and R/06.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
})

mr_path <- file.path(paths$outputs_tables, "mr_estimates_per_method.csv")
if (!file.exists(mr_path)) stop("MR estimates not found. Run R/05 and R/06 first.")
mr_all <- fread(mr_path)

# Use primary MR method (MR-APSS where available, else IVW) for the
# heterogeneity test. Pre-registered choice per plan §5.5.
primary <- mr_all[method %in% c("MR-APSS", "Inverse variance weighted"),
                  .(pair, ancestry, b, se, nsnp,
                    method_used = method[which.max(method == "MR-APSS")]),
                  by = .(pair_outcome = sub("_(EUR|EAS|AFR)$", "", pair))]

# For each canonical pair (without ancestry suffix), test ancestry heterogeneity.
# With 2 ancestries: use the closed-form pairwise Cochran's Q test
#   Q = (b1-b2)^2 / (se1^2 + se2^2);  Q ~ chi-square(1)  → I² = max(0, (Q-1)/Q)
# With ≥ 3 ancestries: use metafor::rma with categorical ancestry moderator.
hetero_rows <- list()
for (po in unique(primary$pair_outcome)) {
  sub_dt <- primary[pair_outcome == po]
  if (nrow(sub_dt) < 2) {
    cat(sprintf("Skip %s — only %d ancestry(ies) present.\n", po, nrow(sub_dt)))
    next
  }
  if (nrow(sub_dt) == 2) {
    b1 <- sub_dt$b[1]; b2 <- sub_dt$b[2]
    s1 <- sub_dt$se[1]; s2 <- sub_dt$se[2]
    Q <- (b1 - b2)^2 / (s1^2 + s2^2)
    Q_pval <- pchisq(Q, df = 1, lower.tail = FALSE)
    I2 <- max(0, (Q - 1) / Q * 100)
    test_kind <- "Pairwise Q (k=2)"
  } else {
    rma_out <- tryCatch(
      metafor::rma(yi = sub_dt$b, sei = sub_dt$se,
                   mods = ~ factor(sub_dt$ancestry),
                   method = "REML", control = list(maxiter = 200)),
      error = function(e) NULL
    )
    if (is.null(rma_out)) {
      cat(sprintf("Skip %s — metafor::rma failed.\n", po)); next
    }
    Q <- rma_out$QE; Q_pval <- rma_out$QEp; I2 <- rma_out$I2
    test_kind <- "Moderator Q (k≥3)"
  }
  hetero_rows[[length(hetero_rows) + 1]] <- data.table(
    pair = po,
    n_ancestries = nrow(sub_dt),
    test_kind = test_kind,
    Cochran_Q = Q,
    Cochran_Q_pval = Q_pval,
    I2_pct = I2,
    method_used = paste(unique(sub_dt$method_used), collapse = "|")
  )
}

hetero_dt <- rbindlist(hetero_rows)

# Multiple testing
hetero_dt[, pval_Bonferroni := pmin(Cochran_Q_pval * .N, 1)]
hetero_dt[, pval_BH := p.adjust(Cochran_Q_pval, method = "BH")]
hetero_dt[, sig_at_Bonferroni_05 := pval_Bonferroni < 0.05]
hetero_dt[, sig_at_BH_05         := pval_BH < 0.05]

out_path <- file.path(paths$outputs_tables, "heterogeneity_meta.csv")
fwrite(hetero_dt, out_path)
cat(sprintf("Wrote %d heterogeneity rows to %s\n", nrow(hetero_dt), out_path))

cat("\n=== Heterogeneity summary (primary H1 test result) ===\n")
print(hetero_dt[order(-I2_pct)])

n_sig_bonf <- sum(hetero_dt$sig_at_Bonferroni_05, na.rm = TRUE)
n_sig_bh   <- sum(hetero_dt$sig_at_BH_05, na.rm = TRUE)
n_total    <- nrow(hetero_dt)
cat(sprintf("\nH1 test: %d / %d pairs significant at Bonferroni 0.05 (target >= 2)\n",
            n_sig_bonf, n_total))
cat(sprintf("Sensitivity: %d / %d pairs significant at BH-FDR 0.05\n",
            n_sig_bh, n_total))

cat("\nGenerate per-pair ancestry-comparison forest plots in R/10_plots.R (or via paper draft).\n")
