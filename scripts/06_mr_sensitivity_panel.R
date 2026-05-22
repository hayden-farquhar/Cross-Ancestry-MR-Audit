# Phase 4 — Sensitivity panel per pair × ancestry.
# Runs IVW-RE, MR-Egger (intercept test), weighted median, weighted mode,
# MR-PRESSO (1000 dists), RadialMR. Appends to
# outputs/tables/mr_estimates_per_method.csv.
#
# Defensive: each method runs in its own tryCatch so a single failure does not
# poison the entire pair × ancestry result.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(MRPRESSO)
  library(RadialMR)
  library(data.table)
})

set.seed(64)

harmonized_files <- list.files(
  paths$processed,
  pattern = "^harmonized_.*\\.rds$", full.names = TRUE
)
if (length(harmonized_files) == 0) {
  stop("No harmonized_*.rds files in data/processed/. Run R/04 first.")
}

new_row <- function(pair, ancestry, method, b = NA_real_, se = NA_real_,
                    pval = NA_real_, nsnp = NA_integer_, note = "") {
  data.frame(
    pair = pair, ancestry = ancestry, method = method,
    b = b, se = se, pval = pval, nsnp = nsnp, note = note,
    stringsAsFactors = FALSE
  )
}

run_sensitivity <- function(harm_path) {
  pair <- sub("^harmonized_(.+)\\.rds$", "\\1", basename(harm_path))
  cat(sprintf("\n[sensitivity] %s\n", pair))
  harm <- readRDS(harm_path)
  harm_kept <- harm[harm$mr_keep, ]
  if (nrow(harm_kept) < 2) {
    cat("  <2 SNPs; skip pair.\n")
    return(NULL)
  }
  ancestry <- sub(".*_(EUR|EAS|AFR)$", "\\1", pair)
  if (!ancestry %in% c("EUR","EAS","AFR")) ancestry <- "EUR"
  n_kept <- nrow(harm_kept)
  rows <- list()

  # ---- IVW / Egger / weighted median / weighted mode (TwoSampleMR::mr) -----
  tryCatch({
    mr_methods <- c("mr_ivw", "mr_egger_regression",
                    "mr_weighted_median", "mr_weighted_mode")
    mr_out <- TwoSampleMR::mr(harm, method_list = mr_methods)
    for (i in seq_len(nrow(mr_out))) {
      rows[[length(rows) + 1]] <- new_row(
        pair, ancestry, mr_out$method[i],
        b = mr_out$b[i], se = mr_out$se[i],
        pval = mr_out$pval[i], nsnp = mr_out$nsnp[i]
      )
    }
    cat(sprintf("  IVW + Egger + median + mode: %d estimates\n", nrow(mr_out)))
  }, error = function(e) {
    cat(sprintf("  IVW/Egger/median/mode FAILED: %s\n", conditionMessage(e)))
  })

  # ---- MR-Egger intercept test ------------------------------------------
  tryCatch({
    egg <- TwoSampleMR::mr_pleiotropy_test(harm)
    if (is.data.frame(egg) && nrow(egg) >= 1) {
      rows[[length(rows) + 1]] <- new_row(
        pair, ancestry, "MR-Egger intercept",
        b = egg$egger_intercept[1], se = egg$se[1],
        pval = egg$pval[1], nsnp = n_kept
      )
      cat(sprintf("  Egger intercept: b=%g  p=%g\n",
                  egg$egger_intercept[1], egg$pval[1]))
    }
  }, error = function(e) {
    cat(sprintf("  Egger intercept FAILED: %s\n", conditionMessage(e)))
  })

  # ---- Within-ancestry IVW heterogeneity (Q + I²) -----------------------
  tryCatch({
    het <- TwoSampleMR::mr_heterogeneity(harm)
    for (i in seq_len(nrow(het))) {
      rows[[length(rows) + 1]] <- new_row(
        pair, ancestry, sprintf("Q (%s)", het$method[i]),
        b = het$Q[i], se = het$Q_df[i],
        pval = het$Q_pval[i], nsnp = n_kept,
        note = "b=Q, se=df"
      )
    }
  }, error = function(e) {
    cat(sprintf("  Heterogeneity FAILED: %s\n", conditionMessage(e)))
  })

  # ---- MR-PRESSO --------------------------------------------------------
  if (n_kept >= 4) {
    tryCatch({
      presso_out <- MRPRESSO::mr_presso(
        BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome", SdExposure = "se.exposure",
        OUTLIERtest = TRUE, DISTORTIONtest = TRUE,
        data = harm_kept, NbDistribution = 1000, SignifThreshold = 0.05
      )
      main <- presso_out$`Main MR results`
      if (is.data.frame(main) && nrow(main) >= 1) {
        rows[[length(rows) + 1]] <- new_row(
          pair, ancestry, "MR-PRESSO (raw)",
          b = main$`Causal Estimate`[1], se = main$Sd[1],
          pval = main$`P-value`[1], nsnp = n_kept
        )
      }
      if (is.data.frame(main) && nrow(main) >= 2) {
        n_outliers <- sum(!is.na(presso_out$`MR-PRESSO results`$`Outlier Test`$Pvalue) &
                            presso_out$`MR-PRESSO results`$`Outlier Test`$Pvalue < 0.05)
        rows[[length(rows) + 1]] <- new_row(
          pair, ancestry, "MR-PRESSO (outlier-corrected)",
          b = main$`Causal Estimate`[2], se = main$Sd[2],
          pval = main$`P-value`[2], nsnp = n_kept - n_outliers,
          note = sprintf("n_outliers=%d", n_outliers)
        )
      }
      # Global pleiotropy test
      gpv <- presso_out$`MR-PRESSO results`$`Global Test`$Pvalue
      if (!is.null(gpv)) {
        rows[[length(rows) + 1]] <- new_row(
          pair, ancestry, "MR-PRESSO Global",
          pval = gpv, nsnp = n_kept,
          note = "Global pleiotropy test"
        )
      }
      cat("  MR-PRESSO complete\n")
    }, error = function(e) {
      cat(sprintf("  MR-PRESSO FAILED: %s\n", conditionMessage(e)))
    })
  }

  # ---- Radial-MR --------------------------------------------------------
  tryCatch({
    radial_dat <- RadialMR::format_radial(
      BXG = harm_kept$beta.exposure, BYG = harm_kept$beta.outcome,
      seBXG = harm_kept$se.exposure, seBYG = harm_kept$se.outcome,
      RSID = harm_kept$SNP
    )
    radial_ivw <- RadialMR::ivw_radial(
      radial_dat, alpha = 0.05 / n_kept, weights = 3
    )
    co <- radial_ivw$coef
    if (is.matrix(co) || is.data.frame(co)) {
      est <- if ("Estimate" %in% colnames(co)) co[1, "Estimate"] else co[1, 1]
      se  <- if ("Std. Error" %in% colnames(co)) co[1, "Std. Error"] else co[1, 2]
      pv  <- if ("Pr(>|t|)" %in% colnames(co)) co[1, "Pr(>|t|)"]
             else if ("Pr(>|z|)" %in% colnames(co)) co[1, "Pr(>|z|)"]
             else NA_real_
      rows[[length(rows) + 1]] <- new_row(
        pair, ancestry, "Radial-IVW",
        b = est, se = se, pval = pv, nsnp = n_kept
      )
      n_radial_outliers <- if (!is.null(radial_ivw$outliers) &&
                                 is.data.frame(radial_ivw$outliers))
                              nrow(radial_ivw$outliers) else 0L
      cat(sprintf("  Radial-IVW: b=%g  p=%g  n_outliers=%d\n",
                  est, pv, n_radial_outliers))
      if (n_radial_outliers > 0) {
        fwrite(radial_ivw$outliers,
               file.path(paths$outputs_tables,
                         sprintf("radial_outliers_%s.csv", pair)))
      }
    }
  }, error = function(e) {
    cat(sprintf("  Radial-IVW FAILED: %s\n", conditionMessage(e)))
  })

  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

all_results <- list()
for (hf in harmonized_files) {
  res <- tryCatch(run_sensitivity(hf),
                  error = function(e) { cat(sprintf("  ERROR top-level: %s\n", conditionMessage(e))); NULL })
  if (!is.null(res)) all_results[[length(all_results) + 1]] <- res
}

if (length(all_results) > 0) {
  all_dt <- rbindlist(lapply(all_results, as.data.table), fill = TRUE)
  out_path <- file.path(paths$outputs_tables, "mr_estimates_per_method.csv")
  fwrite(all_dt, out_path)
  cat(sprintf("\nWrote %d MR estimates across %d pairs to %s\n",
              nrow(all_dt), length(all_results), out_path))
}

cat("\nProceed to R/05 (MR-APSS) or R/07 (coloc) or R/09 (heterogeneity).\n")
