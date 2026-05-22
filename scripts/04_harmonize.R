# Phase 2 — Harmonisation per pair × ancestry, with the §13.1 Steiger filter
# as primary direction-of-causation check and the §12.1 APOE region exclusion
# for lipid → non-AD outcome pairs.
#
# Wraps TwoSampleMR::harmonise_data(action = 2) per pair × ancestry and writes
# one RDS per pair × ancestry to data/processed/.
#
# Run from project root:
#   Rscript R/04_harmonize.R

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

# APOE region in GRCh37 — excluded as primary for lipid → non-AD outcome pairs
# per registration §12.1. Full-genome sensitivity is computed alongside.
APOE_CHR    <- 19L
APOE_BP_LO  <- 44400000L
APOE_BP_HI  <- 46000000L

# --------------------------------------------------------------------------
# Pair × ancestry × outcome map.
# Each entry: pair label; exposure slot (matches R/03 output); outcome source
# tag (matches R/utils.R::read_sumstats); outcome file path; whether the APOE
# region should be excluded (TRUE for lipid → CAD / stroke pairs).
# --------------------------------------------------------------------------
pair_map <- list(
  # Pair 1: LDL-C → CAD
  list(pair = "LDL_CAD_EUR",       exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 56650, outcome_nctrl = 443698),

  # Pair 2: LDL-C → ischemic stroke (AIS primary)
  list(pair = "LDL_AIS_EUR",       exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 34217, outcome_nctrl = 406111),
  # Pair 2 subtype sensitivities
  list(pair = "LDL_LAS_EUR",       exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.3.LAS.EUR.out"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 4373,  outcome_nctrl = 297290),
  list(pair = "LDL_CES_EUR",       exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.4.CES.EUR.out"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 7193,  outcome_nctrl = 355468),
  list(pair = "LDL_SVS_EUR",       exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.5.SVS.EUR.out"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 5386,  outcome_nctrl = 343560),

  # Pair 3: HDL-C → CAD
  list(pair = "HDL_CAD_EUR",       exposure_slot = "GLGC_HDL_EUR",  outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 56650, outcome_nctrl = 443698),

  # Pair 4: Triglycerides → CAD
  list(pair = "TG_CAD_EUR",        exposure_slot = "GLGC_TG_EUR",   outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 56650, outcome_nctrl = 443698),

  # TC → CAD as sensitivity
  list(pair = "TC_CAD_EUR",        exposure_slot = "GLGC_TC_EUR",   outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 56650, outcome_nctrl = 443698),

  # Pair 5: BMI → T2D
  list(pair = "BMI_T2D_EUR",       exposure_slot = "GIANT_BMI_EUR", outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 82878, outcome_nctrl = 403489),

  # Pair 6: Urate → gout
  list(pair = "URATE_GOUT_EUR",    exposure_slot = "CKDGen_URATE_EUR", outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_M13_GOUT.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 12342, outcome_nctrl = 315115),

  # Pair 7: SBP → ischemic stroke
  list(pair = "SBP_AIS_EUR",       exposure_slot = "ICBP_SBP_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 34217, outcome_nctrl = 406111),
  list(pair = "DBP_AIS_EUR",       exposure_slot = "ICBP_DBP_EUR",  outcome_tag = "MEGASTROKE",
       outcome_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 34217, outcome_nctrl = 406111),

  # Pair 9: Lp(a) → calcific aortic valve stenosis (LPA cis-region instruments only — §12.2)
  list(pair = "LPA_CAVS_EUR",      exposure_slot = "SA_LPA_EUR",    outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CAVS_OPERATED.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 12418, outcome_nctrl = 487930),

  # Pair 10: HbA1c → T2D
  list(pair = "HBA1C_T2D_EUR",     exposure_slot = "MAGIC_HBA1C_EUR", outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 82878, outcome_nctrl = 403489),
  list(pair = "HBA1C_T2D_EAS",     exposure_slot = "MAGIC_HBA1C_EAS", outcome_tag = "AGEN_T2D",
       outcome_path = file.path(paths$outcomes, "AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 77418, outcome_nctrl = 356122),
  list(pair = "HBA1C_T2D_AFR",     exposure_slot = "MAGIC_HBA1C_AA", outcome_tag = "FinnGen_R12",
       # NOTE: AFR T2D outcome is DEFERRED (MVP). Using FinnGen as placeholder ONLY for pipeline-flow testing.
       # AFR Pair 10 reported descriptively only per registration §24 (power-floor consideration).
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
       exclude_apoe = FALSE, is_binary = TRUE,  outcome_ncase = 82878, outcome_nctrl = 403489),

  # NC-3: LDL-C → forearm fracture (negative control)
  list(pair = "NC3_LDL_FRACT_EUR", exposure_slot = "GLGC_LDL_EUR",  outcome_tag = "FinnGen_R12",
       outcome_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_ST19_FRACT_FOREA.gz"),
       exclude_apoe = TRUE,  is_binary = TRUE,  outcome_ncase = 4439,  outcome_nctrl = 463106)
)


# --------------------------------------------------------------------------
# Per-pair harmonisation helper
# --------------------------------------------------------------------------
harmonise_one_pair <- function(pm) {
  cat(sprintf("\n[harmonize] %s\n", pm$pair))

  inst_path <- file.path(paths$processed,
                         sprintf("instruments_%s.csv", pm$exposure_slot))
  if (!file.exists(inst_path)) {
    cat(sprintf("  SKIP — instruments not extracted: %s\n", inst_path))
    return(NULL)
  }
  if (!file.exists(pm$outcome_path)) {
    cat(sprintf("  SKIP — outcome file missing: %s\n", pm$outcome_path))
    return(NULL)
  }

  # Load instruments
  exp_dt <- fread(inst_path)
  n_pre <- nrow(exp_dt)

  # APOE-region exclusion (registration §12.1) for lipid → non-AD outcome pairs
  apoe_dropped <- 0L
  if (isTRUE(pm$exclude_apoe)) {
    in_apoe <- (exp_dt$CHR == APOE_CHR &
                  exp_dt$BP >= APOE_BP_LO & exp_dt$BP <= APOE_BP_HI)
    apoe_dropped <- sum(in_apoe, na.rm = TRUE)
    exp_dt <- exp_dt[!in_apoe]
    cat(sprintf("  APOE region excluded — %d SNP(s) dropped; %d remain\n",
                apoe_dropped, nrow(exp_dt)))
  }

  # Format exposure for TwoSampleMR
  exp_df <- TwoSampleMR::format_data(
    as.data.frame(exp_dt), type = "exposure",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "EA", other_allele_col = "NEA",
    eaf_col = "EAF", pval_col = "P", chr_col = "CHR", pos_col = "BP"
  )
  exp_df$id.exposure <- pm$exposure_slot
  exp_df$exposure    <- pm$exposure_slot

  # Read outcome and subset.
  # Primary match: SNP-by-SNP. Fallback: chr:pos (used when outcome SNP IDs are
  # not rsIDs — notably AGEN-T2D which uses `chr:pos_REF/ALT` format).
  out_dt <- read_sumstats(pm$outcome_path, pm$outcome_tag)
  out_sub <- out_dt[SNP %in% exp_df$SNP]
  if (nrow(out_sub) == 0 &&
        all(c("CHR", "BP") %in% names(out_dt)) &&
        all(c("chr.exposure", "pos.exposure") %in% names(exp_df))) {
    cat("  SNP-by-SNP match returned 0; attempting chr:pos fallback.\n")
    out_dt[, .key := paste0(CHR, ":", BP)]
    exp_keys <- paste0(exp_df$chr.exposure, ":", exp_df$pos.exposure)
    out_sub <- out_dt[.key %in% exp_keys]
    if (nrow(out_sub) > 0) {
      # Rewrite outcome SNP IDs to match the exposure (so harmonise_data() can join)
      key2snp <- setNames(exp_df$SNP, paste0(exp_df$chr.exposure, ":", exp_df$pos.exposure))
      out_sub[, SNP := key2snp[.key]]
      out_sub[, .key := NULL]
    }
  }
  cat(sprintf("  Outcome SNPs matched: %d / %d\n", nrow(out_sub), nrow(exp_df)))

  if (nrow(out_sub) < 2) {
    cat(sprintf("  SKIP — fewer than 2 outcome-matched SNPs after subset\n"))
    return(NULL)
  }

  out_df <- TwoSampleMR::format_data(
    as.data.frame(out_sub), type = "outcome",
    snp_col = "SNP", beta_col = "BETA", se_col = "SE",
    effect_allele_col = "EA", other_allele_col = "NEA",
    eaf_col = "EAF", pval_col = "P"
  )
  out_df$id.outcome <- pm$pair
  out_df$outcome    <- pm$pair

  # Harmonise (registration §11.1 — action = 2 drops palindromes with EAF [0.42, 0.58])
  harm <- TwoSampleMR::harmonise_data(
    exposure_dat = exp_df, outcome_dat = out_df, action = 2
  )

  n_post_harm <- sum(harm$mr_keep)
  cat(sprintf("  After action=2 harmonisation: %d SNPs (mr_keep=TRUE)\n", n_post_harm))

  # Steiger filter (registration §13.1 — PRIMARY direction-of-causation check).
  # SNPs where R²_outcome > R²_exposure are dropped from the primary instrument
  # set; the filtered count is reported regardless of whether it triggers the
  # >20% direction-uncertain flag.
  if (n_post_harm >= 2 && all(c("eaf.exposure","beta.exposure","se.exposure",
                                "beta.outcome","se.outcome") %in% names(harm))) {
    samplesize_exp <- if ("samplesize.exposure" %in% names(harm)) harm$samplesize.exposure
                      else max(exp_dt$N, na.rm = TRUE)
    samplesize_out <- pm$outcome_ncase + pm$outcome_nctrl
    n_cases        <- pm$outcome_ncase

    harm$samplesize.exposure <- samplesize_exp
    harm$samplesize.outcome  <- samplesize_out
    harm$ncase.outcome       <- n_cases
    harm$ncontrol.outcome    <- pm$outcome_nctrl
    harm$units.exposure      <- "SD"
    harm$units.outcome       <- if (isTRUE(pm$is_binary)) "log odds" else "SD"

    steiger_out <- tryCatch(
      TwoSampleMR::steiger_filtering(harm),
      error = function(e) {
        cat(sprintf("  Steiger filter ERROR: %s — Steiger flag not applied.\n",
                    conditionMessage(e)))
        NULL
      }
    )
    if (!is.null(steiger_out)) {
      harm <- steiger_out
      n_steiger_pass <- sum(harm$steiger_dir, na.rm = TRUE)
      n_total_kept   <- sum(harm$mr_keep, na.rm = TRUE)
      pct_pass <- ifelse(n_total_kept > 0, 100 * n_steiger_pass / n_total_kept, NA)
      cat(sprintf("  Steiger pass: %d/%d (%.1f%%)\n",
                  n_steiger_pass, n_total_kept, pct_pass))
      # Per registration §13.3: if Steiger pass < 80%, flag direction-uncertain.
      harm$direction_uncertain_flag <- !is.na(pct_pass) && pct_pass < 80
    }
  } else {
    cat("  Steiger filter skipped (insufficient data)\n")
  }

  # Save the harmonised RDS (primary instrument set, after harmonisation + Steiger)
  out_path <- file.path(paths$processed,
                        sprintf("harmonized_%s.rds", pm$pair))
  saveRDS(harm, out_path)
  cat(sprintf("  Saved %s\n", basename(out_path)))

  list(
    pair          = pm$pair,
    n_pre_apoe    = n_pre,
    apoe_dropped  = apoe_dropped,
    n_outcome_match = nrow(out_sub),
    n_post_harm   = n_post_harm,
    n_steiger_pass = if (exists("n_steiger_pass")) n_steiger_pass else NA_integer_,
    direction_uncertain = if ("direction_uncertain_flag" %in% names(harm))
                            unique(harm$direction_uncertain_flag)[1] else NA
  )
}


# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
summary_rows <- list()
for (pm in pair_map) {
  res <- tryCatch(
    harmonise_one_pair(pm),
    error = function(e) {
      cat(sprintf("  ERROR %s: %s\n", pm$pair, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(res)) summary_rows[[length(summary_rows) + 1]] <- res
}

cat("\n=== Harmonisation summary ===\n")
if (length(summary_rows) > 0) {
  sum_dt <- rbindlist(lapply(summary_rows, as.data.table), fill = TRUE)
  print(sum_dt)
  fwrite(sum_dt,
         file.path(paths$outputs_tables, "harmonisation_log.csv"))
  cat(sprintf("\nWrote harmonisation_log.csv (%d pairs)\n", nrow(sum_dt)))
}

cat("\nProceed to R/05_mr_primary_apss.R.\n")
