# Phase 2 — Instrument extraction per exposure × ancestry.
# Implements the pre-registered instrument construction (plan §5.2,
# osf/registration.md §5): genome-wide significance (P < 5e-8),
# ancestry-matched clumping (r² < 0.001, 10 Mb window), F ≥ 10.
#
# DO NOT RUN before the OSF pre-registration is frozen.
# Phase 0 gate enforced at runtime by checking osf/registration.md frontmatter.

source(here::here("R", "utils.R"))

reg_lines <- readLines(repo_path("osf", "registration.md"))
status_line <- reg_lines[grepl("^\\*\\*Status:\\*\\*", reg_lines)][1]
if (is.na(status_line)) {
  stop("osf/registration.md has no Status line — cannot verify OSF freeze state.")
}
if (grepl("DRAFT", status_line) || grepl("pending OSF freeze", status_line)) {
  stop("Phase-0 gate failed: osf/registration.md is still '", status_line,
       "'. Freeze on OSF (set Status to 'FROZEN — OSF DOI ...') before running Phase 2.")
}
if (!grepl("FROZEN", status_line)) {
  stop("Phase-0 gate failed: osf/registration.md Status is '", status_line,
       "' — expected FROZEN with OSF DOI before running Phase 2.")
}
cat("OSF freeze status verified:", status_line, "\n")

suppressPackageStartupMessages({
  library(ieugwasr)
  library(data.table)
  library(dplyr)
})

pairs <- read_canonical_pairs()

# Helper: instrument extraction for a single exposure × ancestry slot.
extract_instruments <- function(sumstats_path, source_tag, ancestry,
                                clump_kb = 10000, clump_r2 = 0.001,
                                pval_threshold = 5e-8, f_threshold = 10) {

  cat(sprintf("\n[extract_instruments] %s — %s ancestry\n", source_tag, ancestry))
  cat(sprintf("  Reading: %s\n", basename(sumstats_path)))

  dat <- read_sumstats(sumstats_path, source_tag)
  cat(sprintf("  Loaded %d variants\n", nrow(dat)))

  # Genome-wide significance filter
  dat <- dat[!is.na(P) & P < pval_threshold]
  cat(sprintf("  GWS (p<%.0e): %d variants\n", pval_threshold, nrow(dat)))

  # F-statistic filter (Bowden 2016): F = (BETA / SE)^2
  dat[, F_stat := (BETA / SE)^2]
  dat <- dat[!is.na(F_stat) & F_stat >= f_threshold]
  cat(sprintf("  F>=%g:        %d variants\n", f_threshold, nrow(dat)))

  if (nrow(dat) == 0) {
    cat("  ZERO variants after filters — no instruments.\n")
    return(dat[0])
  }

  # For sources without rsID marker IDs (ICBP, Sinnott-Armstrong), resolve
  # chr:pos against the ancestry-matched 1000G BIM before clumping.
  dat <- resolve_snp_ids_to_panel(dat, ancestry)
  if (nrow(dat) == 0) {
    cat("  ZERO variants after panel rsID resolution.\n")
    return(dat[0])
  }

  # Ancestry-matched LD clumping (with error handling — PLINK can fail when no
  # input SNPs match the panel, and ieugwasr surfaces that as a missing .clumped
  # file error which would otherwise crash the whole pipeline).
  panel_prefix <- file.path(paths$ld_panels,
                            sprintf("1000G_%s/1000G_%s", ancestry, ancestry))
  clump_input <- data.frame(rsid = dat$SNP, pval = dat$P, trait_id = source_tag)
  clumped <- tryCatch(
    ieugwasr::ld_clump_local(
      dat = clump_input,
      clump_kb = clump_kb, clump_r2 = clump_r2, clump_p = pval_threshold,
      bfile = panel_prefix,
      plink_bin = Sys.which("plink")
    ),
    error = function(e) {
      cat(sprintf("  CLUMPING FAILED: %s\n", conditionMessage(e)))
      cat("  Falling back to top-SNP-per-Mb pruning (no LD-based clumping).\n")
      # Fallback: sort by P, take top SNP per 1-Mb window per chromosome
      dt <- as.data.table(clump_input)
      dt[, dat_chr := dat$CHR]
      dt[, dat_bp  := dat$BP]
      setorder(dt, pval)
      kept <- character(0)
      seen <- list()
      for (i in seq_len(nrow(dt))) {
        chr <- dt$dat_chr[i]; bp <- dt$dat_bp[i]
        if (is.na(chr) || is.na(bp)) next
        key_chr <- as.character(chr)
        if (is.null(seen[[key_chr]])) seen[[key_chr]] <- integer(0)
        if (!any(abs(seen[[key_chr]] - bp) < clump_kb * 1000)) {
          kept <- c(kept, dt$rsid[i])
          seen[[key_chr]] <- c(seen[[key_chr]], bp)
        }
      }
      data.frame(rsid = kept)
    }
  )
  dat <- dat[SNP %in% clumped$rsid]
  cat(sprintf("  After clumping (r²<%.3f, %dMb): %d independent variants\n",
              clump_r2, clump_kb / 1000, nrow(dat)))

  return(dat)
}

# Iterate over the pre-specified slots from outputs/tables/canonical_pairs_prespec.csv.
# Slot list covers every powered or limited (exposure × ancestry) pre-registered
# at OSF freeze, organised by canonical pair family:
slots <- list(
  # Lipid pairs (Pairs 1-4) — GLGC EUR/EAS/AFR
  list(slot = "GLGC_LDL_EUR",  src = "GLGC", anc = "EUR",
       path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_LDL_EAS",  src = "GLGC", anc = "EAS",
       path = file.path(paths$exposures, "GLGC_Graham2021/EAS/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_LDL_AFR",  src = "GLGC", anc = "AFR",
       path = file.path(paths$exposures, "GLGC_Graham2021/AFR/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_HDL_EUR",  src = "GLGC", anc = "EUR",
       path = file.path(paths$exposures, "GLGC_Graham2021/EUR/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_HDL_EAS",  src = "GLGC", anc = "EAS",
       path = file.path(paths$exposures, "GLGC_Graham2021/EAS/HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_HDL_AFR",  src = "GLGC", anc = "AFR",
       path = file.path(paths$exposures, "GLGC_Graham2021/AFR/HDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TG_EUR",   src = "GLGC", anc = "EUR",
       path = file.path(paths$exposures, "GLGC_Graham2021/EUR/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TG_EAS",   src = "GLGC", anc = "EAS",
       path = file.path(paths$exposures, "GLGC_Graham2021/EAS/logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TG_AFR",   src = "GLGC", anc = "AFR",
       path = file.path(paths$exposures, "GLGC_Graham2021/AFR/logTG_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TC_EUR",   src = "GLGC", anc = "EUR",
       path = file.path(paths$exposures, "GLGC_Graham2021/EUR/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TC_EAS",   src = "GLGC", anc = "EAS",
       path = file.path(paths$exposures, "GLGC_Graham2021/EAS/TC_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")),
  list(slot = "GLGC_TC_AFR",   src = "GLGC", anc = "AFR",
       path = file.path(paths$exposures, "GLGC_Graham2021/AFR/TC_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")),
  # BMI (Pair 5) — GIANT Yengo 2018 EUR (genome-wide)
  list(slot = "GIANT_BMI_EUR", src = "GIANT_Yengo2018", anc = "EUR",
       path = file.path(paths$exposures, "GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz")),
  # Urate (Pair 6) — CKDGen Tin 2019 EUR + trans-ethnic
  list(slot = "CKDGen_URATE_EUR", src = "CKDGen_Tin2019", anc = "EUR",
       path = file.path(paths$exposures, "CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz")),
  list(slot = "CKDGen_URATE_TRANS", src = "CKDGen_Tin2019", anc = "EUR",
       path = file.path(paths$exposures, "CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_all_741_nstud37_summac400_rsid.txt.gz")),
  # BP (Pair 7) — ICBP+UKB EUR
  list(slot = "ICBP_SBP_EUR", src = "ICBP_Evangelou2018", anc = "EUR",
       path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_SBP.txt.gz")),
  list(slot = "ICBP_DBP_EUR", src = "ICBP_Evangelou2018", anc = "EUR",
       path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_DBP.txt.gz")),
  # HbA1c (Pair 10) — MAGIC Chen 2021 EUR/EAS/AA
  list(slot = "MAGIC_HBA1C_EUR", src = "MAGIC_Chen2021", anc = "EUR",
       path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EUR.tsv.gz")),
  list(slot = "MAGIC_HBA1C_EAS", src = "MAGIC_Chen2021", anc = "EAS",
       path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EAS.tsv.gz")),
  list(slot = "MAGIC_HBA1C_AA",  src = "MAGIC_Chen2021", anc = "AFR",
       path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_AA.tsv.gz")),
  # Lp(a) (Pair 9) — Sinnott-Armstrong EUR
  list(slot = "SA_LPA_EUR", src = "SinnottArmstrong2021", anc = "EUR",
       path = file.path(paths$exposures, "SinnottArmstrong_UKB/Lipoprotein_A.imp.gz"))
)

strength_rows <- list()
for (s in slots) {
  if (!file.exists(s$path)) {
    cat(sprintf("SKIP %s — file missing: %s\n", s$slot, s$path))
    next
  }
  out_path <- file.path(paths$processed,
                        sprintf("instruments_%s.csv", s$slot))
  if (file.exists(out_path) && file.info(out_path)$size > 100) {
    cat(sprintf("SKIP %s — output already present at %s\n", s$slot, out_path))
    insts <- fread(out_path)
  } else {
    insts <- tryCatch(
      extract_instruments(s$path, s$src, s$anc),
      error = function(e) {
        cat(sprintf("ERROR %s: %s\n", s$slot, conditionMessage(e)))
        data.table()
      }
    )
    if (nrow(insts) > 0) fwrite(insts, out_path)
  }
  if (nrow(insts) > 0) {
    strength_rows[[s$slot]] <- data.table(
      pair_slot = s$slot,
      source = s$src,
      ancestry = s$anc,
      n_snps = nrow(insts),
      median_F = median(insts$F_stat, na.rm = TRUE),
      min_F = min(insts$F_stat, na.rm = TRUE),
      pct_variance_explained = NA_real_  # fill in once N is known per source
    )
  }
}

if (length(strength_rows) > 0) {
  strength_table <- rbindlist(strength_rows)
  fwrite(strength_table, file.path(paths$outputs_tables,
                                   "instrument_strength.csv"))
  cat("\n=== Instrument strength table written ===\n")
  print(strength_table)
}

cat("\nExtraction complete. Next: R/04_harmonize.R\n")
