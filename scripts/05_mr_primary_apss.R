# Phase 3 — Primary MR using MR-APSS per pair × ancestry.
# Pre-registered as the primary method (plan §5.3, osf/registration.md §6) because
# MR-APSS explicitly models sample overlap and polygenicity.
#
# Run from project root:
#   Rscript R/05_mr_primary_apss.R
#
# Architecture:
#   1. For each pair, read raw genome-wide exposure + outcome sumstats.
#   2. MRAPSS::format_data() each (filters to HM3 panel, removes MHC).
#   3. MRAPSS::est_paras() with LDSC to get C (sample-structure intercept matrix)
#      and Omega (polygenic var-covar). LD score files: data/ld_panels/ldscores/<anc>/.
#   4. Build MRdat from our R/04 harmonized_*.rds (pre-registered IV set), looking
#      up per-SNP L2 from est_paras's Mdat to satisfy MRAPSS's expected column set.
#   5. MRAPSS::MRAPSS(MRdat, C, Omega, Cor.SelectionBias=TRUE) → beta, se, pval.
#
# Sumstats files are cached by path (the 6 LDL pairs share one EUR LDL file).
# Cache is bounded to {current exposure, current outcome} to avoid OOM.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(MRAPSS)
  library(data.table)
  library(jsonlite)
})

# Checkpoint — append one JSON line per completed pair so a mid-run crash
# preserves earlier work. Matches the pattern used by R/07.
CHECKPOINT_PATH <- file.path(paths$processed, "mrapss_checkpoint.jsonl")

read_checkpoint <- function() {
  # The JSONL is the "which pairs are done" gate ONLY — toJSON(auto_unbox=TRUE)
  # rounds doubles, so JSONL values lose precision (HDL_CAD_EUR pval 2.08e-5
  # → 0). The canonical full-precision source is the CSV from the previous
  # successful end-of-run write. We use that for `rows` if present.
  if (!file.exists(CHECKPOINT_PATH)) {
    return(list(completed = character(0), rows = list()))
  }
  lines <- readLines(CHECKPOINT_PATH, warn = FALSE)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) return(list(completed = character(0), rows = list()))
  jsonl_completed <- vapply(lines, function(l) fromJSON(l)$pair, character(1),
                            USE.NAMES = FALSE)
  csv_path <- file.path(paths$outputs_tables, "mr_apss_estimates.csv")
  if (file.exists(csv_path)) {
    full <- fread(csv_path)
    # A pair is only "done for resume purposes" if it has BOTH a JSONL line
    # AND a full-precision CSV row. JSONL-only pairs (from a prior run killed
    # before the end-of-script CSV write) get re-computed to preserve CSV
    # integrity. Otherwise the final CSV would silently drop those pairs.
    completed <- intersect(jsonl_completed, full$pair)
    re_run <- setdiff(jsonl_completed, full$pair)
    if (length(re_run)) {
      cat(sprintf("Resume: %d pair(s) in JSONL but missing from CSV — will recompute: %s\n",
                  length(re_run), paste(re_run, collapse = ", ")))
    }
    full <- full[pair %in% completed]
    rows <- split(full, by = "pair", keep.by = TRUE)
    cat(sprintf("Resume: %d JSONL completions, %d CSV-confirmed, %d full-precision rows loaded.\n",
                length(jsonl_completed), length(completed), nrow(full)))
    return(list(completed = completed, rows = unname(rows)))
  }
  completed <- jsonl_completed
  # Fallback: CSV missing (rare). Use lossy JSONL — log a warning.
  warning("Resuming from JSONL only; previously-completed pairs will have ",
          "rounded values until re-computed.")
  rows <- lapply(lines, function(l) as.data.table(fromJSON(l)))
  list(completed = completed, rows = rows)
}

append_checkpoint <- function(row_dt) {
  # toJSON on a data.table emits one object inside an array; we want bare object per line.
  obj <- as.list(row_dt[1, ])
  cat(toJSON(obj, auto_unbox = TRUE, na = "null"), "\n",
      file = CHECKPOINT_PATH, append = TRUE, sep = "")
}

# --------------------------------------------------------------------------
# Pair sources — same structure as R/07. Each entry maps a pair × ancestry
# label to its source-file paths + parser tags + N totals.
# --------------------------------------------------------------------------
pair_sources <- list(
  LDL_CAD_EUR       = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", ncase = 56650, nctrl = 443698),
  LDL_AIS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", ncase = 34217, nctrl = 406111),
  LDL_LAS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.3.LAS.EUR.out"),
                           ancestry = "EUR", ncase = 4373, nctrl = 297290),
  LDL_CES_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.4.CES.EUR.out"),
                           ancestry = "EUR", ncase = 7193, nctrl = 355468),
  LDL_SVS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.5.SVS.EUR.out"),
                           ancestry = "EUR", ncase = 5386, nctrl = 343560),
  HDL_CAD_EUR       = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", ncase = 56650, nctrl = 443698),
  TG_CAD_EUR        = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", ncase = 56650, nctrl = 443698),
  TC_CAD_EUR        = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", ncase = 56650, nctrl = 443698),
  BMI_T2D_EUR       = list(exp_tag = "GIANT_Yengo2018", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
                           ancestry = "EUR", ncase = 82878, nctrl = 403489),
  URATE_GOUT_EUR    = list(exp_tag = "CKDGen_Tin2019", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_M13_GOUT.gz"),
                           ancestry = "EUR", ncase = 12342, nctrl = 315115),
  SBP_AIS_EUR       = list(exp_tag = "ICBP_Evangelou2018", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_SBP.txt.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", ncase = 34217, nctrl = 406111),
  DBP_AIS_EUR       = list(exp_tag = "ICBP_Evangelou2018", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_DBP.txt.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", ncase = 34217, nctrl = 406111),
  LPA_CAVS_EUR      = list(exp_tag = "SinnottArmstrong2021", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "SinnottArmstrong_UKB/Lipoprotein_A.imp.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CAVS_OPERATED.gz"),
                           ancestry = "EUR", ncase = 12418, nctrl = 487930,
                           exp_N = 337000),  # Sinnott UKB Lp(a) effective N ≈ 337k (registration §7.1); per-variant N absent from source
  HBA1C_T2D_EUR     = list(exp_tag = "MAGIC_Chen2021", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EUR.tsv.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
                           ancestry = "EUR", ncase = 82878, nctrl = 403489),
  HBA1C_T2D_EAS     = list(exp_tag = "MAGIC_Chen2021", out_tag = "AGEN_T2D",
                           exp_path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EAS.tsv.gz"),
                           out_path = file.path(paths$outcomes, "AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz"),
                           ancestry = "EAS", ncase = 77418, nctrl = 356122),
  NC3_LDL_FRACT_EUR = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_ST19_FRACT_FOREA.gz"),
                           ancestry = "EUR", ncase = 4439, nctrl = 463106)
)

ld_score_dirs <- list(
  EUR = file.path(paths$ld_panels, "ldscores/EUR"),
  EAS = file.path(paths$ld_panels, "ldscores/EAS"),
  AFR = file.path(paths$ld_panels, "ldscores/AFR")
)

# --------------------------------------------------------------------------
# Sumstats cache (bounded to {current exp, current out}). Same pattern as R/07.
# --------------------------------------------------------------------------
.sumstats_cache <- new.env(parent = emptyenv())
cached_read <- function(path, tag) {
  key <- normalizePath(path, mustWork = FALSE)
  if (!is.null(.sumstats_cache[[key]])) {
    cat(sprintf("  [cache hit] %s\n", basename(path)))
    return(.sumstats_cache[[key]])
  }
  cat(sprintf("  [reading]   %s (tag=%s)\n", basename(path), tag))
  t0 <- Sys.time()
  dt <- read_sumstats(path, tag)
  cat(sprintf("              loaded %d rows in %.1fs\n",
              nrow(dt), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  .sumstats_cache[[key]] <- dt
  dt
}
evict_cache_except <- function(keep_paths) {
  keep_keys <- sapply(keep_paths, function(p) normalizePath(p, mustWork = FALSE))
  for (k in ls(.sumstats_cache)) {
    if (!k %in% keep_keys) rm(list = k, envir = .sumstats_cache)
  }
  invisible(gc(verbose = FALSE))
}

# --------------------------------------------------------------------------
# Adapt our read_sumstats output to MRAPSS::format_data's expected column names.
# Sample-size N is required by format_data; if NA in the source (Sinnott, FinnGen
# outcome), we substitute the registered total N from pair_sources.
# --------------------------------------------------------------------------
format_for_apss <- function(dt, ancestry, label, fallback_N = NULL) {
  # For sources whose SNP column is not rsIDs (ICBP MarkerName, Sinnott
  # chr:pos:ref:alt, AGEN-T2D, Pan-UKB), MRAPSS::format_data's HM3 snplist
  # merge retains zero SNPs and est_paras then dies with "0 (non-NA) cases".
  # Canonicalise to rsIDs via the ancestry-matched 1000G BIM first (same
  # pattern R/07 uses for the coloc merge). Helpers live in R/utils.R.
  ensure_chrpos(dt, ancestry, label = label)
  canonicalise_snp_to_rsid(dt, ancestry, label = label)
  dt <- dt[grepl("^rs[0-9]+$", dt$SNP)]
  df <- as.data.frame(dt)
  df <- df[!is.na(df$BETA) & !is.na(df$SE) & df$SE > 0 & !is.na(df$P), ]
  if (!is.null(fallback_N)) {
    df$N[is.na(df$N)] <- fallback_N
  }
  # Drop rows that still have NA N — format_data needs N for every row
  df <- df[!is.na(df$N), ]
  MRAPSS::format_data(
    dat = df,
    snp_col = "SNP", b_col = "BETA", se_col = "SE", freq_col = "EAF",
    A1_col = "EA", A2_col = "NEA", p_col = "P", n_col = "N",
    min_freq = 0.01  # slightly looser than the 0.05 default to retain Lp(a)-region variants
  )
}

# --------------------------------------------------------------------------
# Per-pair MR-APSS driver
# --------------------------------------------------------------------------
mr_apss_per_pair <- function(pair_label, ps) {
  cat(sprintf("\n[MR-APSS] %s — %s ancestry\n", pair_label, ps$ancestry))

  harm_path <- file.path(paths$processed,
                          sprintf("harmonized_%s.rds", pair_label))
  if (!file.exists(harm_path)) {
    cat(sprintf("  SKIP — harmonised RDS missing: %s\n", basename(harm_path)))
    return(NULL)
  }
  harm <- readRDS(harm_path)
  harm <- harm[harm$mr_keep, ]
  if (nrow(harm) < 5) {
    cat("  SKIP — fewer than 5 instruments after mr_keep filter.\n")
    return(NULL)
  }

  ld_dir <- ld_score_dirs[[ps$ancestry]]
  if (!dir.exists(ld_dir)) {
    cat(sprintf("  SKIP — LD score dir missing: %s\n", ld_dir))
    return(NULL)
  }

  # 1. Read raw sumstats (cached)
  evict_cache_except(c(ps$exp_path, ps$out_path))
  exp_dt <- cached_read(ps$exp_path, ps$exp_tag)
  out_dt <- cached_read(ps$out_path, ps$out_tag)

  # 2. Format for MRAPSS. Outcome cohorts may not report per-variant N; fall
  # back to the total case+ctrl from pair_sources.
  cat("  formatting exposure ...\n")
  exp_fmt <- format_for_apss(exp_dt, ancestry = ps$ancestry, label = "exp",
                              fallback_N = ps$exp_N)  # NULL for sources with per-variant N
  cat(sprintf("              %d HM3 SNPs retained\n", nrow(exp_fmt)))
  cat("  formatting outcome ...\n")
  out_fmt <- format_for_apss(out_dt, ancestry = ps$ancestry, label = "out",
                              fallback_N = ps$ncase + ps$nctrl)
  cat(sprintf("              %d HM3 SNPs retained\n", nrow(out_fmt)))

  # 3. Estimate background parameters via LDSC
  cat("  running est_paras (LDSC) ...\n")
  bg <- tryCatch(
    MRAPSS::est_paras(dat1 = exp_fmt, dat2 = out_fmt,
                      trait1.name = "exposure", trait2.name = "outcome",
                      LDSC = TRUE, ldscore.dir = ld_dir),
    error = function(e) {
      cat(sprintf("  est_paras FAILED: %s\n", conditionMessage(e))); NULL
    }
  )
  if (is.null(bg)) return(NULL)
  cat(sprintf("              C diag = (%.3f, %.3f); off-diag = %.3f\n",
              bg$C[1, 1], bg$C[2, 2], bg$C[1, 2]))
  cat(sprintf("              Omega diag = (%.2e, %.2e); off-diag = %.2e\n",
              bg$Omega[1, 1], bg$Omega[2, 2], bg$Omega[1, 2]))

  # 4. Build MRdat from our R/04 instruments + L2 lookup from bg$dat
  # (bg$dat is the LDSC-harmonised genome-wide table with L2 already attached;
  # see MRAPSS::est_paras body — returns list(dat=, ldsc_res=, C=, Omega=).)
  L2_lookup <- setNames(bg$dat$L2, bg$dat$SNP)
  MRdat <- data.frame(
    SNP      = harm$SNP,
    A1       = toupper(harm$effect_allele.exposure),
    A2       = toupper(harm$other_allele.exposure),
    b.exp    = harm$beta.exposure,
    b.out    = harm$beta.outcome,
    se.exp   = harm$se.exposure,
    se.out   = harm$se.outcome,
    pval.exp = harm$pval.exposure,
    pval.out = harm$pval.outcome,
    L2       = L2_lookup[harm$SNP],
    Threshold = 5e-8   # genome-wide IV selection per registration §11
  )
  retained_before <- nrow(MRdat)
  MRdat <- MRdat[!is.na(MRdat$L2), ]
  cat(sprintf("  MRdat: %d / %d instruments retained after L2 lookup\n",
              nrow(MRdat), retained_before))
  if (nrow(MRdat) < 5) {
    cat("  SKIP — fewer than 5 instruments after L2 lookup.\n")
    return(NULL)
  }

  # 5. MR-APSS
  cat("  fitting MRAPSS ...\n")
  fit <- tryCatch(
    MRAPSS::MRAPSS(MRdat = MRdat, exposure = "exposure", outcome = "outcome",
                   C = bg$C, Omega = bg$Omega, Cor.SelectionBias = TRUE),
    error = function(e) {
      cat(sprintf("  MRAPSS FAILED: %s\n", conditionMessage(e))); NULL
    }
  )
  if (is.null(fit)) return(NULL)

  # MRAPSS::MRAPSS may return beta/se/pvalue as character (formatted internally).
  # Coerce explicitly so downstream sprintf and data.table both behave.
  beta_n <- as.numeric(fit$beta)
  se_n   <- as.numeric(fit$beta.se)
  pval_n <- as.numeric(fit$pvalue)

  # Build the result row BEFORE logging — never lose a successful computation
  # because of a print-formatting error.
  res <- data.table(
    pair = pair_label, ancestry = ps$ancestry, method = "MR-APSS",
    b = beta_n, se = se_n, pval = pval_n, nsnp = nrow(MRdat),
    C_off_diag = bg$C[1, 2], C_diag1 = bg$C[1, 1], C_diag2 = bg$C[2, 2]
  )
  tryCatch(
    cat(sprintf("              beta=%.4f se=%.4f p=%.3g (n=%d)\n",
                beta_n, se_n, pval_n, nrow(MRdat))),
    error = function(e) cat("              (log format failed; result captured)\n")
  )
  res
}

# --------------------------------------------------------------------------
# Driver — iterate over all 16 pairs. Skip any pair already in the checkpoint;
# append each new result immediately so a mid-loop kill preserves earlier work.
# --------------------------------------------------------------------------
ckpt <- read_checkpoint()
if (length(ckpt$completed)) {
  cat(sprintf("Resuming — %d pair(s) already in checkpoint: %s\n",
              length(ckpt$completed), paste(ckpt$completed, collapse = ", ")))
}
out_rows <- ckpt$rows
for (pl in names(pair_sources)) {
  if (pl %in% ckpt$completed) {
    cat(sprintf("\n[MR-APSS] %s — already in checkpoint, skipping.\n", pl))
    next
  }
  res <- tryCatch(
    mr_apss_per_pair(pl, pair_sources[[pl]]),
    error = function(e) {
      cat(sprintf("  ERROR top-level: %s\n", conditionMessage(e))); NULL
    }
  )
  if (!is.null(res)) {
    out_rows[[length(out_rows) + 1]] <- res
    append_checkpoint(res)
  }
}

if (length(out_rows) > 0) {
  mrapss_results <- rbindlist(out_rows, fill = TRUE)
  out_path <- file.path(paths$outputs_tables, "mr_apss_estimates.csv")
  fwrite(mrapss_results, out_path)
  cat(sprintf("\nWrote %d MR-APSS rows to %s\n", nrow(mrapss_results), out_path))

  # Also merge into mr_estimates_per_method.csv (the canonical table consumed
  # by R/10 figures + the manuscript Results section). Replace any prior
  # "MR-APSS" rows so re-runs are idempotent.
  combined_path <- file.path(paths$outputs_tables, "mr_estimates_per_method.csv")
  if (file.exists(combined_path)) {
    existing <- fread(combined_path)
    existing <- existing[method != "MR-APSS"]
    # Keep only the columns shared with the existing table for the append.
    append_cols <- intersect(names(existing), names(mrapss_results))
    combined <- rbind(existing, mrapss_results[, ..append_cols], fill = TRUE)
    fwrite(combined, combined_path)
    cat(sprintf("Merged into %s (now %d rows total)\n",
                combined_path, nrow(combined)))
  }

  cat("\n=== MR-APSS results summary ===\n")
  print(mrapss_results[, .(pair, ancestry, b, se, pval, nsnp)])
} else {
  cat("\nNo MR-APSS rows produced.\n")
}
