# Phase 5 — Colocalisation per pair × ancestry (registration §16).
# For each pair, identify the top-3 (lowest exposure p-value) instrument loci.
# Read 500-kb windows from both exposure and outcome GWAS at those loci.
# Run coloc.abf under three prior parameterisations (primary, stringent, liberal).
# Save per-locus posteriors H0–H4 to outputs/tables/coloc_results.csv.
#
# Run from project root:
#   Rscript R/07_coloc.R

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(coloc)
  library(data.table)
  library(jsonlite)
})

WINDOW_KB <- 500L
N_TOP_LOCI <- 3L
MIN_WINDOW_SNPS <- 30L  # below this, skip coloc for that locus

CKPT_PATH <- file.path(paths$processed, "coloc_checkpoint.jsonl")

# In-memory sumstats cache: keyed by absolute file path. The same GLGC EUR LDL
# file is the exposure for 6 pairs (LDL→CAD/AIS/LAS/CES/SVS/FRACT); without
# caching we re-fread the 3 GB file six times. Env-backed (hash) lookup.
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

# Evict any cached sumstats not in the keep_paths set. Called at the start of
# each pair to bound memory at ~2 files (1 exposure + 1 outcome) — full caching
# would hold ~20 large GWAS files in RAM (>30 GB) which OOMs on this machine.
evict_cache_except <- function(keep_paths) {
  keep_keys <- sapply(keep_paths, function(p) normalizePath(p, mustWork = FALSE))
  for (k in ls(.sumstats_cache)) {
    if (!k %in% keep_keys) {
      rm(list = k, envir = .sumstats_cache)
    }
  }
  invisible(gc(verbose = FALSE))
}

# Checkpoint helpers — write one JSON line per (pair, locus, prior) result.
# On crash, the loop can be resumed by re-running; pairs whose every locus×prior
# row is already in the checkpoint are skipped.
append_checkpoint <- function(rows_dt) {
  if (is.null(rows_dt) || nrow(rows_dt) == 0) return(invisible(NULL))
  con <- file(CKPT_PATH, open = "a")
  on.exit(close(con))
  for (i in seq_len(nrow(rows_dt))) {
    writeLines(toJSON(as.list(rows_dt[i]), auto_unbox = TRUE,
                      digits = NA, na = "null"), con)
  }
}

read_checkpoint <- function() {
  if (!file.exists(CKPT_PATH)) return(data.table())
  lines <- readLines(CKPT_PATH, warn = FALSE)
  if (length(lines) == 0) return(data.table())
  # Tolerate a truncated last line (crash mid-write) by dropping it on parse fail.
  parsed <- lapply(lines, function(ln) {
    tryCatch(fromJSON(ln), error = function(e) NULL)
  })
  parsed <- parsed[!sapply(parsed, is.null)]
  if (length(parsed) == 0) return(data.table())
  rbindlist(lapply(parsed, as.data.table), fill = TRUE)
}

pairs_done_in_checkpoint <- function(ckpt, pair_label) {
  if (nrow(ckpt) == 0) return(FALSE)
  any(ckpt$pair == pair_label)
}

# Fill missing CHR/BP via the 1000G BIM (rsID → chr:pos). MEGASTROKE in
# particular ships without chr/bp columns. Modifies dt in-place (data.table
# reference semantics) so the cache stays up-to-date for reuse across pairs.
ensure_chrpos <- function(dt, ancestry, label = "") {
  na_chr <- sum(is.na(dt$CHR))
  if (na_chr / nrow(dt) < 0.5) return(invisible(dt))
  cat(sprintf("  [chrpos resolve] %s: %d/%d rows lack chr/bp; resolving via 1000G %s BIM\n",
              label, na_chr, nrow(dt), ancestry))
  bim <- panel_chrpos_to_rsid(ancestry)
  m <- match(dt$SNP, bim$SNP)
  resolved <- !is.na(m)
  dt[resolved, CHR := bim$CHR[m[resolved]]]
  dt[resolved, BP  := bim$BP[m[resolved]]]
  cat(sprintf("                   resolved %d/%d rows\n", sum(resolved), nrow(dt)))
  invisible(dt)
}

# Replace non-rsID SNP column with rsID via CHR:BP → BIM lookup. Used for
# ICBP/Sinnott/AGEN/Pan-UKB whose SNP columns are `chr:pos:SNP` or
# `chr:pos:ref:alt`. Assumes source is on GRCh37 (true for all current
# non-rsID sources in this pipeline; 1000G phase 3 PLINK panels are GRCh37).
# Modifies dt in-place. rsIDs are build-agnostic, so this enables build-safe
# merging with FinnGen R12 (GRCh38) and other build-different outcomes.
canonicalise_snp_to_rsid <- function(dt, ancestry, label = "") {
  is_rs <- grepl("^rs[0-9]+$", dt$SNP)
  if (mean(is_rs, na.rm = TRUE) > 0.5) return(invisible(dt))
  cat(sprintf("  [snp canonicalise] %s: %d/%d rows non-rsID; resolving via 1000G %s BIM\n",
              label, sum(!is_rs), nrow(dt), ancestry))
  bim <- panel_chrpos_to_rsid(ancestry)
  needs <- !is_rs & !is.na(dt$CHR) & !is.na(dt$BP)
  dt[, key__ := paste0(CHR, ":", BP)]
  m <- match(dt$key__, bim$key)
  resolved <- needs & !is.na(m)
  dt[resolved, SNP := bim$SNP[m[resolved]]]
  dt[, key__ := NULL]
  cat(sprintf("                     resolved %d/%d rows\n",
              sum(resolved), sum(!is_rs)))
  invisible(dt)
}

# --------------------------------------------------------------------------
# Map each pair × ancestry to its exposure source, outcome source, and the
# parser tags needed for read_sumstats(). The harmonised RDS gives us SNP IDs
# and ancestry; we add the source-file paths for windowed reads.
# --------------------------------------------------------------------------
pair_sources <- list(
  LDL_CAD_EUR       = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 56650, nctrl = 443698),
  LDL_AIS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 34217, nctrl = 406111),
  LDL_LAS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.3.LAS.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 4373,  nctrl = 297290),
  LDL_CES_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.4.CES.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 7193,  nctrl = 355468),
  LDL_SVS_EUR       = list(exp_tag = "GLGC", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.5.SVS.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 5386,  nctrl = 343560),
  HDL_CAD_EUR       = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 56650, nctrl = 443698),
  TG_CAD_EUR        = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 56650, nctrl = 443698),
  TC_CAD_EUR        = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 56650, nctrl = 443698),
  BMI_T2D_EUR       = list(exp_tag = "GIANT_Yengo2018", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 82878, nctrl = 403489),
  URATE_GOUT_EUR    = list(exp_tag = "CKDGen_Tin2019", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_M13_GOUT.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 12342, nctrl = 315115),
  SBP_AIS_EUR       = list(exp_tag = "ICBP_Evangelou2018", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_SBP.txt.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 34217, nctrl = 406111),
  DBP_AIS_EUR       = list(exp_tag = "ICBP_Evangelou2018", out_tag = "MEGASTROKE",
                           exp_path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_DBP.txt.gz"),
                           out_path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 34217, nctrl = 406111),
  LPA_CAVS_EUR      = list(exp_tag = "SinnottArmstrong2021", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "SinnottArmstrong_UKB/Lipoprotein_A.imp.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CAVS_OPERATED.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 12418, nctrl = 487930),
  HBA1C_T2D_EUR     = list(exp_tag = "MAGIC_Chen2021", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EUR.tsv.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_T2D.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 82878, nctrl = 403489),
  HBA1C_T2D_EAS     = list(exp_tag = "MAGIC_Chen2021", out_tag = "AGEN_T2D",
                           exp_path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EAS.tsv.gz"),
                           out_path = file.path(paths$outcomes, "AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz"),
                           ancestry = "EAS", is_binary = TRUE, ncase = 77418, nctrl = 356122),
  NC3_LDL_FRACT_EUR = list(exp_tag = "GLGC", out_tag = "FinnGen_R12",
                           exp_path = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
                           out_path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_ST19_FRACT_FOREA.gz"),
                           ancestry = "EUR", is_binary = TRUE, ncase = 4439,  nctrl = 463106)
)


# Prior parameterisations per registration §16.2
PRIORS <- list(
  primary    = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5),
  stringent  = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-6),
  liberal    = list(p1 = 1e-4, p2 = 1e-4, p12 = 1e-4)
)


# --------------------------------------------------------------------------
# Pick top-N independent loci from instruments file (largest absolute Z; one
# lead per chromosome where multiple SNPs share a chr, taking the smallest p).
# --------------------------------------------------------------------------
top_loci <- function(inst_dt, n = N_TOP_LOCI) {
  dt <- copy(inst_dt)
  dt[, z := abs(BETA / SE)]
  setorder(dt, -z)
  # Take one lead per chromosome window (10 Mb apart) to spread loci
  kept <- list()
  for (i in seq_len(nrow(dt))) {
    chr <- dt$CHR[i]; bp <- dt$BP[i]
    too_close <- any(sapply(kept, function(k) k$CHR == chr && abs(k$BP - bp) < 1e7))
    if (!too_close) {
      kept[[length(kept) + 1]] <- list(SNP = dt$SNP[i], CHR = chr, BP = bp,
                                       P = dt$P[i], BETA = dt$BETA[i],
                                       SE = dt$SE[i])
      if (length(kept) >= n) break
    }
  }
  rbindlist(lapply(kept, as.data.table))
}


# --------------------------------------------------------------------------
# Run coloc per pair × ancestry. Reads both exposure and outcome sumstats once;
# extracts 500-kb windows in-memory for each top locus.
# --------------------------------------------------------------------------
run_pair_coloc <- function(pair_label, ps) {
  cat(sprintf("\n[coloc] %s\n", pair_label))

  # Identify lead loci from instruments
  inst_path <- file.path(
    paths$processed,
    sprintf("instruments_%s.csv",
            sub("_(EUR|EAS|AFR)$", paste0("_", ps$ancestry),
                # Find matching instrument file: pair-label-derived slot name
                # (the instrument-slot label is set at R/03 time and may differ
                # from the pair label; we use harmonised RDS to find SNPs)
                pair_label))
  )
  harm_path <- file.path(paths$processed,
                          sprintf("harmonized_%s.rds", pair_label))
  if (!file.exists(harm_path)) {
    cat(sprintf("  SKIP — harmonised RDS missing: %s\n", basename(harm_path)))
    return(NULL)
  }
  harm <- readRDS(harm_path)
  harm <- harm[harm$mr_keep, ]
  if (nrow(harm) < 2) {
    cat("  SKIP — < 2 SNPs in harmonised set.\n")
    return(NULL)
  }
  # Build the lead-loci table directly from the harmonised set
  lead_tbl <- data.table(
    SNP = harm$SNP, CHR = harm$chr.exposure, BP = harm$pos.exposure,
    P = harm$pval.exposure, BETA = harm$beta.exposure, SE = harm$se.exposure
  )
  lead_tbl <- lead_tbl[!is.na(CHR) & !is.na(BP) & !is.na(P)]
  loci <- top_loci(lead_tbl, n = N_TOP_LOCI)
  if (nrow(loci) == 0) {
    cat("  SKIP — no eligible lead loci with CHR/BP.\n")
    return(NULL)
  }
  cat(sprintf("  Top %d loci: %s\n", nrow(loci),
              paste(loci$SNP, collapse = ", ")))

  # Read exposure & outcome sumstats. The cache is bounded to {this exposure,
  # this outcome} — anything else gets evicted before reading.
  evict_cache_except(c(ps$exp_path, ps$out_path))
  exp_dt <- cached_read(ps$exp_path, ps$exp_tag)
  out_dt <- cached_read(ps$out_path, ps$out_tag)
  ensure_chrpos(exp_dt, ps$ancestry, label = "exp")
  ensure_chrpos(out_dt, ps$ancestry, label = "out")
  canonicalise_snp_to_rsid(exp_dt, ps$ancestry, label = "exp")
  canonicalise_snp_to_rsid(out_dt, ps$ancestry, label = "out")
  cat(sprintf("  exp rows: %d   out rows: %d\n", nrow(exp_dt), nrow(out_dt)))

  results <- list()
  for (i in seq_len(nrow(loci))) {
    locus_chr <- loci$CHR[i]
    locus_bp  <- loci$BP[i]
    win_lo <- locus_bp - WINDOW_KB * 1000
    win_hi <- locus_bp + WINDOW_KB * 1000

    # Exposure window: chr + bp on the exposure's own build (locus_bp comes from
    # the exposure side of harmonised data, so the BP range is valid).
    # Outcome window: chr only — exposure & outcome may be on different builds
    # (GLGC GRCh37 vs FinnGen R12 GRCh38), so BP comparison is unsafe. rsID
    # merge below restricts to the actual overlap.
    exp_win <- exp_dt[!is.na(CHR) & !is.na(BP) &
                       CHR == locus_chr & BP >= win_lo & BP <= win_hi]
    out_win <- out_dt[!is.na(CHR) & CHR == locus_chr]

    # Pre-merge QC: drop rows lacking BETA/SE/P or with non-rsID SNP; dedup by
    # SNP keeping min-P (multi-allelic variants share the rsID and would
    # otherwise trigger coloc.abf's "duplicated snps" rejection).
    exp_win <- exp_win[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P) &
                       grepl("^rs[0-9]+$", SNP)]
    out_win <- out_win[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(P) &
                       grepl("^rs[0-9]+$", SNP)]
    setorder(exp_win, P); exp_win <- exp_win[!duplicated(SNP)]
    setorder(out_win, P); out_win <- out_win[!duplicated(SNP)]

    common_snps <- intersect(exp_win$SNP, out_win$SNP)
    if (length(common_snps) < MIN_WINDOW_SNPS) {
      cat(sprintf("    locus %d %s: only %d common rsIDs (<%d); skip\n",
                  i, loci$SNP[i], length(common_snps), MIN_WINDOW_SNPS))
      next
    }
    exp_win <- exp_win[SNP %in% common_snps]
    out_win <- out_win[SNP %in% common_snps]
    common <- merge(exp_win, out_win, by = "SNP", suffixes = c(".exp", ".out"))
    if (nrow(common) < MIN_WINDOW_SNPS) {
      cat(sprintf("    locus %d %s: only %d SNPs after merge; skip\n",
                  i, loci$SNP[i], nrow(common)))
      next
    }

    # Determine N for the exposure (use the largest reported N, or fall back)
    exp_N <- if ("N.exp" %in% names(common))
               max(common$N.exp, na.rm = TRUE) else 1e5
    if (!is.finite(exp_N) || exp_N <= 0) exp_N <- 1e5

    # coloc.abf needs MAF + N for quant exposures to estimate sdY.
    # Drop rows with NA EAF before computing MAF; clamp to (0, 0.5).
    common <- common[!is.na(EAF.exp) & EAF.exp > 0 & EAF.exp < 1]
    if (nrow(common) < MIN_WINDOW_SNPS) {
      cat(sprintf("    locus %d %s: only %d SNPs with valid EAF; skip\n",
                  i, loci$SNP[i], nrow(common)))
      next
    }
    maf_exp <- pmin(common$EAF.exp, 1 - common$EAF.exp)

    # Build coloc inputs. SNP is the rsID (unique by construction after dedup).
    d_exp <- list(
      beta = common$BETA.exp, varbeta = common$SE.exp^2,
      snp  = common$SNP, position = common$BP.exp,
      type = "quant", N = exp_N, MAF = maf_exp
    )
    out_type <- if (isTRUE(ps$is_binary)) "cc" else "quant"
    s_case <- ps$ncase / (ps$ncase + ps$nctrl)
    d_out <- list(
      beta = common$BETA.out, varbeta = common$SE.out^2,
      snp  = common$SNP, position = common$BP.exp,
      type = out_type
    )
    if (out_type == "cc") {
      d_out$s <- s_case
    } else {
      d_out$N <- (ps$ncase + ps$nctrl)
    }

    primary_H4 <- NA_real_
    for (prior_name in names(PRIORS)) {
      pr <- PRIORS[[prior_name]]
      res <- tryCatch(
        suppressWarnings(coloc::coloc.abf(
          d_exp, d_out, p1 = pr$p1, p2 = pr$p2, p12 = pr$p12)),
        error = function(e) {
          cat(sprintf("    locus %d coloc.abf FAILED: %s\n",
                      i, conditionMessage(e)))
          NULL
        }
      )
      if (is.null(res)) next
      if (prior_name == "primary") primary_H4 <- res$summary[["PP.H4.abf"]]
      results[[length(results) + 1]] <- data.table(
        pair = pair_label, ancestry = ps$ancestry,
        locus_index = i, lead_snp = loci$SNP[i],
        chr = locus_chr, bp = locus_bp, prior = prior_name,
        H0 = res$summary[["PP.H0.abf"]], H1 = res$summary[["PP.H1.abf"]],
        H2 = res$summary[["PP.H2.abf"]], H3 = res$summary[["PP.H3.abf"]],
        H4 = res$summary[["PP.H4.abf"]], nsnps = res$summary[["nsnps"]]
      )
    }
    cat(sprintf("    locus %d %s (%s:%d±%dkb): %d positions; H4(primary)=%.3f\n",
                i, loci$SNP[i], as.character(locus_chr), locus_bp, WINDOW_KB,
                nrow(common), primary_H4))
  }

  if (length(results) == 0) return(NULL)
  rbindlist(results, fill = TRUE)
}


# --------------------------------------------------------------------------
# Driver — resumable from checkpoint
# --------------------------------------------------------------------------
ckpt <- read_checkpoint()
if (nrow(ckpt) > 0) {
  done_pairs <- unique(ckpt$pair)
  cat(sprintf("Checkpoint: %d rows across %d pairs already done (%s)\n",
              nrow(ckpt), length(done_pairs), paste(done_pairs, collapse = ", ")))
}

# Reorder pair_sources to group by exposure path (so the exposure-file cache hits
# maximally on contiguous pairs). Within each group, no ordering preference.
pair_order <- order(sapply(pair_sources, function(p) p$exp_path))
ordered_labels <- names(pair_sources)[pair_order]

for (pl in ordered_labels) {
  if (pairs_done_in_checkpoint(ckpt, pl)) {
    cat(sprintf("[skip] %s — already in checkpoint\n", pl))
    next
  }
  res <- tryCatch(
    run_pair_coloc(pl, pair_sources[[pl]]),
    error = function(e) {
      cat(sprintf("  ERROR top-level: %s\n", conditionMessage(e))); NULL
    }
  )
  # Persist incrementally so a crash mid-loop preserves earlier pairs.
  append_checkpoint(res)
}

# Final assembly: read all checkpoint lines and emit the CSV.
final <- read_checkpoint()
if (nrow(final) > 0) {
  out_path <- file.path(paths$outputs_tables, "coloc_results.csv")
  fwrite(final, out_path)
  cat(sprintf("\nWrote %d coloc rows (across %d pairs) to %s\n",
              nrow(final), length(unique(final$pair)), out_path))
  cat("\n=== Coloc H4 (primary prior) summary per pair ===\n")
  print(final[prior == "primary",
              .(.N, mean_H4 = mean(H4, na.rm = TRUE),
                max_H4 = max(H4, na.rm = TRUE)),
              by = pair])
} else {
  cat("No coloc results produced.\n")
}
