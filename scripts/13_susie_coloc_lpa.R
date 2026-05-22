# Phase 7d — SuSiE-based multi-variant colocalisation for the LPA region.
# Pre-registered §16 sensitivity for the Lp(a) → CAVS pair, motivated by the
# well-described multi-causal-variant architecture of the LPA region that
# violates coloc.abf's single-causal-variant assumption (manuscript §4.4).
#
# Run from project root:
#   Rscript R/13_susie_coloc_lpa.R > logs/13_susie_lpa_$(date +%Y%m%d_%H%M%S).log 2>&1
#
# Window: chr6:160,500,000-161,100,000 (GRCh37) — covers LPA gene + flanking
# variants. The KIV-2 repeat region (~6:161,030,000-161,060,000) drives much
# of the Lp(a) signal but is poorly resolved in array-based GWAS; SuSiE will
# identify the credible sets of LD-independent signals.
#
# Method (Wallace 2020 + SuSiE-coloc):
#   1. Read 500-kb window from Sinnott Lp(a) exposure + FinnGen CAVS outcome.
#   2. Restrict to common SNPs in 1000G EUR LD panel; compute LD via PLINK.
#   3. Fit SuSiE on exposure z-scores → identifies credible sets.
#   4. Fit SuSiE on outcome z-scores → identifies credible sets.
#   5. Per pair of (exposure CS × outcome CS), compute coloc PP.H4 via
#      coloc::coloc.susie(); summarise by maximum H4 over credible-set pairs.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(coloc)
  library(susieR)
  library(data.table)
})

OUT_CSV <- file.path(paths$results_tables, "susie_coloc_lpa.csv")
WINDOW_CHR <- 6L
WINDOW_LO  <- 160500000L
WINDOW_HI  <- 161100000L
EXP_N      <- 361000L   # Sinnott UKB Lp(a) effective N per registration §7.1
OUT_NCASE  <- 12418L    # FinnGen R12 I9_CAVS_OPERATED
OUT_NCTRL  <- 487930L

# --------------------------------------------------------------------------
# 1. Read 500-kb window for exposure + outcome
# --------------------------------------------------------------------------
exp_path <- file.path(paths$exposures, "SinnottArmstrong_UKB/Lipoprotein_A.imp.gz")
out_path <- file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CAVS_OPERATED.gz")

cat("Reading Sinnott Lp(a) exposure ...\n")
exp_dt <- read_sumstats(exp_path, "SinnottArmstrong2021")
cat(sprintf("  exposure full file: %d rows\n", nrow(exp_dt)))
ensure_chrpos(exp_dt, ancestry = "EUR", label = "exp")
canonicalise_snp_to_rsid(exp_dt, ancestry = "EUR", label = "exp")
exp_win <- exp_dt[!is.na(CHR) & !is.na(BP) & CHR == WINDOW_CHR &
                  BP >= WINDOW_LO & BP <= WINDOW_HI &
                  grepl("^rs[0-9]+$", SNP) & !is.na(BETA) & !is.na(SE) & SE > 0]
exp_win <- unique(exp_win, by = "SNP")
cat(sprintf("  exposure window: %d rs-IDed variants in chr%d:%d-%d\n",
            nrow(exp_win), WINDOW_CHR, WINDOW_LO, WINDOW_HI))

cat("Reading FinnGen R12 CAVS outcome ...\n")
out_dt <- read_sumstats(out_path, "FinnGen_R12")
cat(sprintf("  outcome full file: %d rows\n", nrow(out_dt)))
out_win <- out_dt[!is.na(CHR) & !is.na(BP) & CHR == WINDOW_CHR &
                  BP >= WINDOW_LO & BP <= WINDOW_HI &
                  grepl("^rs[0-9]+$", SNP) & !is.na(BETA) & !is.na(SE) & SE > 0]
out_win <- unique(out_win, by = "SNP")
cat(sprintf("  outcome window: %d rs-IDed variants\n", nrow(out_win)))

# --------------------------------------------------------------------------
# 2. Intersect on SNP, build z-scores
# --------------------------------------------------------------------------
common <- intersect(exp_win$SNP, out_win$SNP)
cat(sprintf("Common SNPs in window: %d\n", length(common)))
if (length(common) < 25) {
  stop("Too few common SNPs (<25) for SuSiE; check window/data.")
}
# Sinnott (GRCh37) and FinnGen R12 (GRCh38) merge by rsID; cross-build
# variant-identifier drift reduces the intersection. The LPA region remains
# the relevant analytical target because canonical Lp(a) lead variants
# (rs10455872, rs3798220, KIV-2 tagging variants) are in dbSNP and present in
# both source GWAS rsID sets.

m <- merge(exp_win[SNP %in% common, .(SNP, BP, EA, NEA, EAF, BETA_e = BETA, SE_e = SE, P_e = as.numeric(P))],
           out_win[SNP %in% common, .(SNP, BETA_o = BETA, SE_o = SE, P_o = as.numeric(P), EA_o = EA, NEA_o = NEA)],
           by = "SNP")

# Align outcome to exposure effect-allele frame
same <- m$EA == m$EA_o & m$NEA == m$NEA_o
flip <- m$EA == m$NEA_o & m$NEA == m$EA_o
m <- m[same | flip]
m[flip, BETA_o := -BETA_o]
m[, c("EA_o", "NEA_o") := NULL]
m[, MAF := pmin(EAF, 1 - EAF)]
m <- m[MAF >= 0.01]
cat(sprintf("After allele alignment + MAF >= 0.01: %d SNPs\n", nrow(m)))

m[, z_e := BETA_e / SE_e]
m[, z_o := BETA_o / SE_o]
setorder(m, BP)

# --------------------------------------------------------------------------
# 3. Build LD matrix via PLINK on the 1000G EUR panel
# --------------------------------------------------------------------------
cat("Building LD matrix via PLINK on 1000G EUR ...\n")
ld_bfile <- file.path(paths$ld_panels, "1000G_EUR", "1000G_EUR")
tmp_dir  <- tempfile("susie_lpa_")
dir.create(tmp_dir)
snp_file <- file.path(tmp_dir, "snps.txt")
fwrite(data.table(SNP = m$SNP), snp_file, col.names = FALSE)

ld_out <- file.path(tmp_dir, "ld")
cmd <- sprintf("plink --bfile %s --extract %s --r square --write-snplist --out %s --silent",
               shQuote(ld_bfile), shQuote(snp_file), shQuote(ld_out))
ret <- system(cmd)
if (ret != 0 || !file.exists(paste0(ld_out, ".ld"))) {
  stop("PLINK LD computation failed (exit ", ret, ").")
}

# Read PLINK's space-separated square LD matrix + the snplist (1 SNP per line)
ld_mat <- as.matrix(fread(paste0(ld_out, ".ld")))
ld_snps <- readLines(paste0(ld_out, ".snplist"))
# Subset m and ld_mat to common SNPs (PLINK may have dropped a few not in panel)
keep_ld <- m$SNP %in% ld_snps
m <- m[keep_ld]
ld_idx <- match(m$SNP, ld_snps)
ld_mat <- ld_mat[ld_idx, ld_idx]
rownames(ld_mat) <- colnames(ld_mat) <- m$SNP   # coloc::check_dataset requires named LD matrix
cat(sprintf("LD matrix: %d × %d\n", nrow(ld_mat), ncol(ld_mat)))

# --------------------------------------------------------------------------
# 4. Build coloc dataset lists for the exposure + outcome
# --------------------------------------------------------------------------
exp_ds <- list(
  beta     = m$BETA_e,
  varbeta  = m$SE_e^2,
  snp      = m$SNP,
  position = m$BP,
  type     = "quant",
  N        = EXP_N,
  MAF      = m$MAF,
  LD       = ld_mat
)
out_ds <- list(
  beta     = m$BETA_o,
  varbeta  = m$SE_o^2,
  snp      = m$SNP,
  position = m$BP,
  type     = "cc",
  s        = OUT_NCASE / (OUT_NCASE + OUT_NCTRL),
  N        = OUT_NCASE + OUT_NCTRL,
  MAF      = m$MAF,
  LD       = ld_mat
)
check_dataset(exp_ds, warn.minp = 5e-8)
check_dataset(out_ds, warn.minp = 5e-8)

# --------------------------------------------------------------------------
# 5. Fit SuSiE on each dataset
# --------------------------------------------------------------------------
cat("Fitting SuSiE on exposure ...\n")
S_exp <- runsusie(exp_ds, suffix = "exp", maxit = 1000)
cat(sprintf("  exposure credible sets: %d\n", length(S_exp$sets$cs)))

cat("Fitting SuSiE on outcome ...\n")
S_out <- runsusie(out_ds, suffix = "out", maxit = 1000)
cat(sprintf("  outcome credible sets: %d\n", length(S_out$sets$cs)))

# --------------------------------------------------------------------------
# 6. SuSiE-coloc: per pair of (exposure CS × outcome CS)
# --------------------------------------------------------------------------
cat("Running coloc.susie ...\n")
susie_res <- tryCatch(
  coloc.susie(S_exp, S_out),
  error = function(e) {
    cat(sprintf("  coloc.susie ERROR: %s\n", conditionMessage(e)))
    NULL
  }
)

# --------------------------------------------------------------------------
# 7. Summarise and write output
# --------------------------------------------------------------------------
out_rows <- list()
if (!is.null(susie_res) && !is.null(susie_res$summary)) {
  cat("\n=== coloc.susie summary ===\n")
  print(susie_res$summary)
  for (i in seq_len(nrow(susie_res$summary))) {
    out_rows[[length(out_rows) + 1]] <- data.table(
      pair = "LPA_CAVS_EUR",
      method = "SuSiE-coloc",
      hit1 = susie_res$summary$hit1[i],
      hit2 = susie_res$summary$hit2[i],
      idx1 = susie_res$summary$idx1[i],
      idx2 = susie_res$summary$idx2[i],
      PP.H0.abf = susie_res$summary$PP.H0.abf[i],
      PP.H1.abf = susie_res$summary$PP.H1.abf[i],
      PP.H2.abf = susie_res$summary$PP.H2.abf[i],
      PP.H3.abf = susie_res$summary$PP.H3.abf[i],
      PP.H4.abf = susie_res$summary$PP.H4.abf[i]
    )
  }
} else {
  cat("\nNo coloc.susie summary (likely 0 credible sets in one dataset).\n")
  out_rows[[1]] <- data.table(
    pair = "LPA_CAVS_EUR",
    method = "SuSiE-coloc",
    hit1 = NA_character_, hit2 = NA_character_,
    idx1 = NA_integer_, idx2 = NA_integer_,
    PP.H0.abf = NA_real_, PP.H1.abf = NA_real_,
    PP.H2.abf = NA_real_, PP.H3.abf = NA_real_, PP.H4.abf = NA_real_
  )
}

result_dt <- rbindlist(out_rows, fill = TRUE)
fwrite(result_dt, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

if (nrow(result_dt) > 0 && any(!is.na(result_dt$PP.H4.abf))) {
  max_h4 <- max(result_dt$PP.H4.abf, na.rm = TRUE)
  cat(sprintf("\nMaximum SuSiE-coloc PP.H4 across credible-set pairs: %.4f\n", max_h4))
  if (max_h4 > 0.5) {
    cat("→ Substantial SuSiE-coloc evidence for shared causal variants at LPA — resolves the coloc.abf single-causal-variant misfit.\n")
  } else if (max_h4 > 0.1) {
    cat("→ Modest SuSiE-coloc evidence; SuSiE identified credible sets but H4 below conventional threshold.\n")
  } else {
    cat("→ Low SuSiE-coloc evidence; LPA architecture may not be amenable to coloc-style sharing even at multi-variant scale, or window/LD insufficient.\n")
  }
}
