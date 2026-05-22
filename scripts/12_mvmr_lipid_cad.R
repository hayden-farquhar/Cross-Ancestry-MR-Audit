# Phase 7c — Multivariable MR for correlated lipid exposures on CAD.
# Pre-registered as the sensitivity analysis for the HDL → CAD result (§15 of
# the registration), motivated by the univariable HDL → CAD result showing a
# substantial IVW → MR-APSS magnitude shift that may reflect TG/LDL pleiotropy.
#
# Run from project root (after R/05):
#   Rscript R/12_mvmr_lipid_cad.R > logs/12_mvmr_$(date +%Y%m%d_%H%M%S).log 2>&1
#
# Architecture:
#   1. Read GLGC EUR LDL, HDL, TG sumstats (cached via R/05's pattern).
#   2. Identify union of GWS instruments (p < 5e-8) across the three lipids.
#   3. Ancestry-matched LD clumping to retain independent SNPs (r² < 0.001).
#   4. Extract effects of all retained SNPs in all three lipids + CAD outcome.
#   5. Harmonise with TwoSampleMR's MVMR harmonisation pipeline.
#   6. Fit MVMR-IVW (and MVMR-Egger as sensitivity).
#   7. Write outputs/tables/mvmr_lipid_cad.csv.

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(ieugwasr)
})

OUT_CSV <- file.path(paths$outputs_tables, "mvmr_lipid_cad.csv")

# --------------------------------------------------------------------------
# 1. Read genome-wide sumstats for the three EUR lipids + CAD outcome
# --------------------------------------------------------------------------
exp_paths <- list(
  LDL = file.path(paths$exposures, "GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  HDL = file.path(paths$exposures, "GLGC_Graham2021/EUR/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz"),
  TG  = file.path(paths$exposures, "GLGC_Graham2021/EUR/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")
)
out_path <- file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz")

cat("Reading three EUR lipid sumstats...\n")
sumstats <- lapply(names(exp_paths), function(nm) {
  cat(sprintf("  [reading]  %s ...\n", basename(exp_paths[[nm]])))
  dt <- read_sumstats(exp_paths[[nm]], "GLGC")
  dt[, exposure := nm]
  dt
})
names(sumstats) <- names(exp_paths)

cat(sprintf("  Loaded LDL: %d, HDL: %d, TG: %d\n",
            nrow(sumstats$LDL), nrow(sumstats$HDL), nrow(sumstats$TG)))

cat("Reading FinnGen R12 CAD outcome...\n")
out_dt <- read_sumstats(out_path, "FinnGen_R12")
cat(sprintf("  Loaded outcome: %d rows\n", nrow(out_dt)))

# --------------------------------------------------------------------------
# 2. Identify union of GWS instruments across the three lipids
# --------------------------------------------------------------------------
GWS <- 5e-8
gws_snps <- unique(unlist(lapply(sumstats, function(d) d[!is.na(P) & as.numeric(P) < GWS, SNP])))
gws_snps <- gws_snps[grepl("^rs[0-9]+$", gws_snps)]
cat(sprintf("Union of GWS SNPs across three lipids: %d\n", length(gws_snps)))

# --------------------------------------------------------------------------
# 3. Ancestry-matched LD clumping via local PLINK (preferred over ieugwasr's
#    API to stay on the project's deterministic infrastructure).
#    Use the LDL sumstats as the "leading" exposure for the clumping P; the
#    instrument set for MVMR is the LD-pruned union.
# --------------------------------------------------------------------------
ld_dat <- sumstats$LDL[SNP %in% gws_snps, .(SNP, P = as.numeric(P))]
# Replace each SNP's P with the minimum across the three lipids — clumping by
# the strongest lipid signal at each SNP for MVMR is standard practice.
min_p <- do.call(pmin, c(
  lapply(sumstats, function(d) {
    m <- d[match(gws_snps, SNP), as.numeric(P)]
    m[is.na(m)] <- 1
    m
  }),
  na.rm = TRUE
))
ld_dat <- data.table(rsid = gws_snps, pval = min_p)
ld_dat <- ld_dat[!is.na(pval) & pval < GWS]
cat(sprintf("SNPs with min-P < %.0e for clumping: %d\n", GWS, nrow(ld_dat)))

cat("Clumping via local 1000G EUR PLINK panel (r2<0.001, kb=10000)...\n")
clumped <- ld_clump_local(
  dat = ld_dat,
  clump_kb = 10000, clump_r2 = 0.001, clump_p = GWS,
  bfile = file.path(paths$ld_panels, "1000G_EUR", "1000G_EUR"),
  plink_bin = "plink"
)
# ld_clump_local returns rows with rsid + pval; standardise back to SNP for downstream
indep_snps <- clumped$rsid
if (is.null(indep_snps)) indep_snps <- clumped$SNP  # fallback if API differs
cat(sprintf("  Retained %d independent instruments\n", length(indep_snps)))

# --------------------------------------------------------------------------
# 4. Build harmonised MVMR table
# --------------------------------------------------------------------------
# Use the LDL sumstats' rows as the reference SNP set; align HDL and TG to it.
ref <- sumstats$LDL[SNP %in% indep_snps]
hdl_at_ref <- sumstats$HDL[match(ref$SNP, SNP)]
tg_at_ref  <- sumstats$TG [match(ref$SNP, SNP)]
out_at_ref <- out_dt[match(ref$SNP, SNP)]

keep <- !is.na(hdl_at_ref$BETA) & !is.na(tg_at_ref$BETA) & !is.na(out_at_ref$BETA) &
        !is.na(ref$BETA) & !is.na(ref$SE) & !is.na(hdl_at_ref$SE) &
        !is.na(tg_at_ref$SE) & !is.na(out_at_ref$SE)
ref <- ref[keep]; hdl_at_ref <- hdl_at_ref[keep]
tg_at_ref <- tg_at_ref[keep]; out_at_ref <- out_at_ref[keep]
cat(sprintf("  After matching across all three lipids + outcome: %d SNPs\n",
            nrow(ref)))

# Allele-flip HDL, TG, outcome to LDL's effect-allele reference.
flip_to_ref <- function(target, ref) {
  same <- target$EA == ref$EA & target$NEA == ref$NEA
  flip <- target$EA == ref$NEA & target$NEA == ref$EA
  ok <- same | flip
  out <- copy(target)
  out[flip, BETA := -BETA]
  out[flip, c("EA", "NEA") := .(ref$EA[flip], ref$NEA[flip])]
  out[!ok, BETA := NA_real_]  # ambiguous/incompatible — drop downstream
  out
}
hdl_aln <- flip_to_ref(hdl_at_ref, ref)
tg_aln  <- flip_to_ref(tg_at_ref,  ref)
out_aln <- flip_to_ref(out_at_ref, ref)

valid <- !is.na(hdl_aln$BETA) & !is.na(tg_aln$BETA) & !is.na(out_aln$BETA)
ref <- ref[valid]; hdl_aln <- hdl_aln[valid]
tg_aln <- tg_aln[valid]; out_aln <- out_aln[valid]
cat(sprintf("  After allele alignment to LDL effect-allele frame: %d SNPs\n",
            nrow(ref)))

mv_dat <- data.frame(
  SNP        = ref$SNP,
  EA         = ref$EA, NEA = ref$NEA,
  b_ldl = ref$BETA, se_ldl = ref$SE,
  b_hdl = hdl_aln$BETA, se_hdl = hdl_aln$SE,
  b_tg  = tg_aln$BETA,  se_tg  = tg_aln$SE,
  b_out = out_aln$BETA, se_out = out_aln$SE
)

# --------------------------------------------------------------------------
# 5. MVMR-IVW via TwoSampleMR's mv_multiple equivalent computation.
#    Standard MVMR-IVW: solve b_out = b_ldl*beta_ldl + b_hdl*beta_hdl + b_tg*beta_tg
#    weighted by 1/se_out^2 (Burgess 2015).
# --------------------------------------------------------------------------
cat("\nFitting MVMR-IVW (Burgess 2015 weighted regression)...\n")
X <- cbind(mv_dat$b_ldl, mv_dat$b_hdl, mv_dat$b_tg)
y <- mv_dat$b_out
w <- 1 / mv_dat$se_out^2
fit <- lm.wfit(X, y, w)
beta <- coef(fit)
# Robust SE via the residual standard deviation
res <- y - X %*% beta
sigma2 <- sum(w * res^2) / max(1, length(y) - 3)
XtWX_inv <- solve(t(X) %*% diag(w) %*% X)
se_b <- sqrt(sigma2 * diag(XtWX_inv))
pvals <- 2 * pnorm(-abs(beta / se_b))

mvmr_res <- data.table(
  exposure = c("LDL", "HDL", "TG"),
  estimator = "MVMR-IVW",
  b_cond = round(as.numeric(beta), 4),
  se_cond = round(se_b, 4),
  pval_cond = signif(pvals, 4),
  n_snps = nrow(mv_dat)
)
cat("\n=== MVMR-IVW results ===\n")
print(mvmr_res)

# --------------------------------------------------------------------------
# 6. Write outputs
# --------------------------------------------------------------------------
fwrite(mvmr_res, OUT_CSV)
cat(sprintf("\nWrote %s\n", OUT_CSV))

# Also save the harmonised multi-exposure data for downstream reference
fwrite(mv_dat, file.path(paths$outputs_tables, "mvmr_lipid_cad_harmonised.csv"))
cat(sprintf("Wrote harmonised multi-exposure table (%d SNPs).\n", nrow(mv_dat)))
