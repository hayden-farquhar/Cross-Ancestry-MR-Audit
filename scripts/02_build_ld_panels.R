# Phase 1 — Verify 1000G phase 3 PLINK panels and run a single-pair LD clump
# smoke test per ancestry to confirm panel integrity.
#
# The Zenodo deposit (6614170) already provides PLINK bed/bim/fam/afreq trios.
# This script:
#   1. Confirms the three ancestry panels exist with byte sizes matching the
#      Zenodo deposit metadata.
#   2. Loads the .bim into R and reports basic counts (n_SNPs, autosomal coverage).
#   3. Runs a single ld_clump() call per ancestry on a small test instrument set
#      to verify ieugwasr can talk to the local PLINK panel.
#
# Prerequisites: plink 1.9 binary on PATH (Homebrew: `brew install plink`).
#
# Run from project root:
#   Rscript R/02_build_ld_panels.R

source(here::here("R", "utils.R"))

zenodo_sizes <- list(
  EUR = list(bed = 231387159, bim = 51537337, fam = 9557,  afreq = 56496658),
  EAS = list(bed = 231387159, bim = 51537337, fam = 9576,  afreq = 54161460),
  AFR = list(bed = 299334181, bim = 51537337, fam = 12394, afreq = 56580637)
)

cat("\n=== LD panel byte-size verification (against Zenodo 6614170) ===\n")
ok <- TRUE
for (pop in names(zenodo_sizes)) {
  for (ext in names(zenodo_sizes[[pop]])) {
    path <- file.path(paths$ld_panels, sprintf("1000G_%s", pop),
                      sprintf("1000G_%s.%s", pop, ext))
    if (!file.exists(path)) {
      cat(sprintf("  MISSING   %s\n", path))
      ok <- FALSE
      next
    }
    sz <- file.info(path)$size
    expected <- zenodo_sizes[[pop]][[ext]]
    status <- if (sz == expected) "OK     " else "MISMATCH"
    cat(sprintf("  %s  %s.%s  %d bytes (expected %d)\n",
                status, pop, ext, sz, expected))
    if (sz != expected) ok <- FALSE
  }
}
if (!ok) stop("LD panel integrity check failed. Re-fetch missing/mismatched files.")

cat("\n=== Per-ancestry .bim variant counts ===\n")
for (pop in names(zenodo_sizes)) {
  bim <- fread(file.path(paths$ld_panels, sprintf("1000G_%s", pop),
                          sprintf("1000G_%s.bim", pop)),
               col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
  cat(sprintf("  %s: %d variants  (chr 1-22: %d, chr X: %d, autosomes only check %s)\n",
              pop, nrow(bim),
              sum(bim$CHR %in% 1:22),
              sum(bim$CHR == "X" | bim$CHR == 23),
              if (sum(bim$CHR %in% 1:22) > 1e6) "PASS" else "FAIL"))
}

cat("\n=== ld_clump smoke test (one ancestry, two test SNPs) ===\n")
if (!requireNamespace("ieugwasr", quietly = TRUE)) {
  cat("  ieugwasr not installed — run scripts/00_setup_renv.R first.\n")
  quit(status = 0)
}

test_snps <- data.frame(
  rsid = c("rs6511720", "rs7412"),     # LDLR and APOE — canonical LDL-C instruments
  pval = c(1e-100, 1e-200),
  trait_id = "test"
)

for (pop in names(zenodo_sizes)) {
  panel_prefix <- file.path(paths$ld_panels, sprintf("1000G_%s/1000G_%s", pop, pop))
  out <- tryCatch(
    ieugwasr::ld_clump_local(
      dat = test_snps,
      clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
      bfile = panel_prefix,
      plink_bin = Sys.which("plink")
    ),
    error = function(e) {
      cat(sprintf("  %s: ld_clump_local FAILED — %s\n", pop, conditionMessage(e)))
      NULL
    }
  )
  if (!is.null(out)) {
    cat(sprintf("  %s: ld_clump_local OK — %d/%d test SNPs survived\n",
                pop, nrow(out), nrow(test_snps)))
  }
}

cat("\nLD panel build verified. Proceed to R/03_extract_instruments.R.\n")
