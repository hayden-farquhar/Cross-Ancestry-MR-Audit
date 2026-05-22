# Smoke tests for R/utils.R::read_sumstats() parsers.
# Verifies each source file parses to the expected 10-column schema with
# non-empty BETA, SE, and P columns. Catches schema drift before any pipeline
# stage that consumes downstream.
#
# Run from project root:
#   Rscript tests/test_parsers.R

source(here::here("R", "utils.R"))

cases <- list(
  list(tag = "FinnGen_R12",
       path = file.path(paths$outcomes, "FinnGen_R12/finngen_R12_I9_CHD.gz")),
  list(tag = "GLGC",
       path = file.path(paths$exposures, "GLGC_Graham2021/EAS/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz")),
  list(tag = "MEGASTROKE",
       path = file.path(paths$outcomes, "MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out")),
  list(tag = "AGEN_T2D",
       path = file.path(paths$outcomes, "AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz")),
  list(tag = "SinnottArmstrong2021",
       path = file.path(paths$exposures, "SinnottArmstrong_UKB/Lipoprotein_A.imp.gz")),
  list(tag = "ICBP_Evangelou2018",
       path = file.path(paths$exposures, "ICBP_Evangelou2018/Evangelou_30224653_SBP.txt.gz")),
  list(tag = "MAGIC_Chen2021",
       path = file.path(paths$exposures, "MAGIC_Chen2021/MAGIC1000G_HbA1c_EUR.tsv.gz")),
  list(tag = "CKDGen_Tin2019",
       path = file.path(paths$exposures, "CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz")),
  list(tag = "GIANT_Yengo2018",
       path = file.path(paths$exposures, "GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz")),
  list(tag = "PanUKB_AFR",
       path = file.path(paths$outcomes, "PanUKB_AFR/biomarkers-30790-both_sexes-irnt.tsv.bgz"))
)

expected_cols <- c("SNP", "CHR", "BP", "EA", "NEA", "EAF", "BETA", "SE", "P", "N")
fails <- 0
for (c in cases) {
  if (!file.exists(c$path)) {
    cat(sprintf("  SKIP %-22s — file not present: %s\n", c$tag, c$path))
    next
  }
  res <- tryCatch({
    dt <- read_sumstats(c$path, c$tag)
    stopifnot(all(expected_cols %in% names(dt)))
    n_nonempty <- sum(!is.na(dt$BETA) & !is.na(dt$SE) & !is.na(dt$P))
    if (n_nonempty < 100) stop(sprintf("only %d non-empty SNP rows", n_nonempty))
    "PASS"
  }, error = function(e) sprintf("FAIL — %s", conditionMessage(e)))
  if (substr(res, 1, 4) == "FAIL") fails <- fails + 1
  cat(sprintf("  %s  %-22s\n", res, c$tag))
}

if (fails > 0) {
  cat(sprintf("\n%d parser test(s) failed. Fix R/utils.R before running pipeline.\n", fails))
  quit(status = 1)
}
cat("\nAll parsers PASS.\n")
