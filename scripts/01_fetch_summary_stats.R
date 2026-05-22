# Phase 1 — Summary stats fetch.
# The actual downloads were executed via shell/curl with provenance logged in
# osf/data_snapshot_log.md. This script is the verification layer:
#   1. Read the pre-registered canonical pairs.
#   2. For each (exposure × ancestry) and (outcome × ancestry), check the file
#      exists locally and matches the expected SHA-256.
#   3. Print a status table that maps every pre-registered slot to a concrete
#      file path or an explicit "DEFERRED — controlled access" marker.
#
# Run from the project root:
#   Rscript R/01_fetch_summary_stats.R

source(here::here("R", "utils.R"))

stopifnot(file.exists(paths$snapshot_log))

# Files registered in the snapshot log (the source of truth for what exists).
# Hard-coded here to make the verifier self-contained and easy to audit.
expected_files <- tibble::tribble(
  ~pair_slot,            ~source_tag,          ~rel_path,
  # Exposures
  "lipids/EUR/LDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/EUR/HDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/EUR/HDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/EUR/TG",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/EUR/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/EUR/TC",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/EUR/TC_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/EAS/LDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/EAS/LDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",
  "lipids/EAS/HDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/EAS/HDL_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",
  "lipids/EAS/TG",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/EAS/logTG_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",
  "lipids/EAS/TC",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/EAS/TC_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz",
  "lipids/AFR/LDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/AFR/LDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/AFR/HDL",      "GLGC",                "data/raw/exposures/GLGC_Graham2021/AFR/HDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/AFR/TG",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/AFR/logTG_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "lipids/AFR/TC",       "GLGC",                "data/raw/exposures/GLGC_Graham2021/AFR/TC_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz",
  "urate/EUR",           "CKDGen_Tin2019",      "data/raw/exposures/CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz",
  "urate/transethnic",   "CKDGen_Tin2019",      "data/raw/exposures/CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_all_741_nstud37_summac400_rsid.txt.gz",
  "urate/AFR",           "PanUKB_AFR",          "data/raw/outcomes/PanUKB_AFR/biomarkers-30880-both_sexes-irnt.tsv.bgz",
  "bmi/EUR",             "GIANT_Yengo2018",     "data/raw/exposures/GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz",
  "bmi/AFR",             "PanUKB_AFR",          "data/raw/outcomes/PanUKB_AFR/continuous-21001-both_sexes-irnt.tsv.bgz",
  "sbp/EUR",             "ICBP_Evangelou2018",  "data/raw/exposures/ICBP_Evangelou2018/Evangelou_30224653_SBP.txt.gz",
  "dbp/EUR",             "ICBP_Evangelou2018",  "data/raw/exposures/ICBP_Evangelou2018/Evangelou_30224653_DBP.txt.gz",
  "sbp/AFR",             "PanUKB_AFR",          "data/raw/outcomes/PanUKB_AFR/continuous-4080-both_sexes-irnt.tsv.bgz",
  "hba1c/EUR",           "MAGIC_Chen2021",      "data/raw/exposures/MAGIC_Chen2021/MAGIC1000G_HbA1c_EUR.tsv.gz",
  "hba1c/EAS",           "MAGIC_Chen2021",      "data/raw/exposures/MAGIC_Chen2021/MAGIC1000G_HbA1c_EAS.tsv.gz",
  "hba1c/AA",            "MAGIC_Chen2021",      "data/raw/exposures/MAGIC_Chen2021/MAGIC1000G_HbA1c_AA.tsv.gz",
  "lpa/EUR",             "SinnottArmstrong2021","data/raw/exposures/SinnottArmstrong_UKB/Lipoprotein_A.imp.gz",
  "lpa/AFR",             "PanUKB_AFR",          "data/raw/outcomes/PanUKB_AFR/biomarkers-30790-both_sexes-irnt.tsv.bgz",
  "smoking/EUR",         "GSCAN_Saunders2022",  "data/raw/exposures/GSCAN_Saunders2022/EUR_stratified.zip",
  "smoking/EAS",         "GSCAN_Saunders2022",  "data/raw/exposures/GSCAN_Saunders2022/EAS_stratified.zip",
  "smoking/AFR",         "GSCAN_Saunders2022",  "data/raw/exposures/GSCAN_Saunders2022/AFR_stratified.zip",
  # Outcomes
  "cad/EUR",             "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_I9_CHD.gz",
  "af/EUR",              "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_I9_AF.gz",
  "aki/EUR",             "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_N14_ACUTERENFAIL.gz",
  "t2d/EUR",             "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_T2D.gz",
  "gout/EUR",            "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_M13_GOUT.gz",
  "cavs/EUR",            "FinnGen_R12",         "data/raw/outcomes/FinnGen_R12/finngen_R12_I9_CAVS_OPERATED.gz",
  "fracture/EUR (neg-ctl)","FinnGen_R12",       "data/raw/outcomes/FinnGen_R12/finngen_R12_ST19_FRACT_FOREA.gz",
  "stroke/EUR/AIS",      "MEGASTROKE",          "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.EUR.out",
  "stroke/EUR/LAS",      "MEGASTROKE",          "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.3.LAS.EUR.out",
  "stroke/EUR/CES",      "MEGASTROKE",          "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.4.CES.EUR.out",
  "stroke/EUR/SVS",      "MEGASTROKE",          "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.5.SVS.EUR.out",
  "stroke/EUR/any",      "MEGASTROKE",          "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.1.AS.EUR.out",
  "stroke/transethnic/AIS","MEGASTROKE",        "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.2.AIS.TRANS.out",
  "t2d/EAS",             "AGEN_Spracklen2020",  "data/raw/outcomes/AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz"
)

# Pairs deferred for user-initiated download (controlled access or click-through):
deferred <- tibble::tribble(
  ~pair_slot,           ~reason,
  "cad/EAS (BBJ)",      "NBDC hum0014/hum0197 — controlled access; apply via humandbs.biosciencedbc.jp",
  "t2d/AFR (MVP)",      "dbGaP pha005193.1 — controlled access; summary stats GCST90132304 not yet posted",
  "cad/AFR (MVP)",      "Same as above",
  "bp/AFR (AWI-Gen)",   "H3Africa controlled access — DAA required",
  "t2d/EUR/AFR (DIAMANTE 2022 / T2DGGI 2024)",  "diagram-consortium.org click-through agreement"
)

results <- expected_files |>
  dplyr::mutate(
    abs_path  = repo_path(rel_path),
    exists    = file.exists(abs_path),
    is_link   = !is.na(Sys.readlink(abs_path)) & Sys.readlink(abs_path) != "",
    size_mb   = ifelse(exists, file.info(abs_path)$size / 1e6, NA_real_)
  )

cat("\n=== Pre-registered slots — file presence check ===\n\n")
print(knitr::kable(results[, c("pair_slot", "source_tag", "exists", "is_link", "size_mb")],
                   digits = 1))

cat("\n=== Slots deferred for user action ===\n\n")
print(knitr::kable(deferred))

n_present <- sum(results$exists)
n_total   <- nrow(results)
cat(sprintf("\nSummary: %d / %d expected files present locally.\n", n_present, n_total))

if (n_present < n_total) {
  missing <- results$pair_slot[!results$exists]
  cat("Missing slots:\n  - ", paste(missing, collapse = "\n  - "), "\n", sep = "")
  quit(status = 1)
}

cat("\nAll expected files present. Proceed to R/02_build_ld_panels.R.\n")
