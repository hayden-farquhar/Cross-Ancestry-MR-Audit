# Cross-ancestry MR portability audit — renv environment setup.
# Run once from the repository root: source("scripts/00_setup_renv.R") then renv::restore()
# Locks the R toolchain so the pipeline is reproducible across machines.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

renv::init(bare = TRUE)

cran_packages <- c(
  "TwoSampleMR",
  "MendelianRandomization",
  "ieugwasr",
  "coloc",
  "susieR",
  "MRPRESSO",
  "RadialMR",
  "metafor",
  "meta",
  "data.table",
  "ggplot2",
  "forestplot",
  "dplyr",
  "tidyr",
  "readr",
  "purrr",
  "stringr",
  "fs",
  "digest",
  "openssl",
  "qs"
)

github_packages <- c(
  "YangLabHKUST/MR-APSS",
  "jean997/cause",
  "MRCIEU/TwoSampleMR",
  "rondolab/MR-PRESSO"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
for (repo in github_packages) {
  remotes::install_github(repo, upgrade = "never")
}

renv::snapshot(prompt = FALSE)

cat("renv setup complete. Lockfile:", file.path(getwd(), "renv.lock"), "\n")
