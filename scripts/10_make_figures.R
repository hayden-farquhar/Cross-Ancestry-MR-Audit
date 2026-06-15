# Phase 7 helper — generate forest, funnel, and heterogeneity-comparison plots
# from outputs/tables/*.csv. Writes PNGs to outputs/figures/.
#
# Run from project root:
#   Rscript R/10_make_figures.R

source(here::here("R", "utils.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(data.table)
  library(forestplot)
})

# Per-pair × ancestry forest plot across MR methods.
# Filters to causal-estimate rows (excludes Q heterogeneity, MR-PRESSO Global,
# and MR-Egger intercept — those report different quantities than per-SD log-OR
# and would compress the OR axis if rendered alongside).
.CAUSAL_METHODS <- c(
  "MR-APSS",  # primary estimator — registration §14.1
  "Inverse variance weighted", "MR Egger",
  "Weighted median", "Weighted mode",
  "MR-PRESSO (raw)", "MR-PRESSO (outlier-corrected)",
  "Radial-IVW"
)

draw_pair_forest <- function(pair_name, mr_dt) {
  # Exact ancestry-suffix match (NOT a prefix grepl): a prefix match on e.g.
  # "URATE_GOUT" would also pull in the sensitivity variants URATE_GOUT_TRANS_EUR
  # and URATE_GOUT_EAS2, contaminating the canonical EUR vs EAS forest.
  d <- mr_dt[pair %in% paste0(pair_name, "_", c("EUR", "EAS", "AFR")) &
               method %in% .CAUSAL_METHODS]
  if (nrow(d) == 0) return(invisible(NULL))
  d <- d[!is.na(b) & !is.na(se)]
  if (nrow(d) == 0) return(invisible(NULL))
  d[, ancestry := sub(".*_(EUR|EAS|AFR)$", "\\1", pair)]
  d[, OR  := exp(b)]
  d[, LCI := exp(b - 1.96 * se)]
  d[, UCI := exp(b + 1.96 * se)]

  # Auto-scale x axis: bound at the 5th and 95th percentile of LCI/UCI
  xmin <- min(d$LCI, na.rm = TRUE)
  xmax <- max(d$UCI, na.rm = TRUE)
  xmin <- max(xmin, 0.1)
  xmax <- min(xmax, 10)

  # Friendly axis label: pretty exposure name from the pair label
  pretty_exposure <- sub("_[A-Z]+$", "",
                         sub("^[A-Z0-9]+_", "", pair_name))
  if (pretty_exposure == "") pretty_exposure <- pair_name
  axis_label <- "OR per 1-SD exposure (binary outcome)"

  n_ancestries <- length(unique(d$ancestry))

  p <- ggplot(d, aes(x = OR, y = method, xmin = LCI, xmax = UCI,
                     colour = ancestry)) +
    geom_pointrange() +
    geom_vline(xintercept = 1, linetype = "dashed") +
    facet_wrap(~ ancestry, ncol = 1, scales = "free_y") +
    scale_x_log10(limits = c(xmin * 0.9, xmax * 1.1)) +
    labs(x = axis_label, y = "MR method",
         title = sprintf("%s — per-ancestry × per-method", pair_name),
         subtitle = "Causal estimators only (Q, Global pleiotropy, Egger intercept excluded)") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, colour = "grey40"))

  out_path <- file.path(paths$outputs_figs, sprintf("forest_%s.png", pair_name))
  ggsave(out_path, p, width = 7,
         height = 2 + 2.2 * n_ancestries, dpi = 150)
  invisible(out_path)
}

# Cross-pair heterogeneity summary plot
draw_heterogeneity_summary <- function(het_dt) {
  p <- ggplot(het_dt, aes(x = reorder(pair, I2_pct), y = I2_pct,
                          fill = sig_at_Bonferroni_05)) +
    geom_col() +
    geom_hline(yintercept = 50, linetype = "dashed", colour = "grey") +
    geom_hline(yintercept = 75, linetype = "dashed", colour = "grey") +
    coord_flip() +
    scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
                      name = "Bonferroni p<0.05") +
    labs(x = NULL, y = "Cross-ancestry I² (%)",
         title = "Cross-ancestry heterogeneity per canonical pair") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom")

  out_path <- file.path(paths$outputs_figs, "heterogeneity_summary.png")
  ggsave(out_path, p, width = 8, height = 5, dpi = 150)
  invisible(out_path)
}

# --- Run plots ---
mr_path  <- file.path(paths$outputs_tables, "mr_estimates_per_method.csv")
het_path <- file.path(paths$outputs_tables, "heterogeneity_meta.csv")

if (file.exists(mr_path)) {
  mr_dt <- fread(mr_path)
  # Identify unique pair_root values (strip trailing ancestry suffix)
  pair_roots <- unique(sub("_(EUR|EAS|AFR)$", "", mr_dt$pair))
  pair_roots <- pair_roots[!grepl("_TRANS$|EAS2$", pair_roots)]
  for (pr in pair_roots) draw_pair_forest(pr, mr_dt)
  cat(sprintf("Wrote %d forest plots to %s\n",
              length(pair_roots), paths$outputs_figs))
} else {
  cat("MR estimates not yet computed.\n")
}

if (file.exists(het_path)) {
  het_dt <- fread(het_path)
  draw_heterogeneity_summary(het_dt)
  cat("Wrote heterogeneity_summary.png\n")
}
