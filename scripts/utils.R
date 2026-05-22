# Project utility helpers, sourced by every numbered analysis script.
# Keep this thin: paths, SHA-256 helpers, manifest readers. No analysis here.

suppressPackageStartupMessages({
  library(data.table)
  library(here)
})

repo_path <- function(...) here::here(...)

paths <- list(
  exposures      = repo_path("data", "raw", "exposures"),
  outcomes       = repo_path("data", "raw", "outcomes"),
  ld_panels      = repo_path("data", "ld_panels"),
  processed      = repo_path("data", "processed"),
  reference      = repo_path("data", "reference"),
  outputs_tables = repo_path("outputs", "tables"),
  outputs_figs   = repo_path("outputs", "figures"),
  osf            = repo_path("osf"),
  snapshot_log   = repo_path("osf", "data_snapshot_log.md")
)

read_canonical_pairs <- function() {
  fread(repo_path("outputs", "tables", "canonical_pairs_prespec.csv"))
}

sha256_file <- function(path) {
  # Returns lowercase hex SHA-256 of a file. Streaming; does not load file into memory.
  if (!file.exists(path)) stop("File not found: ", path)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

# Schema-aware reader. Each external source has a unique column convention;
# parsers below verified against actual file headers as of 2026-05-20.
# Returns a data.table with standardised columns: SNP, CHR, BP, EA, NEA, EAF, BETA, SE, P, N
# Caller supplies the source-type tag.
read_sumstats <- function(path, source_type) {
  if (!file.exists(path)) stop("Sumstats file not found: ", path)
  switch(source_type,
    "FinnGen_R12" = {
      raw <- fread(path)
      data.table(
        SNP  = raw$rsids, CHR = raw[[1]], BP = raw$pos,
        EA   = raw$alt, NEA = raw$ref, EAF = raw$af_alt,
        BETA = raw$beta, SE = raw$sebeta, P = raw$pval, N = NA_integer_
      )
    },
    "GLGC" = {
      raw <- fread(path)
      # GLGC 2021 ships in two schemas:
      #   (a) U Mich CSG portal — columns: rsID, CHROM, POS_b37, REF, ALT, N,
      #       N_studies, POOLED_ALT_AF, EFFECT_SIZE, SE, pvalue, pvalue_GC, etc.
      #   (b) EBI GWAS Catalog harmonised — columns: variant_id, chromosome,
      #       base_pair_location, other_allele, effect_allele, n, N_studies,
      #       effect_allele_frequency, beta, standard_error, p_value.
      # Detect via column-name presence and remap.
      if ("rsID" %in% names(raw)) {
        data.table(
          SNP  = raw$rsID, CHR = raw$CHROM, BP = raw$POS_b37,
          EA   = raw$ALT, NEA = raw$REF, EAF = raw$POOLED_ALT_AF,
          BETA = raw$EFFECT_SIZE, SE = raw$SE,
          P    = raw$pvalue_GC,  # genomic-control-corrected raw p
          N    = raw$N
        )
      } else if ("variant_id" %in% names(raw)) {
        data.table(
          SNP  = raw$variant_id, CHR = raw$chromosome, BP = raw$base_pair_location,
          EA   = toupper(raw$effect_allele), NEA = toupper(raw$other_allele),
          EAF  = suppressWarnings(as.numeric(raw$effect_allele_frequency)),
          BETA = suppressWarnings(as.numeric(raw$beta)),
          SE   = suppressWarnings(as.numeric(raw$standard_error)),
          P    = suppressWarnings(as.numeric(raw$p_value)),
          N    = suppressWarnings(as.integer(raw$n))
        )
      } else {
        stop("Unknown GLGC schema. Header columns: ",
             paste(names(raw), collapse = ", "))
      }
    },
    "MEGASTROKE" = {
      raw <- fread(path)  # space-separated; data.table auto-detects
      data.table(
        SNP  = raw$MarkerName, CHR = NA_integer_, BP = NA_integer_,
        EA   = toupper(raw$Allele1), NEA = toupper(raw$Allele2), EAF = raw$Freq1,
        BETA = raw$Effect, SE = raw$StdErr,
        P    = if ("P-value" %in% names(raw)) raw[["P-value"]] else raw$P,
        N    = if ("TotalSampleSize" %in% names(raw)) raw$TotalSampleSize else NA_integer_
      )
    },
    "AGEN_T2D" = {
      raw <- fread(path)
      data.table(
        SNP  = raw$MarkerName, CHR = raw$Chr, BP = raw$Pos,
        EA   = toupper(raw$EA), NEA = toupper(raw$NEA), EAF = raw$EAF,
        BETA = raw$Beta, SE = raw$SE, P = raw$P, N = raw$Neff
      )
    },
    "SinnottArmstrong2021" = {
      raw <- fread(path)
      data.table(
        SNP  = raw$MarkerName, CHR = raw[["#CHROM"]], BP = raw$POS,
        EA   = toupper(raw$ALT), NEA = toupper(raw$REF), EAF = raw$MAF,
        BETA = raw$Effect, SE = raw$StdErr,
        P    = raw[["P-value"]], N = NA_integer_
      )
    },
    "ICBP_Evangelou2018" = {
      # ICBP MarkerName format is `chr:pos:SNP`; extract CHR and BP so that the
      # downstream chr:pos lookup against 1000G can resolve rsIDs.
      raw <- fread(path)  # space-separated
      mn_parts <- tstrsplit(raw$MarkerName, ":", fixed = TRUE)
      data.table(
        SNP  = raw$MarkerName,
        CHR  = suppressWarnings(as.integer(mn_parts[[1]])),
        BP   = suppressWarnings(as.integer(mn_parts[[2]])),
        EA   = toupper(raw$Allele1), NEA = toupper(raw$Allele2),
        EAF  = suppressWarnings(as.numeric(raw$Freq1)),
        BETA = suppressWarnings(as.numeric(raw$Effect)),
        SE   = suppressWarnings(as.numeric(raw$StdErr)),
        P    = suppressWarnings(as.numeric(raw$P)),
        N    = suppressWarnings(as.integer(raw$TotalSampleSize))
      )
    },
    "MAGIC_Chen2021" = {
      raw <- fread(path)
      data.table(
        SNP  = raw$variant, CHR = raw$chromosome, BP = raw$base_pair_location,
        EA   = toupper(raw$effect_allele), NEA = toupper(raw$other_allele),
        EAF  = raw$effect_allele_frequency,
        BETA = raw$beta, SE = raw$standard_error, P = raw$p_value,
        N    = raw$sample_size
      )
    },
    "CKDGen_Tin2019" = {
      # CKDGen files are space-separated and some rows have an empty RSID column
      # (variants without an assigned rs number); use fill=TRUE so fread tolerates
      # the inconsistent column count.
      raw <- fread(path, fill = TRUE)
      data.table(
        SNP  = raw$RSID, CHR = raw$Chr, BP = raw$Pos_b37,
        EA   = toupper(raw$Allele1), NEA = toupper(raw$Allele2),
        EAF  = suppressWarnings(as.numeric(raw$Freq1)),
        BETA = suppressWarnings(as.numeric(raw$Effect)),
        SE   = suppressWarnings(as.numeric(raw$StdErr)),
        P    = if ("P-value" %in% names(raw))
                 suppressWarnings(as.numeric(raw[["P-value"]]))
               else suppressWarnings(as.numeric(raw$P)),
        N    = suppressWarnings(as.integer(raw$n_total_sum))
      )[!is.na(SNP) & SNP != ""]  # drop rows where RSID was missing
    },
    "PanUKB_AFR" = {
      # Pan-UKB per-phenotype files are multi-ancestry — extract AFR columns.
      # No rsid column; construct chr:pos:ref:alt as SNP if needed for matching.
      # File is bgzipped — data.table::fread handles via gzip pipe.
      raw <- fread(cmd = sprintf("gzip -dc '%s'", path))
      data.table(
        SNP  = paste(raw$chr, raw$pos, raw$ref, raw$alt, sep = ":"),
        CHR  = raw$chr,
        BP   = raw$pos,
        EA   = toupper(raw$alt), NEA = toupper(raw$ref),
        EAF  = raw$af_AFR,
        BETA = raw$beta_AFR, SE = raw$se_AFR,
        P    = 10^(-raw$neglog10_pval_AFR),
        N    = NA_integer_
      )
    },
    "GIANT_Yengo2018" = {
      raw <- fread(path)
      # Note: GIANT Yengo uses Freq_Tested_Allele_in_HRS (HRS reference cohort);
      # not a sample-EAF but suitable as a strand-resolving proxy.
      data.table(
        SNP  = raw$SNP, CHR = raw$CHR, BP = raw$POS,
        EA   = toupper(raw$Tested_Allele), NEA = toupper(raw$Other_Allele),
        EAF  = raw$Freq_Tested_Allele_in_HRS,
        BETA = raw$BETA, SE = raw$SE, P = raw$P, N = raw$N
      )
    },
    stop("Unknown source_type: ", source_type,
         ". Add a parser in R/utils.R::read_sumstats().")
  )
}

# Lookup table: chr:pos → rsID, built once from the ancestry panel BIM.
# Used for sources whose marker column is not rs-ID (ICBP, Sinnott-Armstrong).
# Cached per ancestry to avoid re-reading the BIM repeatedly.
.bim_cache <- new.env(parent = emptyenv())

panel_chrpos_to_rsid <- function(ancestry) {
  key <- toupper(ancestry)
  if (!is.null(.bim_cache[[key]])) return(.bim_cache[[key]])
  bim_path <- file.path(
    paths$ld_panels, sprintf("1000G_%s", key),
    sprintf("1000G_%s.bim", key)
  )
  if (!file.exists(bim_path)) {
    stop("BIM file not found for ancestry ", ancestry, ": ", bim_path)
  }
  bim <- fread(
    bim_path, header = FALSE,
    col.names = c("CHR", "SNP", "CM", "BP", "A1", "A2")
  )
  bim[, key := paste0(CHR, ":", BP)]
  .bim_cache[[key]] <- bim
  invisible(bim)
}

# Replace non-rsID instrument SNP IDs with rsIDs from the ancestry-matched
# 1000G BIM (chr:pos match). Drops SNPs without a panel match.
resolve_snp_ids_to_panel <- function(dat, ancestry) {
  stopifnot(all(c("SNP", "CHR", "BP") %in% names(dat)))
  needs_lookup <- !grepl("^rs[0-9]+$", dat$SNP)
  if (!any(needs_lookup)) return(dat)
  cat(sprintf("  Resolving %d non-rsID markers against 1000G %s BIM ...\n",
              sum(needs_lookup), ancestry))
  bim <- panel_chrpos_to_rsid(ancestry)
  dat[, key := paste0(CHR, ":", BP)]
  m <- match(dat$key, bim$key)
  resolved <- !is.na(m)
  dat[resolved, SNP := bim$SNP[m[resolved]]]
  dat <- dat[grepl("^rs", SNP)]  # drop SNPs that still lack an rsID
  dat[, key := NULL]
  cat(sprintf("  Retained %d variants with panel-matched rsIDs\n", nrow(dat)))
  return(dat)
}

# Populate missing CHR/BP via the ancestry-matched BIM (rsID → chr:pos).
# Used for sources whose SNP column carries rsIDs but no chr/bp columns
# (notably MEGASTROKE). No-op if ≥ 50 % of rows already have CHR.
# Mirrors R/07's helper; lifted here so R/05 (and any future caller) can use it.
ensure_chrpos <- function(dt, ancestry, label = "") {
  if (!("CHR" %in% names(dt))) dt[, CHR := NA_integer_]
  if (!("BP"  %in% names(dt))) dt[, BP  := NA_integer_]
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

# Replace non-rsID SNP column with rsID via chr:bp → BIM lookup. Used for
# sources whose SNP column is `chr:pos:SNP` or `chr:pos:ref:alt` (ICBP,
# Sinnott-Armstrong, AGEN, Pan-UKB). Assumes the source is GRCh37 (true for
# all current non-rsID sources here; 1000G phase 3 PLINK panels are GRCh37).
# No-op if > 50 % of rows already look like rsIDs. Modifies dt in-place.
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

log_provenance <- function(file_path, source_type, study, ancestry, notes = NULL) {
  # Append a short single-line provenance summary alongside the main
  # data_snapshot_log.md (which is the human-curated catalogue).
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  sha <- if (file.exists(file_path)) sha256_file(file_path) else NA_character_
  sz  <- if (file.exists(file_path)) file.info(file_path)$size else NA_integer_
  line <- sprintf("%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s",
                  ts, file_path, source_type, study, ancestry,
                  as.integer(sz), sha, ifelse(is.null(notes), "", notes))
  log_path <- repo_path("osf", "data_provenance.tsv")
  if (!file.exists(log_path)) {
    cat("timestamp\tpath\tsource\tstudy\tancestry\tsize_bytes\tsha256\tnotes\n",
        file = log_path)
  }
  cat(line, "\n", file = log_path, append = TRUE)
  invisible(line)
}
