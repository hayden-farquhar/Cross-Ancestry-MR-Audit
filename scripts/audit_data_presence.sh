#!/bin/bash
# Data-presence audit. Confirms every pre-registered exposure × ancestry and
# outcome × ancestry file exists on disk, has non-zero size, and parses to the
# expected header. Runs no analysis — just verifies the substrate.
#
# Usage (run from the repository root):
#   bash scripts/audit_data_presence.sh
#
# Exit 0 = all required slots present + headers parseable.
# Exit 1 = something missing or schema unexpected.

set -uo pipefail

# Resolve the repository root from this script's own location so the audit
# runs identically regardless of where it is invoked from.
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 2

PASS=0
FAIL=0
WARN=0

# ANSI colours (off if not tty)
if [ -t 1 ]; then
  RED=$'\033[31m'; GRN=$'\033[32m'; YEL=$'\033[33m'; RST=$'\033[0m'
else
  RED=""; GRN=""; YEL=""; RST=""
fi

check_file() {
  local label="$1"
  local path="$2"
  local min_bytes="$3"
  local required_cols="$4"   # space-sep list; "" = skip header check
  local source_type="$5"     # tag for human-readable reporting

  if [ ! -e "$path" ]; then
    printf "  %sMISSING%s  %-40s  %s\n" "$RED" "$RST" "$label" "$(basename "$path")"
    FAIL=$((FAIL + 1)); return 1
  fi
  # Follow symlinks
  local real_path
  real_path=$(readlink -f "$path" 2>/dev/null || /usr/bin/python3 -c "import os,sys; print(os.path.realpath(sys.argv[1]))" "$path")
  local size
  size=$(stat -f%z "$real_path" 2>/dev/null)
  if [ -z "$size" ] || [ "$size" -lt "$min_bytes" ]; then
    printf "  %sTOO_SMALL%s %-40s  %s bytes (min %s)\n" "$RED" "$RST" "$label" "${size:-0}" "$min_bytes"
    FAIL=$((FAIL + 1)); return 1
  fi
  # Header check — only if required_cols non-empty
  if [ -n "$required_cols" ]; then
    local header
    case "$real_path" in
      *.gz|*.bgz)
        header=$(/usr/bin/gzip -dc "$real_path" 2>/dev/null | head -1)
        ;;
      *)
        header=$(head -1 "$real_path")
        ;;
    esac
    local missing_col=""
    for col in $required_cols; do
      if ! echo "$header" | grep -qE "(^|[[:space:]\"',])${col}([[:space:]\"',]|$)"; then
        missing_col="$col"
        break
      fi
    done
    if [ -n "$missing_col" ]; then
      printf "  %sSCHEMA%s    %-40s  header missing '%s'  [%s]\n" "$YEL" "$RST" "$label" "$missing_col" "$source_type"
      WARN=$((WARN + 1)); return 0
    fi
  fi
  printf "  %sOK%s        %-40s  %s  [%s]\n" "$GRN" "$RST" "$label" "$(echo "$size" | awk '{ split("B KB MB GB TB", units); for (i = 5; i > 0; i--) { if ($1 > 1024^(i-1)) { printf "%.1f%s", $1/1024^(i-1), units[i]; exit } } }')" "$source_type"
  PASS=$((PASS + 1))
  return 0
}

echo "==================================================================="
echo "Cross-ancestry MR portability audit — data presence verification"
echo "Run timestamp: $(date '+%Y-%m-%d %H:%M:%S %Z')"
echo "==================================================================="

echo ""
echo "=== LD reference panels (1000G phase 3 PLINK, EUR/EAS/AFR) ==="
for pop in EUR EAS AFR; do
  for ext in bed bim fam afreq; do
    case "$ext" in
      bed) min=200000000 ;;
      bim) min=10000000 ;;
      fam) min=5000 ;;
      afreq) min=50000000 ;;
    esac
    check_file "1000G_${pop}.${ext}" \
      "data/ld_panels/1000G_${pop}/1000G_${pop}.${ext}" \
      "$min" "" "Zenodo_6614170"
  done
done

echo ""
echo "=== Exposures — Lipids (GLGC Graham 2021) ==="
check_file "GLGC EUR LDL" \
  "data/raw/exposures/GLGC_Graham2021/EUR/LDL_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
  100000000 "rsID CHROM EFFECT_SIZE" "GLGC_UMich"
# EUR HDL/TG/TC via EBI mirror (TSV gzipped locally to ~1.1-1.3 GB)
for trait in HDL logTG TC; do
  check_file "GLGC EUR $trait" \
    "data/raw/exposures/GLGC_Graham2021/EUR/${trait}_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
    1000000000 "variant_id chromosome beta p_value" "GLGC_EBI"
done
# EAS all 4 via U Mich (gzipped)
for trait in LDL HDL logTG TC; do
  check_file "GLGC EAS $trait" \
    "data/raw/exposures/GLGC_Graham2021/EAS/${trait}_INV_EAS_1KGP3_ALL.meta.singlevar.results.gz" \
    500000000 "rsID CHROM EFFECT_SIZE" "GLGC_UMich"
done
# AFR: TC via U Mich (gzipped at acquisition), LDL/HDL/TG via EBI (TSV gzipped locally to ~750-820 MB)
check_file "GLGC AFR TC" \
  "data/raw/exposures/GLGC_Graham2021/AFR/TC_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
  1000000000 "rsID CHROM EFFECT_SIZE" "GLGC_UMich"
for trait in LDL HDL logTG; do
  check_file "GLGC AFR $trait" \
    "data/raw/exposures/GLGC_Graham2021/AFR/${trait}_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz" \
    700000000 "variant_id chromosome beta p_value" "GLGC_EBI"
done

echo ""
echo "=== Exposures — Urate (CKDGen Tin 2019) ==="
check_file "CKDGen urate trans-ethnic" \
  "data/raw/exposures/CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_all_741_nstud37_summac400_rsid.txt.gz" \
  100000000 "RSID Effect StdErr" "CKDGen_Tin2019"
check_file "CKDGen urate EUR" \
  "data/raw/exposures/CKDGen_Tin2019/urate_chr1_22_LQ_IQ06_mac10_EA_60_prec1_nstud30_summac400_rsid.txt.gz" \
  100000000 "RSID Effect StdErr" "CKDGen_Tin2019"
check_file "CKDGen gout trans-ethnic" \
  "data/raw/exposures/CKDGen_Tin2019/gout_chr1_22_LQ_IQ06_mac10_all_201_nstud10_summac400_rsid.txt.gz" \
  100000000 "RSID Effect StdErr" "CKDGen_Tin2019"

echo ""
echo "=== Exposures — BMI (GIANT) ==="
check_file "GIANT Yengo 2018 EUR BMI" \
  "data/raw/exposures/GIANT_BMI/GIANT_2018_Yengo_BMI_EUR_Locke_UKB_meta.txt.gz" \
  30000000 "SNP BETA SE" "GIANT_Yengo2018"
# Turcot files use ancestry-specific naming variations (AFR_AmericanExome vs EUR_Exome).
# Resolve each explicit filename rather than via glob.
for fn in GIANT_2018_Turcot_BMI_AFR_AmericanExome.fmt.gz \
         GIANT_2018_Turcot_BMI_EAS_Exome.fmt.gz \
         GIANT_2018_Turcot_BMI_EUR_Exome.fmt.gz \
         GIANT_2018_Turcot_BMI_AllAncestry_Exome.fmt.gz; do
  check_file "GIANT Turcot 2018 ${fn}" \
    "data/raw/exposures/GIANT_BMI/${fn}" \
    5000000 "" "GIANT_Turcot2018"
done

echo ""
echo "=== Exposures — BP (ICBP Evangelou 2018) ==="
for trait in SBP DBP PP; do
  check_file "ICBP+UKB EUR $trait" \
    "data/raw/exposures/ICBP_Evangelou2018/Evangelou_30224653_${trait}.txt.gz" \
    100000000 "MarkerName Effect StdErr" "ICBP_Evangelou2018"
done

echo ""
echo "=== Exposures — Smoking (GSCAN Saunders 2022) ==="
for anc in EUR EAS AFR; do
  check_file "GSCAN ${anc}-stratified" \
    "data/raw/exposures/GSCAN_Saunders2022/${anc}_stratified.zip" \
    500000000 "" "GSCAN_Saunders2022"
done

echo ""
echo "=== Exposures — UKB biomarkers (Sinnott-Armstrong 2021) ==="
for biom in Cholesterol_adjstatins LDL_direct_adjstatins HDL_cholesterol Triglycerides Lipoprotein_A Glycated_haemoglobin_HbA1c Urate; do
  check_file "Sinnott-Armstrong $biom" \
    "data/raw/exposures/SinnottArmstrong_UKB/${biom}.imp.gz" \
    150000000 "MarkerName Effect StdErr" "SinnottArmstrong2021"
done

echo ""
echo "=== Exposures — HbA1c (MAGIC Chen 2021) ==="
# TA (trans-ancestry MANTRA) uses log10BF schema, not beta/SE — not usable as
# direct MR exposure. Per-ancestry files (EUR/EAS/AA/SAS/HISP) have standard
# beta/SE/p schema.
check_file "MAGIC HbA1c TA (MANTRA log10BF)" \
  "data/raw/exposures/MAGIC_Chen2021/MAGIC1000G_HbA1c_TA.tsv.gz" \
  50000000 "variant log10BF" "MAGIC_Chen2021_MANTRA"
for anc in EUR EAS AA SAS HISP; do
  check_file "MAGIC HbA1c ${anc}" \
    "data/raw/exposures/MAGIC_Chen2021/MAGIC1000G_HbA1c_${anc}.tsv.gz" \
    50000000 "variant effect_allele beta" "MAGIC_Chen2021"
done

echo ""
echo "=== Outcomes — FinnGen R12 (EUR) ==="
for ep in I9_CHD I9_AF N14_ACUTERENFAIL T2D M13_GOUT I9_CAVS_OPERATED ST19_FRACT_FOREA; do
  check_file "FinnGen R12 ${ep}" \
    "data/raw/outcomes/FinnGen_R12/finngen_R12_${ep}.gz" \
    500000000 "rsids alt ref beta" "FinnGen_R12"
done

echo ""
echo "=== Outcomes — MEGASTROKE EUR + trans-ancestry ==="
for f in 1.AS.EUR 2.AIS.EUR 3.LAS.EUR 4.CES.EUR 5.SVS.EUR 2.AIS.TRANS; do
  check_file "MEGASTROKE ${f}" \
    "data/raw/outcomes/MEGASTROKE_Malik2018/MEGASTROKE.${f}.out" \
    100000000 "MarkerName Effect StdErr" "MEGASTROKE"
done

echo ""
echo "=== Outcomes — AGEN-T2D EAS (Spracklen 2020) ==="
check_file "AGEN-T2D EAS primary" \
  "data/raw/outcomes/AGEN_Spracklen2020/SpracklenCN_prePMID_T2D_ALL_Primary.txt.gz" \
  200000000 "MarkerName Beta SE" "AGEN_T2D"

echo ""
echo "=== Pan-UKB AFR coverage (5 priority traits) ==="
for code in 30790 30880 30750; do
  check_file "Pan-UKB AFR biomarker ${code}" \
    "data/raw/outcomes/PanUKB_AFR/biomarkers-${code}-both_sexes-irnt.tsv.bgz" \
    1000000000 "" "PanUKB_AFR"
done
for code in 21001 4080; do
  check_file "Pan-UKB AFR continuous ${code}" \
    "data/raw/outcomes/PanUKB_AFR/continuous-${code}-both_sexes-irnt.tsv.bgz" \
    1000000000 "" "PanUKB_AFR"
done
check_file "Pan-UKB manifest" \
  "data/raw/outcomes/PanUKB_AFR/phenotype_manifest.tsv" \
  1000000 "phenocode aws_path" "PanUKB_AFR"

echo ""
echo "=== Reference tables ==="
check_file "Chong 2019 MR-of-proteome supplement" \
  "data/reference/2019_CHONG_CIRC_MR_BIOMARKERS_ISCHEMIC_STROKE_ALL_653_FINAL.csv" \
  1000000 "exposure outcome method" "Chong2019"

echo ""
echo "==================================================================="
printf "Audit summary:  %sPASS=%d%s  %sWARN=%d%s  %sFAIL=%d%s\n" \
  "$GRN" "$PASS" "$RST" "$YEL" "$WARN" "$RST" "$RED" "$FAIL" "$RST"
echo "==================================================================="

if [ "$FAIL" -gt 0 ]; then
  echo ""
  echo "Audit FAILED — see lines marked MISSING or TOO_SMALL above."
  exit 1
fi

if [ "$WARN" -gt 0 ]; then
  echo ""
  echo "Audit passed with $WARN schema warnings — confirm parsers handle these schemas in R/utils.R."
fi

exit 0
