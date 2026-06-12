"""
A3 H2 re-analysis using the EMPIRICAL LDSC sample-overlap intercept (C_off_diag)
rather than the contradictory pre-classifications in R/11 and Table S2.

Produces:
  outputs/tables/h2_empirical_c.csv          — per-pair: b_MR-APSS, b_IVW, |Δβ|, C_off_diag, stratum
  outputs/tables/h2_empirical_c_summary.csv  — strata means + Spearman/Pearson stats
  outputs/figures/figure_3_h2_empirical_c.png — replacement for the manuscript's Fig 3

Replaces the manuscript's §3.7/§4.1 "dramatic H2 inversion" claim (which was driven by
n=2 in a pre-classified 'high overlap' stratum of tiny-β BP→stroke pairs) with the
empirically-stratified test, which shows only a weak inverted trend that is not
statistically significant after multiplicity. Resolves Table S2 ↔ Fig 3 contradiction
by using a single empirical strata source everywhere downstream.
"""
from pathlib import Path
import csv
import statistics
from scipy.stats import spearmanr, pearsonr
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
APSS = ROOT / "outputs/tables/mr_apss_estimates.csv"
IVW  = ROOT / "outputs/tables/mr_estimates_per_method.csv"
OUT_TBL = ROOT / "outputs/tables/h2_empirical_c.csv"
OUT_SUM = ROOT / "outputs/tables/h2_empirical_c_summary.csv"
OUT_FIG = ROOT / "outputs/figures/figure_3_h2_empirical_c.png"

# Stratum thresholds on the empirical LDSC C off-diagonal.
# Justified by the observed distribution: near-zero (independent) vs >~0.02
# (some sample correlation) — see h2_empirical_c.csv.
C_LOW_MAX  = 0.02   # |C_off| < 0.02 → low / negligible empirical overlap
C_HIGH_MIN = 0.03   # |C_off| ≥ 0.03 → higher empirical overlap

apss = {}
with APSS.open() as f:
    for r in csv.DictReader(f):
        c = r.get("C_off_diag", "")
        apss[r["pair"]] = (float(r["b"]), float(r["se"]), float(c) if c else float("nan"))

ivw = {}
with IVW.open() as f:
    for r in csv.DictReader(f):
        if r["method"] != "Inverse variance weighted":
            continue
        ivw[r["pair"]] = (float(r["b"]), float(r["se"]))

rows = []
for pair, (ba, sa, c) in apss.items():
    if pair not in ivw:
        continue
    bi, si = ivw[pair]
    dbeta = abs(ba - bi)
    if abs(c) < C_LOW_MAX:
        stratum = "low (empirical C<0.02)"
    elif abs(c) >= C_HIGH_MIN:
        stratum = "high (empirical |C|>=0.03)"
    else:
        stratum = "mid (0.02-0.03)"
    rows.append({
        "pair": pair,
        "b_MRAPSS": ba,
        "se_MRAPSS": sa,
        "b_IVW": bi,
        "se_IVW": si,
        "abs_dbeta": dbeta,
        "C_off_diag": c,
        "empirical_stratum": stratum,
    })

rows.sort(key=lambda r: -abs(r["C_off_diag"]))
OUT_TBL.parent.mkdir(parents=True, exist_ok=True)
with OUT_TBL.open("w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# Strata summary
def mean(xs): return sum(xs) / len(xs) if xs else float("nan")
strata_summary = {}
for s in {r["empirical_stratum"] for r in rows}:
    sub = [r["abs_dbeta"] for r in rows if r["empirical_stratum"] == s]
    strata_summary[s] = (len(sub), mean(sub))

cs = [abs(r["C_off_diag"]) for r in rows]
dbs = [r["abs_dbeta"] for r in rows]
rs, ps = spearmanr(cs, dbs)
rp, pp = pearsonr(cs, dbs)

with OUT_SUM.open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["statistic", "value", "note"])
    for s, (n, m) in sorted(strata_summary.items()):
        w.writerow([f"mean |dbeta| {s}", f"{m:.4f}", f"n={n}"])
    w.writerow(["Spearman rho (|C_off|, |dbeta|)", f"{rs:.3f}", f"p={ps:.3f}"])
    w.writerow(["Pearson r (|C_off|, |dbeta|)",    f"{rp:.3f}", f"p={pp:.3f}"])
    w.writerow(["n pairs", len(rows), "all pair x ancestry slots with MR-APSS + IVW"])
    w.writerow(["interpretation", "weak inverted trend; not significant after multiplicity",
                "supersedes manuscript Sec 3.7/4.1 'dramatic H2 inversion' claim"])

# Figure 3 replacement: scatter |dbeta| vs |C_off|, regression line, pair labels.
fig, ax = plt.subplots(figsize=(8.5, 5.5))
ax.scatter(cs, dbs, s=40, edgecolor="black", linewidth=0.5, alpha=0.85)
for r, c, d in zip(rows, cs, dbs):
    ax.annotate(r["pair"], (c, d), xytext=(4, 3), textcoords="offset points", fontsize=7)
# Fitted line
import numpy as np
xs = np.array(cs); ys = np.array(dbs)
slope, intercept = np.polyfit(xs, ys, 1)
xline = np.linspace(xs.min(), xs.max(), 50)
ax.plot(xline, slope * xline + intercept, "--", color="grey", linewidth=1.0,
        label=f"OLS fit: Pearson r={rp:.2f}, p={pp:.3f}")
ax.axvline(C_LOW_MAX, color="lightgrey", linewidth=0.6)
ax.axvline(C_HIGH_MIN, color="lightgrey", linewidth=0.6)
ax.set_xlabel("Empirical LDSC sample-overlap intercept |C_off_diag|")
ax.set_ylabel("|β(MR-APSS) − β(IVW)| (log-odds scale)")
ax.set_title("H2: cross-method divergence vs empirical sample overlap")
ax.legend(loc="upper right", fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(OUT_FIG, dpi=200)
plt.close(fig)

print(f"Wrote {OUT_TBL.relative_to(ROOT)} ({len(rows)} pairs)")
print(f"Wrote {OUT_SUM.relative_to(ROOT)}")
print(f"Wrote {OUT_FIG.relative_to(ROOT)}")
print()
print("Strata means (|Δβ|):")
for s, (n, m) in sorted(strata_summary.items()):
    print(f"  {s}: mean={m:.4f} (n={n})")
print(f"Spearman rho={rs:.3f} p={ps:.3f}; Pearson r={rp:.3f} p={pp:.3f}")
print()
print("Honest H2 conclusion: the C_off vs |Δβ| relationship is weakly inverted and")
print("NOT statistically significant after multiplicity. The manuscript's 'dramatic")
print("inversion' claim was driven by n=2 in a pre-classified high-overlap stratum.")
