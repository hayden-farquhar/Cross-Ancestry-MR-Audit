"""Main-text Figure 3 (JHG render): leave-one-locus-out decomposition of the
urate->gout cross-ancestry contrast (review-response analysis R/14).

Two panels:
  (A) IVW (mult. RE SE) EUR vs EAS estimates on the per-1-SD-urate log-odds scale
      for Full / drop-ABCG2 / drop-SLC2A9 / drop-both, with the cross-ancestry
      Cochran's Q annotated per scenario (closed-form from the two arms).
  (B) Per-locus single-SNP Wald ratios at ABCG2 and SLC2A9 (EUR vs EAS) against
      the genome-wide IVW reference, showing ABCG2 elevated in BOTH ancestries.

Tier-1 portfolio figure standard: SciencePlots + Okabe-Ito colourblind-safe palette
(EUR = blue #0072B2, EAS = vermilion #D55E00). Vector PDF + 300-dpi PNG.

Reads:  outputs/tables/urate_locus_leaveout.csv, urate_locus_wald.csv
Writes: outputs/figures/figure_loo_urate.{png,pdf}
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

try:
    import scienceplots  # noqa: F401
    plt.style.use(["science", "nature", "no-latex"])
except Exception:
    pass

EUR_C, EAS_C, GREY = "#0072B2", "#D55E00", "#999999"
ROOT = Path(__file__).resolve().parents[1]
lo = pd.read_csv(ROOT / "outputs/tables/urate_locus_leaveout.csv")
wd = pd.read_csv(ROOT / "outputs/tables/urate_locus_wald.csv")

scen = [("full", "Full instrument set"), ("drop_ABCG2", "Drop ABCG2"),
        ("drop_SLC2A9", "Drop SLC2A9"), ("drop_both", "Drop ABCG2 + SLC2A9")]

def arm(anc, sc):
    r = lo[(lo.ancestry == anc) & (lo.scenario == sc)].iloc[0]
    return r.b, r.se

fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.4, 4.0), gridspec_kw={"width_ratios": [1.55, 1]})
OFF = 0.18

# ---- Panel A: leave-one-locus-out forest ----
ys = (np.arange(len(scen)) * 1.0)[::-1]
for y, (sc, lab) in zip(ys, scen):
    bE, sE = arm("EUR", sc); bA, sA = arm("EAS", sc)
    Q = (bA - bE) ** 2 / (sE ** 2 + sA ** 2)            # closed-form cross-ancestry Q
    for (b, s, c, off) in [(bE, sE, EUR_C, OFF), (bA, sA, EAS_C, -OFF)]:
        axA.plot([b - 1.96 * s, b + 1.96 * s], [y + off, y + off], color=c, lw=1.6, solid_capstyle="butt", clip_on=False)
        axA.scatter([b], [y + off], s=38, color=c, edgecolor="black", lw=0.4, marker="s", zorder=4, clip_on=False)
    axA.text(3.05, y, f"Q = {Q:.0f}   ratio {bA/bE:.2f}×", va="center", ha="left", fontsize=7.0, color="black")

axA.axvline(0, color="black", ls="--", lw=0.6, alpha=0.7)
axA.set_yticks(ys); axA.set_yticklabels([l for _, l in scen], fontsize=8.5)
axA.set_ylim(-0.55, len(scen) - 0.45)
axA.set_xlim(-0.2, 4.5); axA.set_xticks([0, 1, 2, 3])
axA.set_xlabel("β (log-odds gout per 1-SD urate)", fontsize=8.5)
axA.set_title("A  Leave-one-locus-out (IVW)", fontsize=9.5, loc="left")
axA.scatter([], [], color=EUR_C, marker="s", s=38, label="European")
axA.scatter([], [], color=EAS_C, marker="s", s=38, label="East Asian")
axA.legend(fontsize=7.5, loc="lower right", frameon=False, bbox_to_anchor=(1.0, -0.02))

# ---- Panel B: per-locus Wald ratios vs genome-wide reference ----
gw = {"EUR": arm("EUR", "full")[0], "EAS": arm("EAS", "full")[0]}
loci = ["ABCG2", "SLC2A9"]
yb = (np.arange(len(loci)) * 1.0)[::-1]
for y, loc in zip(yb, loci):
    for anc, c, off in [("EUR", EUR_C, OFF), ("EAS", EAS_C, -OFF)]:
        r = wd[(wd.ancestry == anc) & (wd.locus == loc)].iloc[0]
        axB.plot([r.wald_b - 1.96 * r.wald_se, r.wald_b + 1.96 * r.wald_se], [y + off, y + off],
                 color=c, lw=1.6, solid_capstyle="butt")
        axB.scatter([r.wald_b], [y + off], s=38, color=c, edgecolor="black", lw=0.4, marker="o", zorder=4)
axB.axvline(gw["EUR"], color=EUR_C, ls=":", lw=0.9, alpha=0.9)
axB.axvline(gw["EAS"], color=EAS_C, ls=":", lw=0.9, alpha=0.9)
axB.set_yticks(yb); axB.set_yticklabels(loci, fontsize=8.5)
axB.set_ylim(-0.55, len(loci) - 0.45)
axB.set_xlim(0.4, 3.3); axB.set_xticks([1, 2, 3])
axB.set_xlabel("Single-SNP Wald β (log-odds)", fontsize=8.5)
axB.set_title("B  Per-locus Wald ratios", fontsize=9.5, loc="left")
axB.text(gw["EUR"], len(loci) - 0.62, "GW\nEUR", color=EUR_C, fontsize=6.2, ha="center", va="top")
axB.text(gw["EAS"], len(loci) - 0.62, "GW\nEAS", color=EAS_C, fontsize=6.2, ha="center", va="top")

fig.tight_layout(pad=0.6)
out_png = ROOT / "outputs/figures/figure_loo_urate.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(ROOT / "outputs/figures/figure_loo_urate.pdf", bbox_inches="tight")
print(f"wrote {out_png} and .pdf")
