"""Generate the new main-text Figure 2: urate->gout cross-ancestry forest plot.

Replaces the prior representative-panels figure (which depicted LDL->LAS,
HDL->CAD, HbA1c->T2D EAS — chosen for the old K=1 + H2-inversion headline).

Shows three MR-APSS rows on the per-1-SD-urate log-odds scale:
  1. URATE_GOUT_EUR        — canonical EUR pair (EUR urate -> FinnGen gout)
  2. URATE_GOUT_TRANS_EUR  — A3 identical-exposure sensitivity (trans-urate -> FinnGen gout)
  3. URATE_GOUT_EAS        — EAS pair (trans-urate -> Major 2024 EAS gout)

Plus annotated Q-test results for the confounded and de-confounded comparisons.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

ROOT = Path(__file__).resolve().parents[1]
OUT  = ROOT / "outputs/figures/figure_2_urate_gout_cross_ancestry.png"

# (label, b, se, n, role)  — values from outputs/tables/mr_apss_estimates.csv
rows = [
    ("URATE_GOUT_EUR\n(EUR-urate → FinnGen gout, EUR)\nCanonical pair",
     1.21267, 0.12375, 63, "primary_eur"),
    ("URATE_GOUT_TRANS_EUR\n(trans-urate → FinnGen gout, EUR)\nA3 identical-exposure sensitivity",
     1.28889, 0.12181, 90, "sensitivity"),
    ("URATE_GOUT_EAS\n(trans-urate → Major 2024 gout, EAS)\nEAS arm",
     2.29088, 0.25644, 88, "primary_eas"),
]

# Q-tests
Q_confounded = 14.34
p_confounded = 1.5e-4
I2_confounded = 93.0

Q_deconfounded = 12.46
p_deconfounded = 4.2e-4
I2_deconfounded = 92.0

# Layout
fig, ax = plt.subplots(figsize=(11.5, 5.5))

y_positions = list(range(len(rows), 0, -1))  # top to bottom
colors = {"primary_eur": "#1f77b4", "sensitivity": "#7f7f7f", "primary_eas": "#d62728"}

for y, (lab, b, se, n, role) in zip(y_positions, rows):
    lo = b - 1.96*se
    hi = b + 1.96*se
    color = colors[role]
    ax.plot([lo, hi], [y, y], color=color, linewidth=2.0, solid_capstyle="butt", zorder=3)
    ax.scatter([b], [y], s=180, color=color, edgecolor="black", linewidth=0.7, zorder=4, marker="s")
    ax.text(b, y - 0.16, f"β={b:.3f} (95%CI {lo:.2f}–{hi:.2f}); n={n}",
            ha="center", va="top", fontsize=8.5, color=color)

# x-axis on log-OR scale
ax.axvline(0, color="black", linestyle="--", linewidth=0.7, alpha=0.7, zorder=2)
ax.set_xlabel("β (log-odds gout per 1-SD urate); vertical line marks null (β = 0)", fontsize=10)
# secondary axis for OR
def b2or(b): return np.exp(b)
def or2b(o): return np.log(o)
secax = ax.secondary_xaxis("top", functions=(b2or, or2b))
secax.set_xlabel("Odds ratio (OR) per 1-SD urate", fontsize=10)
secax.set_xticks([1, 2, 3, 5, 10, 20])
secax.set_xticklabels(["1", "2", "3", "5", "10", "20"])

ax.set_yticks(y_positions)
ax.set_yticklabels([r[0] for r in rows], fontsize=9)
ax.set_xlim(-0.1, 3.3)
ax.set_ylim(0.3, len(rows) + 1.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.tick_params(axis="x", labelsize=9)

# Q-test annotation box on the right side
q_text = (
    "Cross-ancestry Cochran's Q (MR-APSS primary):\n\n"
    f"  Confounded test (URATE_GOUT_EUR vs URATE_GOUT_EAS):\n"
    f"    Q = {Q_confounded:.1f}, p = {p_confounded:.1e}, I² = {I2_confounded:.0f}%\n\n"
    f"  De-confounded test (URATE_GOUT_TRANS_EUR vs URATE_GOUT_EAS):\n"
    f"    Q = {Q_deconfounded:.1f}, p = {p_deconfounded:.1e}, I² = {I2_deconfounded:.0f}%\n\n"
    f"  Discrepancy survives identical-exposure de-confound;\n"
    f"  ≈ 13% of Q attributable to exposure-instrument set."
)
ax.text(0.98, 0.98, q_text, transform=ax.transAxes, ha="right", va="top",
        fontsize=8.0, family="monospace",
        bbox=dict(facecolor="#fff8e1", edgecolor="#bdbdbd", boxstyle="round,pad=0.5"))

fig.tight_layout()
fig.savefig(OUT, dpi=300, bbox_inches="tight")
plt.close(fig)
print(f"Wrote {OUT.relative_to(ROOT)}")
