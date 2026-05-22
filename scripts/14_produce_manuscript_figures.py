#!/usr/bin/env python3
"""Produce the three headline manuscript figures.

Figure 1: PRISMA 2020 flow diagram (pair × ancestry attrition).
Figure 2: Headline forest plot — three selected pairs (HDL→CAD, LDL→LAS, HBA1C→T2D cross-ancestry).
Figure 3: H2 sample-overlap strip plot showing |β_APSS − β_IVW| per slot.

All outputs to outputs/figures/.
"""

from pathlib import Path
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "results" / "figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# ----------------------------------------------------------------------------
# Figure 1 — PRISMA 2020 flow diagram
# ----------------------------------------------------------------------------
def figure_1_prisma():
    fig, ax = plt.subplots(figsize=(8.5, 9.5))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 14)
    ax.axis("off")

    def box(x, y, w, h, text, fc="#eef2f7", ec="#2a3f5f", fontsize=10):
        b = FancyBboxPatch(
            (x - w / 2, y - h / 2),
            w, h,
            boxstyle="round,pad=0.08,rounding_size=0.15",
            linewidth=1.2, edgecolor=ec, facecolor=fc,
        )
        ax.add_patch(b)
        ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, wrap=True)

    def arrow(x1, y1, x2, y2):
        a = FancyArrowPatch(
            (x1, y1), (x2, y2),
            arrowstyle="->,head_length=8,head_width=6",
            linewidth=1.2, color="#2a3f5f",
        )
        ax.add_patch(a)

    def side_note(x, y, text):
        ax.text(x, y, text, ha="left", va="center", fontsize=9, style="italic", color="#6b7c93")

    # Title
    ax.text(5, 13.5, "Figure 1. PRISMA 2020 pair-selection flow",
            ha="center", va="center", fontsize=12, fontweight="bold")

    # Box 1 — initial candidates
    box(5, 12, 6.5, 1.0,
        "Initial candidate slots\n10 canonical pairs × 3 ancestries\n+ 3 negative-control slots = 33 slots")

    arrow(5, 11.5, 5, 10.6)
    side_note(8.6, 11.0, "Pair-selection inclusion criteria\n(§2.2): ≥3 prior MR papers,\n≥20 GWS instruments, ≥3,000\ncases in ≥2 ancestries, policy\nrelevance.")

    # Box 2 — pre-registered powered
    box(5, 10.1, 6.5, 1.0,
        "Pre-registered POWERED or LIMITED\n21 slots (after §2.7 power floor)")

    arrow(5, 9.6, 5, 8.7)
    side_note(8.6, 9.1, "Power floor (§2.7):\noutcome cases ≥ 3,000 AND\ndetectable OR within [1/1.5, 1.5]\nat 80% power, α = 0.05/K.\nExcluded: 8 AFR slots\n+ 1 EAS-lipid slot")

    # Box 3 — controlled-access excluded
    box(5, 8.2, 6.5, 1.2,
        "DEFERRED non-EUR outcome cohorts removed\n5 EAS lipid → cardiovascular slots\n(BBJ outcome data controlled-access at OSF freeze)")

    arrow(5, 7.5, 5, 6.6)
    side_note(8.6, 7.1, "DEFERRED slots cannot enter\nthe H1 denominator post-hoc\nper registration §17.5 even if\nthe data become available later.")

    # Box 4 — final analytical
    box(5, 6.1, 6.5, 1.2,
        "Final analytical set\n16 slots: 11 EUR + 1 EAS + 1 EUR neg-control\n(LDL × 5, HDL/TG/TC → CAD, BMI → T2D,\nurate → gout, BP × 2 → AIS, Lp(a) → CAVS,\nHbA1c → T2D EUR + EAS, NC-3 LDL → fracture)")

    arrow(5, 5.5, 5, 4.7)

    # Box 5 — H1 testing universe
    box(5, 4.1, 6.5, 1.2,
        "H1 testing universe = 1 pair\nHbA1c → T2D only (EUR + EAS)\nVerdict: transport-conditional\nQ = 0.006, p = 0.94, |Δβ| = 0.073 (>0.05 strict threshold)",
        fc="#ffe7d6", ec="#c97a4a")

    arrow(5, 3.5, 5, 2.5)

    # Box 6 — H2 descriptive
    box(5, 1.9, 6.5, 1.2,
        "H2 descriptive: pleiotropy, not overlap\nMean |β_APSS − β_IVW|: no-overlap 0.123 (n=10),\nmoderate 0.121 (n=4), high 0.003 (n=2).\nDirectional pre-spec falsified.",
        fc="#d6e9ff", ec="#3a6fa5")

    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure_1_prisma_flow.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Figure 2 — Headline forest plot (3 panels)
# ----------------------------------------------------------------------------
def figure_2_forest():
    methods = [
        "MR-APSS", "Inverse variance weighted", "MR Egger",
        "Weighted median", "Weighted mode",
        "MR-PRESSO (raw)", "MR-PRESSO (outlier-corrected)",
        "Radial-IVW",
    ]
    method_labels = {
        "MR-APSS": "MR-APSS (primary)",
        "Inverse variance weighted": "IVW",
        "MR Egger": "MR-Egger",
        "Weighted median": "Weighted median",
        "Weighted mode": "Weighted mode",
        "MR-PRESSO (raw)": "MR-PRESSO (raw)",
        "MR-PRESSO (outlier-corrected)": "MR-PRESSO (corrected)",
        "Radial-IVW": "Radial-IVW",
    }

    estimates = {}
    with open(ROOT / "results" / "tables" / "mr_estimates_per_method.csv") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                b = float(r["b"])
                se = float(r["se"])
            except (ValueError, KeyError):
                continue
            estimates.setdefault(r["pair"], {})[r["method"]] = (b, se)

    headline = [
        ("LDL_LAS_EUR", "LDL → large-artery stroke (EUR)"),
        ("HDL_CAD_EUR", "HDL → coronary artery disease (EUR)"),
        ("HBA1C_T2D_EAS", "HbA1c → type 2 diabetes (EAS)"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5), sharey=False)
    import math

    for ax, (pair, label) in zip(axes, headline):
        rows = []
        for m in methods:
            if m in estimates.get(pair, {}):
                b, se = estimates[pair][m]
                rows.append((method_labels[m], b, se))
        if not rows:
            ax.text(0.5, 0.5, f"No data for {pair}", ha="center", va="center", transform=ax.transAxes)
            continue

        y_positions = list(range(len(rows), 0, -1))
        ORs = [math.exp(b) for _, b, _ in rows]
        LCIs = [math.exp(b - 1.96 * se) for _, b, se in rows]
        UCIs = [math.exp(b + 1.96 * se) for _, b, se in rows]
        names = [r[0] for r in rows]

        for y, name, or_, lci, uci in zip(y_positions, names, ORs, LCIs, UCIs):
            color = "#c0392b" if name.startswith("MR-APSS") else "#2c3e50"
            ax.errorbar(or_, y,
                        xerr=[[or_ - lci], [uci - or_]],
                        fmt="o", color=color, ecolor=color, capsize=3, markersize=6)

        ax.axvline(1, linestyle="--", color="#888", linewidth=0.8)
        ax.set_yticks(y_positions)
        ax.set_yticklabels(names, fontsize=9)
        ax.set_xlabel("OR per 1-SD exposure (log scale)", fontsize=10)
        ax.set_xscale("log")

        all_lci = [v for v in LCIs if v > 0]
        all_uci = [v for v in UCIs if math.isfinite(v)]
        if all_lci and all_uci:
            xmin = max(0.1, min(all_lci) * 0.9)
            xmax = min(20, max(all_uci) * 1.1)
            ax.set_xlim(xmin, xmax)

        ax.set_title(label, fontsize=10, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Figure 2. Headline forest plots: MR-APSS vs sensitivity panel for three representative pairs",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure_2_headline_forest.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# Figure 3 — H2 strip plot |β_APSS − β_IVW| coloured by overlap stratum
# ----------------------------------------------------------------------------
def figure_3_h2_strip():
    rows = []
    with open(ROOT / "results" / "tables" / "mr_apss_vs_ivw.csv") as f:
        reader = csv.DictReader(f)
        for r in reader:
            try:
                rows.append({
                    "pair": r["pair"],
                    "ancestry": r["ancestry"],
                    "overlap": r["overlap"],
                    "abs_diff_b": float(r["abs_diff_b"]),
                })
            except (ValueError, KeyError):
                continue

    if not rows:
        raise SystemExit("No rows in mr_apss_vs_ivw.csv; cannot draw figure 3")

    rows.sort(key=lambda r: (-r["abs_diff_b"]))

    # x positions: cluster by stratum on a horizontal axis with jitter
    stratum_xpos = {"no": 1.0, "moderate": 2.0, "high": 3.0}
    stratum_color = {"no": "#1f77b4", "moderate": "#ff7f0e", "high": "#2ca02c"}
    stratum_label = {"no": "No / minimal\n(< 5%)", "moderate": "Moderate\n(5–30%)", "high": "High\n(> 30%)"}

    fig, ax = plt.subplots(figsize=(8.5, 6))
    import random
    random.seed(2026)

    plotted = {"no": [], "moderate": [], "high": []}
    for r in rows:
        x_base = stratum_xpos[r["overlap"]]
        x = x_base + random.uniform(-0.2, 0.2)
        ax.scatter(x, r["abs_diff_b"],
                   s=70,
                   color=stratum_color[r["overlap"]],
                   alpha=0.8, edgecolor="white", linewidth=1)
        plotted[r["overlap"]].append(r["abs_diff_b"])

        # Label outliers (top 2 in each stratum + the high-stratum points)
        if r["overlap"] == "high" or r["abs_diff_b"] > 0.18:
            label = r["pair"].replace("_EUR", "").replace("_EAS", "").replace("_", " → ")
            ax.annotate(label, (x, r["abs_diff_b"]),
                        xytext=(5, 5), textcoords="offset points",
                        fontsize=8, color="#444")

    # Mean lines per stratum
    for strat, values in plotted.items():
        if not values:
            continue
        mean = sum(values) / len(values)
        ax.plot([stratum_xpos[strat] - 0.3, stratum_xpos[strat] + 0.3], [mean, mean],
                color=stratum_color[strat], linewidth=3, alpha=0.6)
        ax.text(stratum_xpos[strat], mean + 0.012,
                f"mean = {mean:.3f}\n(n = {len(values)})",
                ha="center", va="bottom", fontsize=9,
                color=stratum_color[strat], fontweight="bold")

    ax.set_xticks(list(stratum_xpos.values()))
    ax.set_xticklabels([stratum_label[s] for s in ["no", "moderate", "high"]], fontsize=10)
    ax.set_xlabel("Pre-classified sample-overlap stratum", fontsize=11)
    ax.set_ylabel("|β_MR-APSS − β_IVW|", fontsize=11)
    ax.set_title("Figure 3. H2 descriptive: MR-APSS vs IVW divergence inverts the pre-specified gradient.\n"
                 "High-overlap pairs (BP → AIS, n=2) show essentially no divergence; no-overlap pairs (n=10) average 0.123.",
                 fontsize=10.5, loc="left", pad=12)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle=":", linewidth=0.6, alpha=0.5)
    ax.set_ylim(bottom=-0.02)

    plt.tight_layout()
    plt.savefig(OUT_DIR / "figure_3_h2_strip.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    figure_1_prisma()
    print(f"Wrote {OUT_DIR / 'figure_1_prisma_flow.png'}")
    figure_2_forest()
    print(f"Wrote {OUT_DIR / 'figure_2_headline_forest.png'}")
    figure_3_h2_strip()
    print(f"Wrote {OUT_DIR / 'figure_3_h2_strip.png'}")
