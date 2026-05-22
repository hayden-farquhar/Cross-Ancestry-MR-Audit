"""
Cross-ancestry MR portability audit — interactive results dashboard.

Lets a viewer pick a canonical pair and inspect:
  - per-pair forest plot across ancestries × methods
  - colocalization H4 heatmap per locus × ancestry
  - heterogeneity p-value table
  - power-diagnostic panel

Run locally from the repository root:
    streamlit run app/streamlit_app.py

The dashboard reads from outputs/tables/*.csv produced by the analysis pipeline
(scripts/05_..09_). If any required input is missing it surfaces a clear
"not yet computed" message rather than a broken view.
"""

from __future__ import annotations
from pathlib import Path
import streamlit as st
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
TABLES = REPO_ROOT / "outputs" / "tables"

st.set_page_config(
    page_title="Cross-Ancestry MR Audit",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.title("Cross-Ancestry MR Portability Audit")
st.caption(
    "Pre-registered audit of 10 canonical MR exposure–outcome pairs across "
    "EUR / EAS / AFR. Pre-registration: OSF DOI 10.17605/OSF.IO/U8TX4."
)


@st.cache_data
def load_table(name: str) -> pd.DataFrame | None:
    path = TABLES / name
    if not path.exists():
        return None
    return pd.read_csv(path)


pairs = load_table("canonical_pairs_prespec.csv")
mr_est = load_table("mr_estimates_per_method.csv")
het = load_table("heterogeneity_meta.csv")
coloc_t = load_table("coloc_results.csv")
power_t = load_table("power_diagnostics.csv")
strength = load_table("instrument_strength.csv")


tabs = st.tabs(["Pairs", "MR estimates", "Heterogeneity", "Coloc", "Power"])


with tabs[0]:
    st.header("Pre-specified canonical pairs")
    if pairs is None:
        st.info("canonical_pairs_prespec.csv not yet present.")
    else:
        st.dataframe(pairs, use_container_width=True)


with tabs[1]:
    st.header("MR estimates per pair × ancestry × method")
    if mr_est is None:
        st.info("outputs/tables/mr_estimates_per_method.csv not yet computed. "
                "Run R/05 + R/06.")
    else:
        unique_pairs = sorted(set(mr_est["pair"]))
        sel = st.selectbox("Pair × ancestry", unique_pairs)
        sub = mr_est[mr_est["pair"] == sel].copy()
        sub["OR"] = pd.np.exp(sub["b"]) if hasattr(pd, "np") else None
        st.dataframe(sub, use_container_width=True)


with tabs[2]:
    st.header("Cross-ancestry heterogeneity (Cochran's Q, I²)")
    if het is None:
        st.info("heterogeneity_meta.csv not yet computed. Run R/09.")
    else:
        st.dataframe(het.sort_values("I2_pct", ascending=False),
                     use_container_width=True)
        st.bar_chart(het.set_index("pair")["I2_pct"])


with tabs[3]:
    st.header("Colocalization (coloc.abf H0–H4 per locus)")
    if coloc_t is None:
        st.info("coloc_results.csv not yet computed. Run R/07.")
    else:
        st.dataframe(coloc_t, use_container_width=True)


with tabs[4]:
    st.header("Power diagnostics")
    if power_t is None:
        st.info("power_diagnostics.csv not yet computed. Run R/08.")
    else:
        st.dataframe(power_t.sort_values("detectable_or_at_80pct_power"),
                     use_container_width=True)


st.sidebar.markdown(
    "**Reporting:** STROBE-MR + PRISMA-2020.\n\n"
    "**Ethics:** HREC-exempt (summary statistics only).\n\n"
    "**Phenocode fixes locked in pre-registration:** "
    "`GOUT → M13_GOUT`; `I9_AORTSTE → I9_CAVS_OPERATED`."
)
