#!/usr/bin/env python3
"""
ADMET Analysis Workflow: pkCSM + SwissADME Integration
=======================================================
Author  : Sana Begum
Role    : Research Analyst, Era's Lucknow Medical College & Hospital
Email   : begumsana686@gmail.com
GitHub  : https://github.com/Sanabegum09

Description
-----------
This workflow integrates ADMET predictions from two widely used web servers:
  1. pkCSM  — https://biosig.lab.uq.edu.au/pkcsm/
     Predicts: Absorption (Caco-2, HIA, P-gp), Distribution (BBB, PPB, VDss),
               Metabolism (CYP inhibition/substrate), Excretion (renal OCT2),
               Toxicity (AMES mutagenicity, hERG, hepatotoxicity, LD50, skin sensitisation)

  2. SwissADME — https://www.swissadme.ch/
     Predicts: Physicochemical properties, Lipinski/Veber/Egan/Ghose rules,
               Bioavailability radar, GI absorption, BBB permeant, P-gp substrate,
               CYP inhibition (2C19, 2C9, 2D6, 3A4, 1A2), Water solubility (ESOL/Ali)

Usage in Computational Drug Discovery (Sana Begum's workflow)
--------------------------------------------------------------
  STEP 1: Run molecular docking (AutoDock Vina) → get top hits
  STEP 2: Prepare SMILES → submit to pkCSM & SwissADME (batch)
  STEP 3: Download CSV results from both servers
  STEP 4: Run this script to parse, integrate, filter, and visualise

Script Usage
------------
  python admet_pkcsm_swissadme.py \\
      --pkcsm   pkcsm_results.csv \\
      --swissadme swissadme_results.csv \\
      --smiles  compounds.smi \\
      --output  admet_integrated/ \\
      --plot

Prepare Input for Web Servers
------------------------------
  For pkCSM:
      - Go to https://biosig.lab.uq.edu.au/pkcsm/prediction
      - Select 'Batch' tab → paste SMILES (one per line, optionally with name)
      - Download: 'Download results as CSV'

  For SwissADME:
      - Go to https://www.swissadme.ch/
      - Paste SMILES list (one per line)
      - Run → Export → 'Download results (CSV)'

Requirements
------------
    pip install pandas numpy matplotlib seaborn rdkit

References
----------
    Pires DE et al. pkCSM: Predicting Small-Molecule Pharmacokinetic and
    Toxicity Properties Using Graph-Based Signatures. J Med Chem. 2015;58(9):4066-72.
    doi:10.1021/acs.jmedchem.5b00104

    Daina A et al. SwissADME: a free web tool to evaluate pharmacokinetics,
    drug-likeness and medicinal chemistry friendliness of small molecules.
    Sci Rep. 2017;7:42717. doi:10.1038/srep42717
"""

import os, sys, argparse, logging, warnings
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
warnings.filterwarnings("ignore")

# Optional RDKit for SMILES prep utility
try:
    from rdkit import Chem
    from rdkit.Chem import Draw, Descriptors
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False

logging.basicConfig(level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ── Colour scheme ────────────────────────────────────────────
BLUE   = "#1A5276"
RED    = "#C0392B"
GREEN  = "#1E8449"
ORANGE = "#CA6F1E"
GREY   = "#ABB2B9"


# ════════════════════════════════════════════════════════════
# STEP 0 — PREPARE INPUT FOR WEB SERVERS
# ════════════════════════════════════════════════════════════

def prepare_smiles_from_sdf(sdf_file: str, out_smi: str) -> pd.DataFrame:
    """
    Extract SMILES from an SDF file and save formatted input for pkCSM/SwissADME.

    Parameters
    ----------
    sdf_file : str   Path to SDF file (e.g. top-20 docking hits)
    out_smi  : str   Output SMILES file path

    Returns
    -------
    pd.DataFrame : Table of compound names and SMILES
    """
    if not RDKIT_AVAILABLE:
        log.warning("RDKit not installed — install via: pip install rdkit")
        return pd.DataFrame()

    suppl = Chem.SDMolSupplier(sdf_file, sanitize=True, removeHs=True)
    records = []
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        name   = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Compound_{i+1}"
        smiles = Chem.MolToSmiles(mol)
        mw     = round(Descriptors.MolWt(mol), 2)
        records.append({"Compound": name, "SMILES": smiles, "MW": mw})

    df = pd.DataFrame(records)
    with open(out_smi, "w") as fh:
        for _, row in df.iterrows():
            fh.write(f"{row['SMILES']} {row['Compound']}\n")
    log.info(f"Prepared {len(df)} SMILES → {out_smi}")
    log.info("Submit this file to pkCSM (batch) and SwissADME for ADMET analysis")
    return df


# ════════════════════════════════════════════════════════════
# STEP 1 — PARSE pkCSM RESULTS
# ════════════════════════════════════════════════════════════

# pkCSM standard column map (their CSV headers)
PKCSM_COLUMNS = {
    # Absorption
    "Caco2_Permeability":          "Caco2 permeability",      # log Papp (cm/s)  > -5.15 = good
    "Intestinal_Absorption_HIA":   "Human Intestinal Absorption",  # % absorbed, >30% = good
    "Pgp_Inhibitor":               "P-glycoprotein I inhibitor",
    "Pgp_Substrate":               "P-glycoprotein II substrate",
    # Distribution
    "VDss":                        "VDss (human)",            # L/kg, 0.04–20 = drug-like
    "BBB_Permeability":            "BBB permeability",        # log BB, >-1 = permeable
    "CNS_Permeability":            "CNS permeability",        # log PS, >-2 = CNS active
    "Plasma_Protein_Binding":      "Fraction unbound (human)",
    # Metabolism
    "CYP1A2_Inhibitor":            "CYP1A2 inhibitor",
    "CYP2C19_Inhibitor":           "CYP2C19 inhibitor",
    "CYP2C9_Inhibitor":            "CYP2C9 inhibitor",
    "CYP2D6_Inhibitor":            "CYP2D6 inhibitor",
    "CYP3A4_Inhibitor":            "CYP3A4 inhibitor",
    "CYP2D6_Substrate":            "CYP2D6 substrate",
    "CYP3A4_Substrate":            "CYP3A4 substrate",
    # Excretion
    "Renal_OCT2_Substrate":        "Renal OCT2 substrate",
    "Total_Clearance":             "Total Clearance",         # log ml/min/kg
    # Toxicity
    "AMES_Mutagenicity":           "AMES toxicity",           # Binary
    "Max_Tolerated_Dose":          "Max. tolerated dose (human)", # log mg/kg/day
    "hERG_Blockade":               "hERG I inhibitor",        # Cardiac safety
    "hERG_Blockade_II":            "hERG II inhibitor",
    "Oral_Rat_Acute_Toxicity_LD50":"Oral Rat Acute Toxicity (LD50)", # mol/kg
    "Oral_Rat_Chronic_Toxicity":   "Oral Rat Chronic Toxicity (LOAEL)",
    "Hepatotoxicity":              "Hepatotoxicity",
    "Skin_Sensitisation":          "Skin Sensitisation",
    "T_Pyriformis_Toxicity":       "T. pyriformis toxicity",
    "Minnow_Toxicity":             "Minnow toxicity",
}

def parse_pkcsm(csv_path: str) -> pd.DataFrame:
    """
    Parse pkCSM batch prediction CSV output.

    The CSV is downloaded directly from https://biosig.lab.uq.edu.au/pkcsm/prediction
    after batch submission.
    """
    df = pd.read_csv(csv_path)
    log.info(f"pkCSM: loaded {len(df)} compounds, {len(df.columns)} properties")

    # Standardise compound name column
    name_col = next((c for c in df.columns if "name" in c.lower() or "compound" in c.lower()
                     or "smiles" in c.lower()), df.columns[0])
    df = df.rename(columns={name_col: "Compound"})

    # Rename known columns
    reverse_map = {v: k for k, v in PKCSM_COLUMNS.items()}
    df = df.rename(columns={c: reverse_map.get(c, c) for c in df.columns})

    df["Source"] = "pkCSM"
    return df


# ════════════════════════════════════════════════════════════
# STEP 2 — PARSE SwissADME RESULTS
# ════════════════════════════════════════════════════════════

# SwissADME standard column map
SWISSADME_COLUMNS = {
    "SMILES":                    "SMILES",
    "MW":                        "Molecular weight",
    "LogP_SwissADME":            "Log Po/w (iLOGP)",       # or XLOGP3 or WLOGP
    "HBD":                       "H-bond donors",
    "HBA":                       "H-bond acceptors",
    "TPSA":                      "Topological Polar Surface Area (TPSA)",
    "Rotatable_Bonds":           "Rotatable bonds",
    "Lipinski_Pass":             "Lipinski #violations",    # 0 = pass
    "Veber_Pass":                "Veber #violations",
    "GI_Absorption":             "GI absorption",           # High / Low
    "BBB_Permeant":              "BBB permeant",            # Yes/No
    "Pgp_Substrate_Swiss":       "P-gp substrate",
    "CYP1A2_Inhibitor_Swiss":    "CYP1A2 inhibitor",
    "CYP2C19_Inhibitor_Swiss":   "CYP2C19 inhibitor",
    "CYP2C9_Inhibitor_Swiss":    "CYP2C9 inhibitor",
    "CYP2D6_Inhibitor_Swiss":    "CYP2D6 inhibitor",
    "CYP3A4_Inhibitor_Swiss":    "CYP3A4 inhibitor",
    "Water_Solubility_ESOL":     "Water solubility (ESOL)",  # log mol/L
    "Water_Solubility_Class":    "Water solubility class",   # Soluble/Insoluble
    "Bioavailability_Score":     "Bioavailability Score",    # 0.17, 0.55, 0.85
    "Drug_Likeness":             "Drug-likeness",
    "Medicinal_Chemistry":       "Medicinal Chemistry",
    "PAINS_Alerts_Swiss":        "PAINS alerts",
    "Brenk_Alerts_Swiss":        "Brenk alerts",
    "LeadLikeness_Swiss":        "Lead-likeness",
}

def parse_swissadme(csv_path: str) -> pd.DataFrame:
    """
    Parse SwissADME CSV output.

    Downloaded from http://www.swissadme.ch/ after batch run → Export CSV.
    """
    df = pd.read_csv(csv_path)
    log.info(f"SwissADME: loaded {len(df)} compounds, {len(df.columns)} properties")

    # Standardise compound identifier
    name_col = next((c for c in df.columns if "name" in c.lower() or
                     "compound" in c.lower() or "id" in c.lower()), df.columns[0])
    df = df.rename(columns={name_col: "Compound"})

    reverse_map = {v: k for k, v in SWISSADME_COLUMNS.items()}
    df = df.rename(columns={c: reverse_map.get(c, c) for c in df.columns})

    # Convert Lipinski violations: 0 violations = Pass
    if "Lipinski_Pass" in df.columns:
        df["Lipinski_Pass"] = df["Lipinski_Pass"].apply(
            lambda x: "Pass" if str(x).strip() in ("0", "0.0", "0 violation") else "Fail"
        )

    df["Source"] = "SwissADME"
    return df


# ════════════════════════════════════════════════════════════
# STEP 3 — INTEGRATE & FILTER
# ════════════════════════════════════════════════════════════

def integrate_admet(pkcsm_df: pd.DataFrame, swiss_df: pd.DataFrame) -> pd.DataFrame:
    """
    Merge pkCSM and SwissADME results on compound name.

    Returns a single integrated DataFrame with key properties from both servers.
    """
    # Select key pkCSM columns (those that actually exist)
    pkcsm_keep = ["Compound",
                  "Caco2_Permeability", "Intestinal_Absorption_HIA",
                  "Pgp_Substrate", "BBB_Permeability", "VDss",
                  "CYP1A2_Inhibitor", "CYP2C9_Inhibitor", "CYP2C19_Inhibitor",
                  "CYP2D6_Inhibitor", "CYP3A4_Inhibitor",
                  "AMES_Mutagenicity", "hERG_Blockade", "hERG_Blockade_II",
                  "Hepatotoxicity", "Skin_Sensitisation",
                  "Oral_Rat_Acute_Toxicity_LD50", "Max_Tolerated_Dose"]
    pkcsm_keep = [c for c in pkcsm_keep if c in pkcsm_df.columns]

    # Select key SwissADME columns
    swiss_keep = ["Compound",
                  "MW", "LogP_SwissADME", "HBD", "HBA", "TPSA", "Rotatable_Bonds",
                  "Lipinski_Pass", "GI_Absorption", "BBB_Permeant",
                  "CYP2D6_Inhibitor_Swiss", "CYP3A4_Inhibitor_Swiss",
                  "Water_Solubility_ESOL", "Water_Solubility_Class",
                  "Bioavailability_Score", "PAINS_Alerts_Swiss",
                  "Drug_Likeness", "LeadLikeness_Swiss"]
    swiss_keep = [c for c in swiss_keep if c in swiss_df.columns]

    merged = pd.merge(
        pkcsm_df[pkcsm_keep],
        swiss_df[swiss_keep],
        on="Compound", how="outer", suffixes=("_pkcsm", "_swiss")
    )
    log.info(f"Integrated: {len(merged)} compounds with combined ADMET properties")
    return merged


def apply_drug_filters(df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply standard drug-likeness and safety filters using both pkCSM and SwissADME data.

    Filter criteria (drug-like profile):
      ✓  Lipinski Ro5 = Pass (SwissADME)
      ✓  GI Absorption = High (SwissADME)
      ✓  AMES mutagenicity = Non-mutagen (pkCSM)
      ✓  hERG blockade I = Non-inhibitor (pkCSM)  [cardiac safety]
      ✓  Hepatotoxicity = Non-hepatotoxic (pkCSM)
      ✓  PAINS alerts = 0 (SwissADME)
    """
    df = df.copy()
    filters = {}

    if "Lipinski_Pass" in df.columns:
        filters["F1_Lipinski"] = df["Lipinski_Pass"].astype(str).str.lower().isin(["pass", "yes", "true"])

    if "GI_Absorption" in df.columns:
        filters["F2_GI_Absorption"] = df["GI_Absorption"].astype(str).str.lower().isin(["high"])

    if "AMES_Mutagenicity" in df.columns:
        filters["F3_AMES_Safe"] = df["AMES_Mutagenicity"].astype(str).str.lower().isin(
            ["no", "non-mutagen", "false", "0"])

    if "hERG_Blockade" in df.columns:
        filters["F4_hERG_Safe"] = df["hERG_Blockade"].astype(str).str.lower().isin(
            ["no", "non-inhibitor", "false", "0"])

    if "Hepatotoxicity" in df.columns:
        filters["F5_Liver_Safe"] = df["Hepatotoxicity"].astype(str).str.lower().isin(
            ["no", "non-hepatotoxic", "false", "0"])

    if "PAINS_Alerts_Swiss" in df.columns:
        filters["F6_PAINS_Clean"] = df["PAINS_Alerts_Swiss"].astype(str).isin(["0", "0.0"])

    for name, mask in filters.items():
        df[name] = mask

    if filters:
        all_pass = pd.DataFrame(filters).all(axis=1)
        df["Overall_Pass"] = all_pass
        df["Filters_Passed"] = pd.DataFrame(filters).sum(axis=1).astype(str) + "/" + str(len(filters))
    else:
        df["Overall_Pass"] = True
        df["Filters_Passed"] = "N/A"

    log.info(f"Filter results: {df['Overall_Pass'].sum()}/{len(df)} compounds passed all filters")
    return df


# ════════════════════════════════════════════════════════════
# STEP 4 — VISUALISATION (publication-quality)
# ════════════════════════════════════════════════════════════

def plot_bioavailability_radar(df: pd.DataFrame, output_dir: str) -> None:
    """
    Bioavailability radar chart for top 5 compounds (SwissADME-style).
    Properties: MW≤500, LogP≤5, HBD≤5, HBA≤10, TPSA≤140, RotBonds≤10
    Normalised 0–1 for radar display.
    """
    props = ["MW", "LogP_SwissADME", "HBD", "HBA", "TPSA", "Rotatable_Bonds"]
    labels = ["MW\n(≤500)", "LogP\n(≤5)", "HBD\n(≤5)", "HBA\n(≤10)", "TPSA\n(≤140)", "RotBonds\n(≤10)"]
    limits = [500, 5, 5, 10, 140, 10]

    df_plot = df[[c for c in props if c in df.columns]].dropna()
    if df_plot.empty:
        log.warning("No data available for radar plot (check column names match SwissADME output)")
        return

    top5 = df_plot.head(5)
    N = len([c for c in props if c in df.columns])
    if N < 3:
        return

    available_props = [p for p in props if p in df.columns]
    available_labels = [labels[i] for i, p in enumerate(props) if p in df.columns]
    available_limits = [limits[i] for i, p in enumerate(props) if p in df.columns]

    angles = np.linspace(0, 2*np.pi, len(available_props), endpoint=False).tolist()
    angles += angles[:1]

    fig, axes = plt.subplots(1, min(5, len(top5)), figsize=(4*min(5, len(top5)), 4),
                              subplot_kw=dict(polar=True))
    if len(top5) == 1:
        axes = [axes]

    colors = [BLUE, RED, GREEN, ORANGE, "#8E44AD"]
    for i, (idx, row) in enumerate(top5.iterrows()):
        ax   = axes[i]
        vals = [row[p]/lim for p, lim in zip(available_props, available_limits)]
        vals += vals[:1]
        ax.plot(angles, vals, color=colors[i], linewidth=2)
        ax.fill(angles, vals, color=colors[i], alpha=0.2)
        ax.set_thetagrids(np.degrees(angles[:-1]), available_labels, fontsize=8)
        ax.set_ylim(0, 1.2)
        ax.axhline(y=1.0, color="red", linestyle="--", linewidth=0.8, alpha=0.7)
        compound_name = df.iloc[idx]["Compound"] if "Compound" in df.columns else f"Cpd {i+1}"
        ax.set_title(str(compound_name)[:15], fontsize=9, fontweight="bold", pad=15, color=BLUE)

    plt.suptitle("Bioavailability Radar — Top Compounds\n(Red dashed = drug-like boundary)",
                 fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    out = os.path.join(output_dir, "bioavailability_radar.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"Radar chart saved: {out}")


def plot_admet_heatmap(df: pd.DataFrame, output_dir: str) -> None:
    """
    Categorical heatmap of ADMET pass/fail properties across all compounds.
    Green = Pass/Safe · Red = Fail/Unsafe · Grey = N/A
    """
    binary_cols = {
        "Lipinski_Pass":      "Lipinski Ro5",
        "GI_Absorption":      "GI Absorption (High)",
        "BBB_Permeant":       "BBB Permeant",
        "AMES_Mutagenicity":  "AMES Safe ✓",
        "hERG_Blockade":      "hERG Safe ✓",
        "Hepatotoxicity":     "Liver Safe ✓",
        "Skin_Sensitisation": "Skin Safe ✓",
        "Pgp_Substrate":      "P-gp Substrate",
        "CYP3A4_Inhibitor":   "CYP3A4 Inhibitor",
        "CYP2D6_Inhibitor":   "CYP2D6 Inhibitor",
    }
    available = {v: k for k, v in binary_cols.items() if k in df.columns}
    if not available:
        log.warning("No binary columns found for heatmap — check that pkCSM/SwissADME CSVs are correctly parsed")
        return

    safe_positive = {"Lipinski Ro5", "GI Absorption (High)", "BBB Permeant"}
    safe_negative  = {"AMES Safe ✓", "hERG Safe ✓", "Liver Safe ✓", "Skin Safe ✓"}

    def encode(col_name, series):
        """1=green (safe/pass), 0=red (fail/unsafe), 0.5=grey (unknown)"""
        result = []
        for v in series:
            s = str(v).strip().lower()
            if col_name in safe_negative:
                # For these, "No/False/0" means safe
                result.append(1 if s in ("no","non-inhibitor","non-mutagen","non-hepatotoxic",
                                          "non-sensitizer","false","0","0.0") else 0)
            else:
                result.append(1 if s in ("yes","pass","high","true","1","1.0") else 0)
        return result

    matrix_data = {}
    for label, col in available.items():
        matrix_data[label] = encode(label, df[col])

    matrix = pd.DataFrame(matrix_data, index=df["Compound"] if "Compound" in df.columns else range(len(df)))
    matrix = matrix.head(20)  # max 20 compounds for readability

    fig, ax = plt.subplots(figsize=(max(8, len(matrix.columns)*1.2), max(5, len(matrix)*0.45)))
    cmap = plt.cm.colors.LinearSegmentedColormap.from_list("admet", [RED, ORANGE, GREEN], N=2)
    sns.heatmap(matrix, ax=ax, cmap=cmap, vmin=0, vmax=1,
                linewidths=0.5, linecolor="white",
                cbar_kws={"ticks": [0.25, 0.75], "label": ""},
                annot=False)
    cbar = ax.collections[0].colorbar
    cbar.set_ticklabels(["Fail / Unsafe", "Pass / Safe"])
    ax.set_xlabel("", fontsize=11)
    ax.set_ylabel("Compound", fontsize=11)
    ax.set_title("Integrated ADMET Profile — pkCSM + SwissADME\n(Author: Sana Begum)",
                 fontsize=12, fontweight="bold", color=BLUE)
    plt.xticks(rotation=35, ha="right", fontsize=9)
    plt.yticks(rotation=0, fontsize=9)
    plt.tight_layout()
    out = os.path.join(output_dir, "admet_heatmap.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"ADMET heatmap saved: {out}")


def plot_toxicity_summary(df: pd.DataFrame, output_dir: str) -> None:
    """Bar chart showing percentage of compounds passing each safety filter."""
    tox_cols = {
        "AMES_Mutagenicity":  "AMES Non-mutagenic",
        "hERG_Blockade":      "hERG Non-blocker",
        "Hepatotoxicity":     "Non-hepatotoxic",
        "Skin_Sensitisation": "Non-sensitiser",
    }
    available = {k: v for k, v in tox_cols.items() if k in df.columns}
    if not available:
        return

    rates = {}
    for col, label in available.items():
        safe = df[col].astype(str).str.lower().isin(["no","false","0","0.0",
            "non-inhibitor","non-mutagen","non-hepatotoxic","non-sensitizer"])
        rates[label] = round(100 * safe.sum() / len(df), 1)

    fig, ax = plt.subplots(figsize=(8, 4))
    bars = ax.barh(list(rates.keys()), list(rates.values()),
                   color=[GREEN if v >= 70 else ORANGE if v >= 40 else RED for v in rates.values()],
                   edgecolor="white", height=0.5)
    for bar, val in zip(bars, rates.values()):
        ax.text(bar.get_width() + 1, bar.get_y() + bar.get_height()/2,
                f"{val}%", va="center", fontsize=10, fontweight="bold")
    ax.set_xlim(0, 115)
    ax.set_xlabel("% of Compounds Passing", fontsize=11)
    ax.set_title("Toxicity Safety Summary (pkCSM)", fontsize=12, fontweight="bold", color=BLUE)
    ax.axvline(x=70, color=BLUE, linestyle="--", alpha=0.5, label="70% threshold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = os.path.join(output_dir, "toxicity_summary.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    plt.close()
    log.info(f"Toxicity summary saved: {out}")


# ════════════════════════════════════════════════════════════
# REPORTING
# ════════════════════════════════════════════════════════════

def generate_report(df: pd.DataFrame, output_dir: str) -> None:
    """Save full integrated results and a filtered drug-like hits table."""
    all_path  = os.path.join(output_dir, "admet_integrated_full.csv")
    pass_path = os.path.join(output_dir, "admet_drug_like_hits.csv")
    df.to_csv(all_path, index=False)

    if "Overall_Pass" in df.columns:
        hits = df[df["Overall_Pass"] == True]
    else:
        hits = df

    hits.to_csv(pass_path, index=False)

    log.info(f"\n{'='*60}")
    log.info(f"  ADMET ANALYSIS COMPLETE (pkCSM + SwissADME)")
    log.info(f"{'='*60}")
    log.info(f"  Total compounds analysed : {len(df)}")
    log.info(f"  Drug-like hits (all pass): {len(hits)}")
    log.info(f"  Full results CSV         : {all_path}")
    log.info(f"  Filtered hits CSV        : {pass_path}")
    log.info(f"{'='*60}\n")


# ════════════════════════════════════════════════════════════
# MAIN
# ════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(
        description="Integrate pkCSM + SwissADME ADMET results",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Requirements")[0]
    )
    p.add_argument("--pkcsm",      help="pkCSM batch CSV output")
    p.add_argument("--swissadme",  help="SwissADME CSV output")
    p.add_argument("--smiles",     help="SDF or SMILES file (to prepare server inputs)")
    p.add_argument("--output",     default="admet_results/", help="Output directory")
    p.add_argument("--plot",       action="store_true", help="Generate all plots")
    p.add_argument("--prepare",    action="store_true",
                   help="Prepare SMILES for server submission only (use with --smiles)")
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    log.info("=" * 60)
    log.info("  ADMET Workflow: pkCSM + SwissADME  — Sana Begum")
    log.info("=" * 60)

    # Mode 1: Prepare inputs for web servers
    if args.prepare:
        if not args.smiles:
            log.error("Provide --smiles <file.sdf> with --prepare")
            sys.exit(1)
        out_smi = os.path.join(args.output, "compounds_for_servers.smi")
        prepare_smiles_from_sdf(args.smiles, out_smi)
        log.info(f"\nNext steps:\n"
                 f"  1. pkCSM  → https://biosig.lab.uq.edu.au/pkcsm/prediction\n"
                 f"             Upload {out_smi} → Download CSV\n"
                 f"  2. SwissADME → https://www.swissadme.ch/\n"
                 f"             Paste SMILES list → Export CSV\n"
                 f"  3. Re-run: python admet_pkcsm_swissadme.py "
                 f"--pkcsm pkcsm.csv --swissadme swissadme.csv --output {args.output} --plot")
        return

    # Mode 2: Integrate results
    if not args.pkcsm and not args.swissadme:
        log.error("Provide at least one of --pkcsm or --swissadme CSV result files.")
        log.info("Use --prepare to generate server input files from your SDF/SMILES first.")
        sys.exit(1)

    dfs = []
    pkcsm_df, swiss_df = pd.DataFrame(), pd.DataFrame()
    if args.pkcsm:
        pkcsm_df = parse_pkcsm(args.pkcsm)
        dfs.append(pkcsm_df)
    if args.swissadme:
        swiss_df = parse_swissadme(args.swissadme)
        dfs.append(swiss_df)

    if pkcsm_df.empty and not swiss_df.empty:
        merged = swiss_df
    elif swiss_df.empty and not pkcsm_df.empty:
        merged = pkcsm_df
    else:
        merged = integrate_admet(pkcsm_df, swiss_df)

    merged = apply_drug_filters(merged)
    generate_report(merged, args.output)

    if args.plot:
        plot_admet_heatmap(merged, args.output)
        plot_bioavailability_radar(merged, args.output)
        plot_toxicity_summary(merged, args.output)
        log.info("All plots generated.")


if __name__ == "__main__":
    main()
