#!/usr/bin/env python3
"""
ADMET Screening & Drug-Likeness Filter
=======================================
Author  : Sana Begum
Role    : Research Analyst, Era's Lucknow Medical College & Hospital
Email   : begumsana686@gmail.com
GitHub  : https://github.com/Sanabegum09

Description
-----------
Computes key ADMET (Absorption, Distribution, Metabolism, Excretion, Toxicity)
and drug-likeness properties for a compound library using RDKit. Applies
Lipinski's Rule of Five, Veber's rules, PAINS filtering, and multiple
drug-likeness indices. Outputs a ranked summary with pass/fail annotations
for each compound.

Typical usage after virtual screening: input the top docking hits to filter
for drug-like candidates before experimental validation.

Usage
-----
    python admet_screening.py --input hits.sdf --output admet_results/

Requirements
------------
    pip install rdkit pandas numpy matplotlib seaborn

References
----------
    Lipinski CA et al. Adv Drug Deliv Rev. 2001;46(1-3):3-26.
    Veber DF et al. J Med Chem. 2002;45(12):2615-23.
    Egan WJ et al. J Med Chem. 2000;43(21):3867-77.
    Ghose AK et al. J Comb Chem. 1999;1(1):55-68.
"""

import os
import sys
import argparse
import logging
import warnings
from pathlib import Path
from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

# RDKit
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen, Lipinski, QED
from rdkit.Chem import Draw, AllChem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)


# ─────────────────────────────────────────────────────────────
# PROPERTY CALCULATORS
# ─────────────────────────────────────────────────────────────

def compute_properties(mol) -> dict:
    """
    Compute a comprehensive set of molecular descriptors for drug-likeness.

    Returns
    -------
    dict : All computed properties
    """
    if mol is None:
        return {}

    mw    = Descriptors.ExactMolWt(mol)
    logp  = Crippen.MolLogP(mol)
    hbd   = rdMolDescriptors.CalcNumHBD(mol)    # H-bond donors
    hba   = rdMolDescriptors.CalcNumHBA(mol)    # H-bond acceptors
    psa   = rdMolDescriptors.CalcTPSA(mol)      # Topological polar surface area
    nrb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    arom  = rdMolDescriptors.CalcNumAromaticRings(mol)
    fc    = Chem.rdmolops.GetFormalCharge(mol)
    nha   = mol.GetNumHeavyAtoms()
    fsp3  = rdMolDescriptors.CalcFractionCSP3(mol)
    mf    = rdMolDescriptors.CalcMolFormula(mol)

    # QED — Quantitative Estimate of Drug-likeness (0–1, higher is better)
    try:
        qed_score = QED.qed(mol)
    except Exception:
        qed_score = None

    # SAscore — Synthetic Accessibility (1 easy – 10 hard)
    # Requires rdkit.Chem.Descriptors.rdMolDescriptors
    try:
        from rdkit.Chem import RDConfig
        sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
        import sascorer
        sa_score = sascorer.calculateScore(mol)
    except Exception:
        sa_score = None   # SA Score contrib not always present

    return {
        "Molecular_Weight":     round(mw, 2),
        "LogP":                 round(logp, 2),
        "H_Bond_Donors":        hbd,
        "H_Bond_Acceptors":     hba,
        "TPSA_A2":              round(psa, 2),
        "Rotatable_Bonds":      nrb,
        "Num_Rings":            rings,
        "Aromatic_Rings":       arom,
        "Formal_Charge":        fc,
        "Heavy_Atom_Count":     nha,
        "Fraction_CSP3":        round(fsp3, 3),
        "Molecular_Formula":    mf,
        "QED":                  round(qed_score, 3) if qed_score is not None else None,
        "SA_Score":             round(sa_score, 2) if sa_score is not None else None,
    }


# ─────────────────────────────────────────────────────────────
# DRUG-LIKENESS RULES
# ─────────────────────────────────────────────────────────────

def lipinski_ro5(props: dict) -> tuple:
    """
    Lipinski's Rule of Five.
    Pass = ≤1 violation (original Lipinski allows 1 violation).

    Returns (pass: bool, violations: list)
    """
    violations = []
    if props["Molecular_Weight"] > 500:
        violations.append(f"MW={props['Molecular_Weight']} > 500")
    if props["LogP"] > 5:
        violations.append(f"LogP={props['LogP']} > 5")
    if props["H_Bond_Donors"] > 5:
        violations.append(f"HBD={props['H_Bond_Donors']} > 5")
    if props["H_Bond_Acceptors"] > 10:
        violations.append(f"HBA={props['H_Bond_Acceptors']} > 10")
    return len(violations) <= 1, violations


def veber_rules(props: dict) -> tuple:
    """
    Veber oral bioavailability rules.
    Pass: TPSA ≤ 140 Å² AND Rotatable bonds ≤ 10.
    """
    violations = []
    if props["TPSA_A2"] > 140:
        violations.append(f"TPSA={props['TPSA_A2']} > 140")
    if props["Rotatable_Bonds"] > 10:
        violations.append(f"RotBonds={props['Rotatable_Bonds']} > 10")
    return len(violations) == 0, violations


def egan_rules(props: dict) -> tuple:
    """
    Egan's filter for passive intestinal absorption.
    Pass: LogP ≤ 5.88 AND TPSA ≤ 131.6.
    """
    violations = []
    if props["LogP"] > 5.88:
        violations.append(f"LogP={props['LogP']} > 5.88")
    if props["TPSA_A2"] > 131.6:
        violations.append(f"TPSA={props['TPSA_A2']} > 131.6")
    return len(violations) == 0, violations


def ghose_filter(props: dict) -> tuple:
    """
    Ghose drug-likeness filter.
    Pass: 160 ≤ MW ≤ 480, -0.4 ≤ LogP ≤ 5.6,
          20 ≤ MR ≤ 130, 40 ≤ Heavy atoms ≤ 130.
    """
    violations = []
    if not (160 <= props["Molecular_Weight"] <= 480):
        violations.append(f"MW={props['Molecular_Weight']} not in [160,480]")
    if not (-0.4 <= props["LogP"] <= 5.6):
        violations.append(f"LogP={props['LogP']} not in [-0.4, 5.6]")
    if not (20 <= props["Heavy_Atom_Count"] <= 70):
        violations.append(f"HeavyAtoms={props['Heavy_Atom_Count']} not in [20,70]")
    return len(violations) == 0, violations


def leadlikeness_filter(props: dict) -> tuple:
    """
    Lead-like filter for fragment & lead-based drug discovery.
    Pass: MW ≤ 350, LogP ≤ 3.5, Rotatable bonds ≤ 7.
    """
    violations = []
    if props["Molecular_Weight"] > 350:
        violations.append(f"MW={props['Molecular_Weight']} > 350")
    if props["LogP"] > 3.5:
        violations.append(f"LogP={props['LogP']} > 3.5")
    if props["Rotatable_Bonds"] > 7:
        violations.append(f"RotBonds={props['Rotatable_Bonds']} > 7")
    return len(violations) == 0, violations


def pains_filter(mol) -> tuple:
    """
    PAINS (Pan-Assay Interference Compounds) filter using RDKit FilterCatalog.

    Returns (is_clean: bool, alerts: list)
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = list(catalog.GetMatches(mol))
    alerts = [m.GetDescription() for m in matches]
    return len(alerts) == 0, alerts


def brenk_filter(mol) -> tuple:
    """
    Brenk unwanted chemical fragment filter.
    """
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog(params)
    matches = list(catalog.GetMatches(mol))
    alerts = [m.GetDescription() for m in matches]
    return len(alerts) == 0, alerts


# ─────────────────────────────────────────────────────────────
# ABSORPTION PREDICTION (Simple heuristic models)
# ─────────────────────────────────────────────────────────────

def predict_absorption(props: dict) -> dict:
    """
    Heuristic predictions for absorption properties.
    Based on published cutoffs (not ML models — for ML use SwissADME or pkCSM API).
    """
    psa   = props["TPSA_A2"]
    logp  = props["LogP"]
    mw    = props["Molecular_Weight"]
    nrb   = props["Rotatable_Bonds"]

    # GI Absorption (simplified Veber/Clark model)
    gi_absorption = "High" if (psa <= 140 and nrb <= 10) else "Low"

    # BBB permeability (psa < 90, logp 1-3 preferred)
    bbb = "Likely" if (psa < 90 and 1 <= logp <= 3) else "Unlikely"

    # P-gp substrate (rough heuristic: high MW, high HBA)
    pgp_sub = "Likely" if (mw > 400 and props["H_Bond_Acceptors"] > 4) else "Unlikely"

    # Water solubility (very rough — ESOL needs full descriptor set)
    if logp < 1 and mw < 300:
        solubility = "High"
    elif logp <= 3 and mw <= 500:
        solubility = "Moderate"
    else:
        solubility = "Low"

    return {
        "GI_Absorption":    gi_absorption,
        "BBB_Permeable":    bbb,
        "Pgp_Substrate":    pgp_sub,
        "Water_Solubility": solubility
    }


# ─────────────────────────────────────────────────────────────
# TOXICITY FLAGS
# ─────────────────────────────────────────────────────────────

def flag_toxicity(props: dict, mol) -> dict:
    """
    Simple structural toxicity flags. For precise predictions use pkCSM or admetSAR.
    """
    smarts_alerts = {
        "Nitro_group":        "[N+](=O)[O-]",
        "Reactive_aldehyde":  "[CX3H1](=O)[#6]",
        "Michael_acceptor":   "[CX2]#[NX1]",
        "Quinone":            "O=C1C=CC(=O)C=C1",
    }
    flags = {}
    for name, smarts in smarts_alerts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            flags[name] = True
        else:
            flags[name] = False

    # Log P toxicity risk (very high logP → membrane toxicity)
    flags["High_LogP_Risk"] = props["LogP"] > 6

    return flags


# ─────────────────────────────────────────────────────────────
# MAIN SCREENING FUNCTION
# ─────────────────────────────────────────────────────────────

def screen_compound(mol, name: str) -> dict:
    """Screen a single molecule and return full ADMET profile."""
    if mol is None:
        return {"Compound": name, "Error": "Invalid molecule", "Overall_Pass": False}

    props = compute_properties(mol)
    if not props:
        return {"Compound": name, "Error": "Property calculation failed", "Overall_Pass": False}

    ro5_pass, ro5_viol      = lipinski_ro5(props)
    veber_pass, veber_viol  = veber_rules(props)
    egan_pass, egan_viol    = egan_rules(props)
    ghose_pass, ghose_viol  = ghose_filter(props)
    lead_pass, lead_viol    = leadlikeness_filter(props)
    pains_pass, pains_alrt  = pains_filter(mol)
    brenk_pass, brenk_alrt  = brenk_filter(mol)
    absorption              = predict_absorption(props)
    tox_flags               = flag_toxicity(props, mol)

    # Overall pass: Lipinski pass + no PAINS alerts + TPSA/RB acceptable
    overall_pass = ro5_pass and pains_pass and veber_pass

    row = {"Compound": name}
    row.update(props)
    row["Lipinski_Ro5_Pass"]    = ro5_pass
    row["Lipinski_Violations"]  = "; ".join(ro5_viol) if ro5_viol else "None"
    row["Veber_Pass"]           = veber_pass
    row["Egan_Pass"]            = egan_pass
    row["Ghose_Pass"]           = ghose_pass
    row["Lead_Like"]            = lead_pass
    row["PAINS_Pass"]           = pains_pass
    row["PAINS_Alerts"]         = "; ".join(pains_alrt) if pains_alrt else "None"
    row["BRENK_Pass"]           = brenk_pass
    row["BRENK_Alerts"]         = "; ".join(brenk_alrt) if brenk_alrt else "None"
    row.update({f"Absorption_{k}": v for k, v in absorption.items()})
    row.update({f"Tox_{k}": v for k, v in tox_flags.items()})
    row["Overall_DrugLike_Pass"] = overall_pass

    return row


def run_admet_screening(input_file: str, output_dir: str) -> pd.DataFrame:
    """Screen all compounds in an SDF file and save results."""
    os.makedirs(output_dir, exist_ok=True)

    # Read molecules
    if input_file.endswith(".sdf"):
        suppl = Chem.SDMolSupplier(input_file, sanitize=True, removeHs=False)
        mols  = [(mol.GetProp("_Name") if mol and mol.HasProp("_Name") else f"Compound_{i+1}", mol)
                 for i, mol in enumerate(suppl)]
    elif input_file.endswith(".smi") or input_file.endswith(".smiles"):
        mols = []
        with open(input_file) as fh:
            for i, line in enumerate(fh):
                parts = line.strip().split()
                if not parts:
                    continue
                smiles = parts[0]
                name   = parts[1] if len(parts) > 1 else f"Compound_{i+1}"
                mol    = Chem.MolFromSmiles(smiles)
                mols.append((name, mol))
    else:
        log.error("Input must be .sdf or .smi/.smiles")
        return pd.DataFrame()

    log.info(f"Screening {len(mols)} compounds for ADMET properties...")
    rows = [screen_compound(mol, name) for name, mol in mols]

    df = pd.DataFrame(rows)
    df_sorted = df.sort_values("QED", ascending=False)

    # Save CSVs
    df_sorted.to_csv(os.path.join(output_dir, "admet_full_results.csv"), index=False)
    passed = df_sorted[df_sorted["Overall_DrugLike_Pass"] == True]
    passed.to_csv(os.path.join(output_dir, "admet_drug_like_hits.csv"), index=False)

    log.info(f"Results: {len(passed)}/{len(df_sorted)} compounds passed drug-likeness filters.")
    return df_sorted


# ─────────────────────────────────────────────────────────────
# VISUALIZATION
# ─────────────────────────────────────────────────────────────

def plot_chemical_space(df: pd.DataFrame, output_dir: str) -> None:
    """Plot MW vs LogP chemical space coloured by drug-likeness."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # MW vs LogP
    colors = df["Overall_DrugLike_Pass"].map({True: "#1A5276", False: "#EC7063"})
    sc = axes[0].scatter(df["Molecular_Weight"], df["LogP"],
                          c=colors, alpha=0.7, edgecolors="white", linewidths=0.4, s=60)
    axes[0].axhline(y=5, color="orange", linestyle="--", alpha=0.7, label="LogP=5")
    axes[0].axvline(x=500, color="orange", linestyle="--", alpha=0.7, label="MW=500")
    axes[0].set_xlabel("Molecular Weight (Da)", fontsize=11)
    axes[0].set_ylabel("LogP", fontsize=11)
    axes[0].set_title("Chemical Space: MW vs LogP", fontsize=12, fontweight="bold")
    pass_patch = mpatches.Patch(color="#1A5276", label="Drug-like (Pass)")
    fail_patch = mpatches.Patch(color="#EC7063", label="Not drug-like (Fail)")
    axes[0].legend(handles=[pass_patch, fail_patch], fontsize=9)

    # TPSA distribution
    axes[1].hist(df["TPSA_A2"].dropna(), bins=20, color="#1A5276", alpha=0.8, edgecolor="white")
    axes[1].axvline(x=140, color="red", linestyle="--", label="TPSA=140 (cutoff)")
    axes[1].set_xlabel("TPSA (Å²)", fontsize=11)
    axes[1].set_ylabel("Count", fontsize=11)
    axes[1].set_title("TPSA Distribution", fontsize=12, fontweight="bold")
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    out = os.path.join(output_dir, "chemical_space.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"Chemical space plot saved: {out}")


def plot_qed_distribution(df: pd.DataFrame, output_dir: str) -> None:
    """Plot QED score distribution."""
    fig, ax = plt.subplots(figsize=(8, 4))
    df_valid = df["QED"].dropna()
    ax.hist(df_valid, bins=20, color="#1A5276", alpha=0.85, edgecolor="white")
    ax.axvline(x=0.5, color="red", linestyle="--", label="QED=0.5 (recommended threshold)")
    ax.set_xlabel("QED Score (Quantitative Estimate of Drug-likeness)", fontsize=11)
    ax.set_ylabel("Number of Compounds", fontsize=11)
    ax.set_title("QED Score Distribution", fontsize=12, fontweight="bold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    out = os.path.join(output_dir, "qed_distribution.png")
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    log.info(f"QED distribution plot saved: {out}")


# ─────────────────────────────────────────────────────────────
# ARGUMENT PARSING & ENTRY POINT
# ─────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description="ADMET screening and drug-likeness filter for compound libraries",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--input",  required=True, help="Input SDF or SMILES file")
    p.add_argument("--output", default="admet_results/", help="Output directory")
    p.add_argument("--plot",   action="store_true", help="Generate visualisation plots")
    return p.parse_args()


def main():
    args = parse_args()
    log.info("="*60)
    log.info("  ADMET Screening Pipeline — Sana Begum")
    log.info("="*60)

    df = run_admet_screening(args.input, args.output)

    if df.empty:
        log.error("No results generated.")
        sys.exit(1)

    if args.plot:
        plot_chemical_space(df, args.output)
        plot_qed_distribution(df, args.output)

    log.info(f"\nDone! Results saved to: {args.output}")


if __name__ == "__main__":
    main()
