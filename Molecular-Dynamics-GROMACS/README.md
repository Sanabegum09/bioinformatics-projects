# Project 4 — Molecular Dynamics Simulation using GROMACS

## Overview

Complete and reproducible Molecular Dynamics (MD) Simulation workflow for investigating structural stability, conformational dynamics, and temporal interaction behavior of protein–ligand complexes under physiological conditions. The pipeline encompasses system preparation, energy minimization, NVT/NPT equilibration, production simulation, and comprehensive trajectory analysis.

This workflow has been validated through computational protocols in published drug discovery research on SARS-CoV-2 Main Protease (Mpro) and DNA Methyltransferase (DNMT) inhibitors, providing critical validation of docking predictions prior to binding free-energy calculations.

---

## Workflow Architecture

```
Protein–Ligand Complex (from Docking)
    ↓
Topology & Force Field Generation
    ↓
Ligand Parameterization (GAFF/CGenFF)
    ↓
Simulation Box Definition & Solvation (TIP3P Water)
    ↓
Ionic Strength Neutralization (Na+/Cl− ions)
    ↓
Steepest Descent Energy Minimization
    ↓
NVT Equilibration (Temperature equilibration)
    ↓
NPT Equilibration (Pressure & Temperature)
    ↓
Production MD Simulation (100–200 ns)
    ↓
Trajectory Analysis (RMSD · RMSF · Rg · SASA · H-bonds)
    ↓
Publication-quality Analysis Plots & Stability Assessment
```

---

## Software & Tools

| Tool | Purpose |
|------|---------|
| **GROMACS 2020+** | Molecular dynamics simulation engine |
| **SwissParam** | Ligand topology/parameter generation |
| **CGenFF** | Ligand parameterization (CHARMM force field) |
| **PyMOL** | Structure visualization and inspection |
| **Python (MDAnalysis)** | Trajectory analysis and data processing |

---

## Key Features

✅ Complete MD simulation workflow from preparation to analysis  
✅ Automated force field parameterization  
✅ Physiological conditions (298 K, 1 bar, 0.15 M ionic strength)  
✅ Energy minimization with convergence verification  
✅ Equilibration protocols (NVT → NPT)  
✅ Long-timescale production simulations (100–200 ns)  
✅ Comprehensive trajectory analysis suite  
✅ Integration with ADMET and MM-PBSA workflows  

---

## Trajectory Analysis

**Metrics Calculated:**
- Root Mean Square Deviation (RMSD) — Structural stability
- Root Mean Square Fluctuation (RMSF) — Per-residue flexibility
- Radius of Gyration (Rg) — Protein compactness
- Solvent Accessible Surface Area (SASA) — Surface exposure
- Hydrogen Bonds — Interaction persistence
- Secondary Structure (DSSP) — α-helix/β-sheet stability
- Principal Component Analysis (PCA) — Conformational motions

---

## Output Files

- `em.gro` — Energy-minimized structure
- `md.xtc` — Complete MD trajectory (100–200 ns)
- `rmsd.png` — Structural deviation over time
- `rmsf.png` — Per-residue flexibility profile
- `hbonds.png` — Protein-ligand hydrogen bond stability
- `trajectory_analysis.csv` — Summary statistics

---

## Applications

- Lead Compound Validation
- Protein Stability Analysis
- Binding Mode Assessment
- Drug Discovery
- Biomolecular Simulations
- Molecular Mechanism Studies

---

## Research Relevance

Demonstrates expertise in molecular dynamics simulations for validating structure-based drug discovery results, enabling informed prioritization of candidates for binding free-energy calculations and experimental validation.

---

## Next Steps

1. **Project 5** — MM-PBSA Binding Free Energy Analysis
2. Experimental Validation
3. Multi-target Screening

---

## References

- GROMACS Manual: https://manual.gromacs.org/
- SwissParam: https://www.swissparam.ch/
- CGenFF: https://cgenff.paramchem.org/
