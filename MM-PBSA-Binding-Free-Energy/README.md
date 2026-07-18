# Project 5 — MM/PBSA Binding Free Energy Analysis

## Overview

Complete workflow for calculating protein–ligand binding free energy using the Molecular Mechanics/Poisson–Boltzmann Surface Area (MM-PBSA) method. The pipeline integrates molecular dynamics trajectories from GROMACS with MM-PBSA calculations to quantify binding affinity, decompose energetic contributions by residue, and identify key intermolecular interactions driving ligand binding.

This workflow has been implemented in published computational drug discovery research on SARS-CoV-2 Main Protease and DNA methyltransferase inhibitors, providing quantitative validation of docking predictions and enabling rigorous ranking of lead compounds.

---

## Workflow Architecture

```
Completed MD Trajectory (md.xtc from GROMACS)
    ↓
Trajectory Processing & Frame Selection
    ↓
Index File Generation (Protein · Ligand · System)
    ↓
MM-PBSA Parameter Configuration
    ↓
Energy Calculation Phase
    ├─ Molecular Mechanics (MM) Energy
    ├─ Poisson-Boltzmann (PB) Solvation
    ├─ Non-Polar Solvation (SASA)
    └─ Entropic Contribution (optional)
    ↓
Per-Frame Free Energy Calculation
    ↓
Trajectory-Averaged Binding Free Energy (ΔG_bind)
    ↓
Residue Decomposition Analysis
    ↓
Publication-Ready ΔG Estimation & Visualization
```

---

## Thermodynamic Basis

### Binding Free Energy Equation

**ΔG_bind = ΔE_MM + ΔG_solvation − TΔS_conf**

Where:
- **ΔE_MM** = Molecular mechanics energy (bonded + non-bonded)
- **ΔG_solvation** = ΔG_PB + ΔG_SASA (electrostatic + non-polar)
- **TΔS_conf** = Entropic contribution from conformational changes

---

## Software & Tools

| Tool | Purpose |
|------|---------|
| **g_mmpbsa** | GROMACS interface for MM-PBSA calculations |
| **MMPBSA.py (AmberTools)** | Alternative MM-PBSA implementation |
| **GROMACS** | Trajectory processing and energy extraction |
| **Poisson-Boltzmann Solver** | Solvation energy estimation |
| **Python** | Data analysis and visualization |

---

## Key Features

✅ Thermodynamically rigorous binding free energy calculation  
✅ Integration of MM, PB, and SASA energy components  
✅ Per-frame energy decomposition and averaging  
✅ Residue-level contribution analysis  
✅ Hot-spot identification for structure-based design  
✅ Convergence diagnostics and statistical analysis  
✅ Publication-quality energy decomposition visualizations  

---

## Interpretation

### Binding Free Energy Scale

- **Excellent Binders:** ΔG_bind < −10 kcal/mol  
- **Good Binders:** ΔG_bind −7 to −10 kcal/mol  
- **Moderate Binders:** ΔG_bind −5 to −7 kcal/mol  
- **Weak Binders:** ΔG_bind > −5 kcal/mol  

### Energy Decomposition

**Favorable Components:**
- ΔE_MM: Enthalpic gain from interactions
- ΔG_PB: Electrostatic stabilization
- ΔG_SASA: Hydrophobic burial

**Hot Spots:** Residues with ΔG < −1.0 kcal/mol represent optimal mutation targets

---

## Output Files

- `binding_free_energy.csv` — **Final ΔG_bind values (per frame and averaged)**
- `decomp_results.dat` — Per-residue energy decomposition
- `binding_energy_plot.png` — ΔG convergence over trajectory
- `energy_components_bar.png` — Stacked bar: MM + PB + SASA
- `decomp_heatmap.png` — Per-residue contribution heatmap
- `hot_spots.csv` — Identified binding hot-spot residues
- `mmpbsa_summary.txt` — **Publication-ready ΔG_bind summary**

---

## Applications

- Lead Ranking by calculated binding affinity
- Binding Mode Validation
- Structure-Based Optimization
- Affinity Prediction and Correlation
- Multi-target Profiling
- Drug-Target Selectivity Assessment

---

## Research Relevance

Implements rigorous MM-PBSA methodology to calculate binding free energies from MD trajectories, providing quantitative validation of docking predictions and enabling evidence-based prioritization of lead compounds. Residue decomposition analysis enables structure-guided optimization by identifying key interaction anchors and potential mutation sites.

---

## Best Practices

✓ Use consistent force fields throughout  
✓ Ensure adequate MD equilibration (100–200 ns minimum)  
✓ Decimate frames appropriately for efficiency  
✓ Validate calculations against experimental binding data  
✓ Report confidence intervals from trajectory statistics  

---

## References

- Kollman et al. (2000). "Calculating structures and free energies of complex molecules"
- MMPBSA.py: https://ambermd.org/
- g_mmpbsa: https://www.gromacs.org/
