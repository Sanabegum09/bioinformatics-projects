# Project 1 — Molecular Docking & Virtual Screening Pipeline

## Overview

This project presents a complete and reproducible Structure-Based Drug Discovery (SBDD) workflow for identifying potential lead compounds against therapeutic protein targets. The pipeline integrates protein preparation, ligand library generation, molecular docking, virtual screening, interaction analysis, and lead prioritization using AutoDock Vina, Discovery Studio Visualizer, and PyMOL.

This workflow has been developed and validated through computational strategies applied in published research on phytochemical inhibitors against SARS-CoV-2 Main Protease and DNA methylotransferases (DNMTs) for head and neck cancer, demonstrating efficacy in identifying high-affinity lead compounds prior to molecular dynamics validation and MM-PBSA binding free-energy calculations.

---

## Workflow Architecture

```
Target Protein (PDB) → Protein Preparation → Ligand Library → Energy Minimization 
    ↓
PDBQT Conversion → Grid Box Definition → AutoDock Vina Screening → Docking Ranking
    ↓
Protein–Ligand Interaction Analysis → Lead Prioritization → Publication-ready Output
```

---

## Software & Tools

| Tool | Purpose |
|------|---------|
| **AutoDock Vina** | Molecular docking and virtual screening |
| **AutoDock Tools (ADT)** | Protein and ligand PDBQT preparation |
| **Discovery Studio Visualizer** | Protein-ligand interaction analysis |
| **PyMOL** | Structure visualization and analysis |
| **UCSF Chimera** | Protein preparation and refinement |
| **Open Babel** | File format conversion and energy minimization |

---

## Key Features

✅ End-to-end molecular docking workflow  
✅ Automated protein and ligand preparation  
✅ High-throughput virtual screening capability  
✅ Comprehensive protein-ligand interaction analysis  
✅ Publication-quality visualization output  
✅ Seamless integration with ADMET, MD, and MM-PBSA pipelines  

---

## Applications

- Computer-Aided Drug Design (CADD)
- Structure-Based Drug Discovery
- Virtual Screening
- Natural Product Screening
- Antiviral Drug Discovery
- Cancer Drug Discovery

---

## Research Relevance

Reflects validated methodologies from published drug discovery studies targeting SARS-CoV-2 Main Protease and DNA methyltransferase inhibitors for cancer therapy.

---

## Next Steps

1. **Project 2** — ADMET Prediction Pipeline
2. **Project 4** — Molecular Dynamics Simulation (GROMACS)
3. **Project 5** — MM-PBSA Binding Free Energy Analysis
