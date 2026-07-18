# Project 2 — ADMET Prediction & Drug-likeness Assessment Pipeline

## Overview

Comprehensive pharmacokinetic and toxicological assessment workflow for evaluating drug-likeness potential of lead compounds. The pipeline integrates automated ADMET screening using pkCSM and SwissADME servers, generating predictions for absorption, distribution, metabolism, excretion, and toxicity properties essential for prioritizing candidates.

---

## Workflow Architecture

```
Top-Ranked Docking Hits (SDF/SMILES) 
    ↓
SMILES/MOL2 Format Preparation & Cleaning
    ↓
Batch Submission to pkCSM Server (ADMET Properties)
    ↓
Batch Submission to SwissADME Server (Drug-likeness)
    ↓
Automated Results Parsing & Integration
    ↓
Publication-ready ADMET Report & Drug-like Hit Selection
```

---

## Web Servers & Predictions

### pkCSM Server
**URL:** https://biosig.lab.uq.edu.au/pkcsm/

| Prediction Category | Properties Assessed |
|---|---|
| **Absorption** | Caco-2 permeability · HIA · P-glycoprotein substrate/inhibitor |
| **Distribution** | VDss · BBB permeability · CNS permeability |
| **Metabolism** | CYP substrate/inhibition panel |
| **Excretion** | Renal OCT2 substrate · Total clearance |
| **Toxicity** | AMES mutagenicity · hERG blockade · Hepatotoxicity · LD50 · Skin sensitization |

### SwissADME Server
**URL:** https://www.swissadme.ch/

| Prediction Category | Properties Assessed |
|---|---|
| **Physicochemical** | MW · LogP · H-bond donors/acceptors · Rotatable bonds |
| **Drug-likeness** | Lipinski · Veber · Egan · Ghose rules · QED |
| **Bioavailability** | GI absorption · BBB permeant · PAINS/Brenk alerts |

---

## Key Features

✅ Automated batch submission to industry-standard servers  
✅ Comprehensive ADMET predictions (A·D·M·E·T)  
✅ Drug-likeness assessment (Lipinski, Veber, Egan, Ghose)  
✅ Integrated bioavailability prediction  
✅ Publication-quality heatmap visualizations  
✅ Systematic lead compound prioritization  

---

## Applications

- Lead Optimization
- Drug Safety Assessment
- Bioavailability Prediction
- ADMET Filtering
- Multi-target Screening
- Translational Drug Discovery

---

## Output Files

- `pkcsm_results.csv` — Raw ADMET predictions
- `swissadme_results.csv` — Physicochemical properties
- `merged_admet.csv` — Integrated profile
- `admet_heatmap.png` — ADMET property visualization
- `admet_drug_like_hits.csv` — **Filtered drug-like candidates**

---

## Research Relevance

Implements established ADMET assessment methodologies enabling systematic evaluation of pharmacokinetic properties essential for lead optimization.

---

## References

- pkCSM: https://biosig.lab.uq.edu.au/pkcsm/
- SwissADME: https://www.swissadme.ch/
