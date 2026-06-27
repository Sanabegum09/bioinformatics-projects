# 💊 ADMET Screening & Drug-Likeness Filter

Computes ADMET properties and drug-likeness scores for compound libraries using **RDKit**. Use this after virtual screening to prioritize hits for experimental validation.

## Properties Computed
| Category | Properties |
|----------|-----------|
| Physicochemical | MW, LogP, HBD, HBA, TPSA, Rotatable bonds, Rings, Formal charge, Fsp3 |
| Drug-likeness | Lipinski Ro5, Veber rules, Egan filter, Ghose filter, Lead-likeness |
| Quality | QED score (0–1), SA Score (1–10) |
| Flags | PAINS alerts, BRENK unwanted fragments |
| Absorption | GI absorption, BBB permeability, P-gp substrate, water solubility |
| Toxicity | Nitro groups, Michael acceptors, reactive aldehydes, high-LogP risk |

## Quick Start
```bash
pip install -r requirements.txt

# Screen an SDF file
python admet_screening.py --input hits.sdf --output admet_results/ --plot

# Screen a SMILES file
python admet_screening.py --input compounds.smi --output admet_results/ --plot
```

## Output
| File | Description |
|------|-------------|
| `admet_results/admet_full_results.csv` | All compounds with every property |
| `admet_results/admet_drug_like_hits.csv` | Filtered drug-like candidates |
| `admet_results/chemical_space.png` | MW vs LogP coloured by drug-likeness |
| `admet_results/qed_distribution.png` | QED score distribution |

---
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College | [GitHub](https://github.com/Sanabegum09)
