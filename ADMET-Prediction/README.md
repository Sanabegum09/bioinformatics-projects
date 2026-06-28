# 💊 ADMET Analysis: pkCSM + SwissADME Workflow

**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College & Hospital  
**Email:** begumsana686@gmail.com | **GitHub:** [@Sanabegum09](https://github.com/Sanabegum09)

This workflow replicates the ADMET analysis used in Sana Begum's published research:
> *"Phytochemicals as DNMTs inhibitors for targeted therapy in head and neck cancer: An in-silico study."*
> Results in Chemistry 18 (2025) 102783. [doi:10.1016/j.rechem.2025.102783](https://doi.org/10.1016/j.rechem.2025.102783)

---

## Servers Used

| Server | URL | Predicts |
|--------|-----|---------|
| **pkCSM** | https://biosig.lab.uq.edu.au/pkcsm/ | Absorption, Distribution, Metabolism, Excretion, Toxicity |
| **SwissADME** | https://www.swissadme.ch/ | Physicochemical, Lipinski/Veber/Egan rules, Bioavailability radar |

---

## Full ADMET Workflow

### Step 1 — Prepare compound SMILES
```bash
python admet_pkcsm_swissadme.py \
    --smiles  top_docking_hits.sdf \
    --output  admet_results/ \
    --prepare
# → generates admet_results/compounds_for_servers.smi
```

### Step 2 — Submit to pkCSM
1. Go to https://biosig.lab.uq.edu.au/pkcsm/prediction
2. Click **Batch** tab
3. Paste the SMILES from `compounds_for_servers.smi`
4. Click **Predict** → wait → **Download results as CSV**
5. Save as `pkcsm_results.csv`

### Step 3 — Submit to SwissADME
1. Go to https://www.swissadme.ch/
2. Paste the SMILES list (one per line, without names)
3. Click **Run** → **Export** → Download CSV
4. Save as `swissadme_results.csv`

### Step 4 — Integrate & visualise
```bash
python admet_pkcsm_swissadme.py \
    --pkcsm     pkcsm_results.csv \
    --swissadme swissadme_results.csv \
    --output    admet_results/ \
    --plot
```

---

## Properties Predicted

### pkCSM
| Category | Properties |
|----------|-----------|
| **Absorption** | Caco-2 permeability, HIA (%), P-gp substrate/inhibitor |
| **Distribution** | VDss, BBB permeability, CNS permeability, PPB |
| **Metabolism** | CYP1A2/2C9/2C19/2D6/3A4 inhibition & substrate |
| **Excretion** | Renal OCT2 substrate, Total clearance |
| **Toxicity** | AMES mutagenicity, hERG I/II, Hepatotoxicity, LD50, Skin sensitisation |

### SwissADME
| Category | Properties |
|----------|-----------|
| **Physicochemical** | MW, LogP (iLOGP/XLOGP3/WLOGP), HBD, HBA, TPSA, Rotatable bonds |
| **Drug-likeness** | Lipinski Ro5, Veber, Egan, Ghose, Lead-likeness |
| **Absorption** | GI absorption (High/Low), BBB permeant (Yes/No) |
| **Metabolism** | CYP1A2/2C19/2C9/2D6/3A4 inhibition |
| **Solubility** | ESOL, Ali, SiRMS (log mol/L + class) |
| **Medicinal chemistry** | PAINS alerts, Brenk alerts, QED, Bioavailability score |

---

## Outputs

| File | Description |
|------|-------------|
| `admet_integrated_full.csv` | All compounds with all properties merged |
| `admet_drug_like_hits.csv` | Filtered drug-like candidates |
| `admet_heatmap.png` | Pass/fail heatmap across all safety filters |
| `bioavailability_radar.png` | SwissADME-style radar for top 5 hits |
| `toxicity_summary.png` | % of compounds passing each safety filter |

---

## Drug-Likeness Filter Criteria Applied
```
✓ Lipinski Ro5 Pass     (MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10)
✓ GI Absorption = High  (Veber: TPSA ≤ 140, RotBonds ≤ 10)
✓ AMES mutagenicity = No
✓ hERG I inhibitor = No (cardiac safety)
✓ Hepatotoxicity = No
✓ PAINS alerts = 0
```

---

## Requirements
```bash
pip install pandas numpy matplotlib seaborn rdkit
```

---

## References
- Pires DE et al. pkCSM. *J Med Chem.* 2015;58(9):4066. [doi:10.1021/acs.jmedchem.5b00104](https://doi.org/10.1021/acs.jmedchem.5b00104)
- Daina A et al. SwissADME. *Sci Rep.* 2017;7:42717. [doi:10.1038/srep42717](https://doi.org/10.1038/srep42717)
