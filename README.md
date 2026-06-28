<div align="center">

<img src="https://capsule-render.vercel.app/api?type=waving&color=gradient&customColorList=12&height=220&section=header&text=Bioinformatics%20Research%20Portfolio&fontSize=38&fontColor=FFFFFF&fontAlignY=40&desc=Sana%20Begum%20%7C%20Computational%20Drug%20Discovery%20%7C%20Structural%20Bioinformatics%20%7C%20NGS%20Analysis&descAlignY=62&descSize=15&descColor=D6EAF8&animation=fadeIn" width="100%"/>

</div>

<div align="center">

[![LinkedIn](https://img.shields.io/badge/LinkedIn-0A66C2?style=for-the-badge&logo=linkedin&logoColor=white)](https://www.linkedin.com/in/sana-begum-964a57215)
[![ResearchGate](https://img.shields.io/badge/ResearchGate-00CCBB?style=for-the-badge&logo=researchgate&logoColor=white)](https://www.researchgate.net/profile/Sana-Begum)
[![Email](https://img.shields.io/badge/Email-EA4335?style=for-the-badge&logo=gmail&logoColor=white)](mailto:begumsana686@gmail.com)
[![GitHub followers](https://img.shields.io/github/followers/Sanabegum09?style=for-the-badge&color=1A5276)](https://github.com/Sanabegum09)

</div>

---

## 👩‍🔬 About This Repository

This is the research portfolio of **Sana Begum**, Bioinformatics Research Analyst at [Era's Lucknow Medical College & Hospital](https://www.eralucknow.com/), Lucknow. It contains complete, reproducible bioinformatics pipelines built from **3+ years of active research** spanning computational drug discovery and multi-omics analysis.

| | |
|---|---|
| 🏛 **Position** | Research Analyst, Era's Lucknow Medical College & Hospital |
| 🎓 **Education** | M.Sc. Bioinformatics, University of Allahabad *(CGPA: 9.42 / 10)* |
| 🏆 **Qualification** | ICAR-JRF 2022 · CRET-22 |
| 📄 **Publications** | 2 first-author · 2 under review |
| 🔬 **Expertise** | Structural Bioinformatics · Molecular Docking · MD Simulation · NGS |
| 💊 **Disease areas** | Cancer · Neurodegeneration · Infectious disease · Metabolic · Autoimmune |

---

## 📄 Publications

> ### [1] Phytochemicals as DNMTs Inhibitors in Head & Neck Cancer
> **Begum S**, Fatima G, Singh P, Naeem A, Fatima N, Siddiqui Z.
> *"Phytochemicals as DNMTs inhibitors for targeted therapy in head and neck cancer: An in-silico study."*
> **Results in Chemistry** 18 (2025) 102783.
>
> [![DOI](https://img.shields.io/badge/DOI-10.1016/j.rechem.2025.102783-1A5276?style=flat-square)](https://doi.org/10.1016/j.rechem.2025.102783)
> ![Tools](https://img.shields.io/badge/Tools-AutoDock_Vina_·_GROMACS_·_SwissADME_·_pkCSM_·_MM--PBSA-117A65?style=flat-square)

> ### [2] SARS-CoV-2 Main Protease Inhibitor Identification
> **Begum S**, Singh V, Kumari P, Som A.
> *"Identification of the most potent bioactive natural compound as main protease inhibitor of SARS-CoV-2: Molecular docking, molecular dynamics simulations and MM-PBSA studies."* (2024)
>
> [![ResearchGate](https://img.shields.io/badge/ResearchGate-View_Publication-00CCBB?style=flat-square&logo=researchgate)](https://www.researchgate.net/publication/384256260)
> ![Tools](https://img.shields.io/badge/Tools-AutoDock_Vina_·_GROMACS_·_SwissADME_·_pkCSM_·_PyMOL-117A65?style=flat-square)

> ### [3] Isatin Derivatives as ALK Inhibitors in NSCLC *(Under Review)*
> Design, synthesis and biological evaluation of Isatin derivatives as probable ALK inhibitors against non-small cell lung cancer.

> ### [4] Violacein as BMI1 Inhibitor in Bladder Cancer *(In Preparation)*
> Computational evaluation of Violacein and its intermediates as BMI1 inhibitors in bladder cancer.

---

## 🗂️ Repository Contents

> **Note:** Only folders listed below currently contain scripts. Other areas are in active development.

### 🔬 [Molecular-Docking/](./Molecular-Docking/)
*Automated virtual screening pipeline — core of all published drug discovery work*

| File | Description |
|------|-------------|
| [`autodock_vina_pipeline.py`](./Molecular-Docking/autodock_vina_pipeline.py) | Batch AutoDock Vina virtual screening with hit ranking & CSV export |
| [`SARS-CoV2-Mpro-Study/`](./Molecular-Docking/SARS-CoV2-Mpro-Study/) | Complete workflow from published SARS-CoV-2 Mpro inhibitor study (2024) |

```python
# Screen 500 compounds against a target in one command
python autodock_vina_pipeline.py \
    --receptor target.pdbqt \
    --ligands_dir library_pdbqt/ \
    --center_x -26.5 --center_y 15.2 --center_z -16.8 \
    --output screening_results/ --top_n 20
```

---

### 💊 [ADMET-Prediction/](./ADMET-Prediction/)
*Replicates the exact ADMET workflow used in published research (Results in Chemistry, 2025)*

**Servers used:** [pkCSM](https://biosig.lab.uq.edu.au/pkcsm/) · [SwissADME](https://www.swissadme.ch/)

| File | Description |
|------|-------------|
| [`admet_screening.py`](./ADMET-Prediction/admet_screening.py) | Parses & integrates pkCSM + SwissADME CSV results; generates heatmap, radar, toxicity plots |
| [`requirements.txt`](./ADMET-Prediction/requirements.txt) | Python dependencies |

```
Workflow:
  top docking hits → prepare SMILES → submit to pkCSM + SwissADME
  → download CSVs → python admet_screening.py → integrated report + plots
```

**Properties predicted:** Caco-2 · HIA · P-gp · VDss · BBB · CYP1A2/2C9/2C19/2D6/3A4 · AMES · hERG · Hepatotoxicity · Lipinski Ro5 · TPSA · QED · PAINS · Bioavailability radar

---

### 🧬 [NGS-Analysis/](./NGS-Analysis/)
*Multi-type NGS analysis using Alzheimer's disease public datasets*

| Sub-folder | Script | Dataset | Analysis |
|-----------|--------|---------|----------|
| [`Alzheimers-RNA-seq/`](./NGS-Analysis/Alzheimers-RNA-seq/) | `alzheimers_rnaseq_pipeline.R` | GSE110226 (GEO) | DESeq2, KEGG/GO enrichment, AD hallmark genes |
| [`Alzheimers-miRNA-seq/`](./NGS-Analysis/Alzheimers-miRNA-seq/) | `alzheimers_mirna_pipeline.R` | GSE46579 (GEO) | DESeq2, miR-29/132/146a spotlight, target export |
| [`Alzheimers-DNA-seq/`](./NGS-Analysis/Alzheimers-DNA-seq/) | `alzheimers_dnaseq_pipeline.sh` + `.R` | ADNI WES | GATK4, APOE ε4 genotyping, ANNOVAR annotation |

```r
# RNA-seq: auto-downloads GSE110226 from GEO and runs full analysis
Rscript NGS-Analysis/Alzheimers-RNA-seq/alzheimers_rnaseq_pipeline.R \
    --output results_rnaseq/
```

---

## 🛠️ Tools & Software

### Computational Drug Discovery
![AutoDock Vina](https://img.shields.io/badge/AutoDock_Vina-1A5276?style=flat-square)
![GROMACS](https://img.shields.io/badge/GROMACS-1A5276?style=flat-square)
![gmx_MMPBSA](https://img.shields.io/badge/gmx__MMPBSA-1A5276?style=flat-square)
![MODELLER](https://img.shields.io/badge/MODELLER-1A5276?style=flat-square)
![AlphaFold2](https://img.shields.io/badge/AlphaFold2-1A5276?style=flat-square)
![SWISS-MODEL](https://img.shields.io/badge/SWISS--MODEL-1A5276?style=flat-square)
![SwissADME](https://img.shields.io/badge/SwissADME-1A5276?style=flat-square)
![pkCSM](https://img.shields.io/badge/pkCSM-1A5276?style=flat-square)

### Visualisation
![PyMOL](https://img.shields.io/badge/PyMOL-117A65?style=flat-square)
![UCSF Chimera](https://img.shields.io/badge/UCSF_Chimera-117A65?style=flat-square)
![Discovery Studio](https://img.shields.io/badge/Discovery_Studio-117A65?style=flat-square)
![LigPlot+](https://img.shields.io/badge/LigPlot+-117A65?style=flat-square)
![Grace](https://img.shields.io/badge/Grace-117A65?style=flat-square)

### NGS & Genomics
![FastQC](https://img.shields.io/badge/FastQC-8E44AD?style=flat-square)
![HISAT2](https://img.shields.io/badge/HISAT2-8E44AD?style=flat-square)
![STAR](https://img.shields.io/badge/STAR-8E44AD?style=flat-square)
![BWA--MEM2](https://img.shields.io/badge/BWA--MEM2-8E44AD?style=flat-square)
![GATK4](https://img.shields.io/badge/GATK4-8E44AD?style=flat-square)
![Samtools](https://img.shields.io/badge/Samtools-8E44AD?style=flat-square)
![DESeq2](https://img.shields.io/badge/DESeq2-8E44AD?style=flat-square)
![ANNOVAR](https://img.shields.io/badge/ANNOVAR-8E44AD?style=flat-square)

### Programming
![Python](https://img.shields.io/badge/Python_3-3776AB?style=flat-square&logo=python&logoColor=white)
![R](https://img.shields.io/badge/R_4+-276DC3?style=flat-square&logo=r&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-4EAA25?style=flat-square&logo=gnubash&logoColor=white)
![Linux](https://img.shields.io/badge/Linux_HPC-FCC624?style=flat-square&logo=linux&logoColor=black)

---

## 🔬 Drug Discovery Workflow

```
                 ┌─────────────────────────────────────────────┐
                 │       COMPUTATIONAL DRUG DISCOVERY          │
                 │            Sana Begum — 10+ targets         │
                 └─────────────────────────────────────────────┘
                                      │
         ┌──────────────┬─────────────┼───────────────┬──────────────┐
         ▼              ▼             ▼               ▼              ▼
   Target Prep    Compound Lib   Docking (Vina)  ADMET Filter   Hit Selection
   PDB / MODELLER  PubChem/ZINC  Virtual Screen  pkCSM+SwissADME  Top 3–5 hits
   AlphaFold2      Literature    Rank by ΔG      Lipinski/hERG   for MD
         │              │             │               │              │
         └──────────────┴─────────────┴───────────────┴──────────────┘
                                      │
                          ┌───────────▼───────────┐
                          │  MD Simulation (100ns) │
                          │  GROMACS · CHARMM36    │
                          │  NVT → NPT → Production│
                          └───────────┬───────────┘
                                      │
                          ┌───────────▼───────────┐
                          │  MM-PBSA / gmx_MMPBSA  │
                          │  ΔGbind · RMSD · RMSF  │
                          │  Rg · H-bonds · SASA   │
                          └───────────┬───────────┘
                                      │
                          ┌───────────▼───────────┐
                          │  Interaction Analysis  │
                          │  PyMOL · LigPlot+      │
                          │  Discovery Studio      │
                          └───────────┬───────────┘
                                      │
                               📄 Publication
```

---

## 🦠 Disease Areas Covered

| Domain | Targets / Disease |
|--------|------------------|
| **Oncology** | Head & Neck Cancer (DNMT) · NSCLC (ALK) · Bladder Cancer (BMI1) · Prostate Cancer · Melanoma |
| **Infectious Disease** | SARS-CoV-2 (Mpro) · Dengue |
| **Neurodegeneration** | Alzheimer's Disease · Parkinson's Disease · Autism Spectrum Disorder |
| **Metabolic** | Diabetes · Thyroid Disorders |
| **Autoimmune** | Rheumatoid Arthritis · Psoriasis |
| **Epigenetics** | DNMT1/3A/3B inhibitors · Histone modification |

---

## ⚡ Quick Start

```bash
# Clone the repository
git clone https://github.com/Sanabegum09/bioinformatics-projects.git
cd bioinformatics-projects

# Install Python dependencies (docking + ADMET)
pip install -r Molecular-Docking/requirements.txt
pip install -r ADMET-Prediction/requirements.txt

# Install R dependencies (NGS)
Rscript -e "BiocManager::install(c('GEOquery','DESeq2','EnhancedVolcano','clusterProfiler','org.Hs.eg.db'))"
Rscript -e "install.packages(c('ggplot2','pheatmap','ggrepel','dplyr','viridis','optparse'))"
```

---

## 📊 GitHub Stats

<div align="center">

<img height="160" src="https://github-readme-stats.vercel.app/api?username=Sanabegum09&show_icons=true&title_color=1A5276&icon_color=1A5276&border_radius=10&count_private=false&hide_border=false"/>
<img height="160" src="https://github-readme-stats.vercel.app/api/top-langs/?username=Sanabegum09&layout=compact&title_color=1A5276&border_radius=10&langs_count=6&hide_border=false"/>

</div>

---

## 📫 Contact & Collaboration

I welcome collaborations in computational drug discovery, multi-omics analysis, and bioinformatics pipeline development.

<div align="center">

| | |
|---|---|
| 📧 **Email** | begumsana686@gmail.com |
| 💼 **LinkedIn** | [linkedin.com/in/sana-begum-964a57215](https://www.linkedin.com/in/sana-begum-964a57215) |
| 🔬 **ResearchGate** | [researchgate.net/profile/Sana-Begum](https://www.researchgate.net/profile/Sana-Begum) |
| 🏛 **Institution** | Era's Lucknow Medical College & Hospital, Lucknow |

</div>

---

<div align="center">

*"Computational approaches to accelerate drug discovery across disease frontiers."*

<img src="https://capsule-render.vercel.app/api?type=waving&color=gradient&customColorList=12&height=100&section=footer" width="100%"/>

</div>
