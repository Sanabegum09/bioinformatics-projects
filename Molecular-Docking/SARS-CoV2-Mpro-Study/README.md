# 🦠 SARS-CoV-2 Main Protease (Mpro) Inhibitor Study
## In-Silico Identification of Phytochemical Inhibitors

**Author:** Sana Begum  
**Affiliation:** Research Analyst, Era's Lucknow Medical College & Hospital, Lucknow  
**Email:** begumsana686@gmail.com  
**GitHub:** [@Sanabegum09](https://github.com/Sanabegum09)  

---

## 📋 Overview

This repository documents the complete computational workflow for the study:

> **"Identification of the most potent bioactive natural compound as main protease inhibitor of SARS-CoV-2: Molecular docking, molecular dynamics simulations and MM-PBSA studies"**  
> Begum S, Singh V, Kumari P, Som A. (2024) — [ResearchGate](https://www.researchgate.net/publication/384256260)

The SARS-CoV-2 main protease (Mpro / 3CLpro) is an attractive antiviral drug target because it is essential for viral polyprotein processing and has no close human homolog. This study performs a systematic in-silico screen of natural phytochemicals to identify potent Mpro inhibitors.

---

## 🔬 Workflow

```
1. Target Preparation
   └── PDB: 6LU7 (Mpro with N3 inhibitor)
       └── Remove ligand & water → add polar H → assign Gasteiger charges

2. Compound Library
   └── Phytochemical database (~500 compounds)
       └── Literature-curated natural products with antiviral activity
       └── Converted to PDBQT with Open Babel

3. Virtual Screening (AutoDock Vina)
   └── Grid box: active site (His41, Cys145)
   └── Exhaustiveness = 32 (for production)
   └── Top 20 hits shortlisted

4. ADMET Profiling
   └── SwissADME, pkCSM, admetSAR
   └── Lipinski Ro5 compliance
   └── BBB, GI absorption, P-gp substrate

5. Molecular Dynamics Simulation (GROMACS)
   └── Force field: CHARMM36 (protein) + CGenFF (ligands)
   └── Simulation time: 100 ns
   └── Periodic boundary conditions (TIP3P water)
   └── NPT ensemble, temperature: 300 K, pressure: 1 bar

6. MM-PBSA Free Energy Calculation
   └── Binding free energy: ΔG_bind = ΔG_complex − (ΔG_protein + ΔG_ligand)
   └── gmx_MMPBSA tool
   └── Per-residue decomposition analysis

7. Results & Comparison
   └── Ranked by binding affinity, ΔG_bind (MM-PBSA), RMSD stability
```

---

## 🖥️ System Setup & Requirements

### Software
| Tool | Version | Use |
|------|---------|-----|
| GROMACS | 2022+ | MD simulation |
| AutoDock Vina | 1.2.3 | Molecular docking |
| Open Babel | 3.1+ | Format conversion |
| PyMOL | 2.5+ | Visualisation |
| UCSF Chimera | 1.16+ | Structure prep |
| Python | 3.8+ | Analysis scripts |
| gmx_MMPBSA | 1.6+ | MM-PBSA calculations |

### Python dependencies
```bash
pip install numpy pandas matplotlib seaborn biopython rdkit
pip install gmx_MMPBSA
```

---

## 📁 Repository Structure

```
sars-cov2-mpro-study/
├── README.md                    # This file
├── data/
│   ├── 6LU7_receptor.pdb        # Crystal structure (from PDB)
│   ├── 6LU7_prepared.pdbqt      # Prepared receptor
│   └── phytochemical_library/   # Screened compound library
├── docking/
│   ├── vina_config.txt          # Grid box configuration
│   ├── screening_results.csv    # All docking scores
│   └── top20_hits/              # PDBQT files for top hits
├── md_simulation/
│   ├── topol.top                # GROMACS topology
│   ├── em.mdp                   # Energy minimisation parameters
│   ├── nvt.mdp                  # NVT equilibration parameters
│   ├── npt.mdp                  # NPT equilibration parameters
│   └── prod.mdp                 # Production MD parameters
├── mmpbsa/
│   ├── mmpbsa.in                # MM-PBSA input file
│   └── FINAL_RESULTS_MMPBSA.dat # Final binding energies
├── scripts/
│   ├── prepare_receptor.sh      # Receptor preparation script
│   ├── run_screening.py         # Docking pipeline (uses ../molecular-docking-pipeline/)
│   ├── analyse_trajectory.py    # RMSD, RMSF, Rg analysis
│   └── plot_results.py          # Publication-quality figures
└── results/
    ├── top_hits_admet.csv       # ADMET of top hits
    ├── md_analysis/             # Trajectory analysis outputs
    └── figures/                 # Publication figures
```

---

## ⚙️ Running the Analysis

### Step 1: Receptor Preparation
```bash
# Download from RCSB PDB
wget https://files.rcsb.org/download/6LU7.pdb

# Remove HETATM (ligand/water) and save receptor only
grep "^ATOM" 6LU7.pdb > 6LU7_receptor.pdb

# Add H and convert to PDBQT using AutoDockTools (prepare_receptor4.py)
python prepare_receptor4.py -r 6LU7_receptor.pdb -o data/6LU7_prepared.pdbqt \
    -A hydrogens -U nphs_lps_waters_nonstdres
```

### Step 2: Virtual Screening
```bash
# Using the molecular docking pipeline from this repository
python ../molecular-docking-pipeline/autodock_vina_pipeline.py \
    --receptor data/6LU7_prepared.pdbqt \
    --ligands_dir data/phytochemical_library/ \
    --config docking/vina_config.txt \
    --output docking/results/ \
    --top_n 20
```

### Step 3: ADMET Filtering
```bash
python ../admet-screening/admet_screening.py \
    --input docking/top20_hits.sdf \
    --output results/admet/ \
    --plot
```

### Step 4: MD Simulation (GROMACS)
```bash
# Full workflow — see md_simulation/ for .mdp files

# Energy minimisation
gmx grompp -f md_simulation/em.mdp -c complex.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# NVT equilibration (100 ps)
gmx grompp -f md_simulation/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
gmx mdrun -v -deffnm nvt

# NPT equilibration (100 ps)
gmx grompp -f md_simulation/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -n index.ndx -o npt.tpr
gmx mdrun -v -deffnm npt

# Production MD (100 ns)
gmx grompp -f md_simulation/prod.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o prod.tpr
gmx mdrun -v -deffnm prod -ntmpi 1 -ntomp 8
```

### Step 5: MM-PBSA
```bash
gmx_MMPBSA -O -i mmpbsa/mmpbsa.in \
    -cs prod.tpr -ct prod.xtc \
    -ci index.ndx -cg 1 13 \
    -cp topol.top -o FINAL_RESULTS_MMPBSA.dat \
    -eo FINAL_RESULTS_MMPBSA.csv
```

---

## 📊 Key Results

| Compound | Docking Score (kcal/mol) | ΔG_bind MM-PBSA (kcal/mol) | RMSD Stability |
|----------|--------------------------|---------------------------|----------------|
| [Best hit] | —   | — | Stable (< 2.0 Å) |
| N3 (control) | -7.2 | — | Stable |

> Full results available in `results/` after running the pipeline.

---

## 📚 Key References

1. Jin Z et al. *Structure of Mpro from SARS-CoV-2 and discovery of its inhibitors.* Nature. 2020;582:289–293. [PDB: 6LU7]
2. Trott O, Olson AJ. *AutoDock Vina: improving the speed and accuracy of docking.* J Comput Chem. 2010;31(2):455–461.
3. Begum S, Singh V, Kumari P, Som A. *Identification of the most potent bioactive natural compound as Mpro inhibitor of SARS-CoV-2.* ResearchGate (2024). https://www.researchgate.net/publication/384256260

---

## 📫 Contact

**Sana Begum**  
Research Analyst | Era's Lucknow Medical College & Hospital  
📧 begumsana686@gmail.com  
🔗 [LinkedIn](https://www.linkedin.com/in/sana-begum-964a57215) | [ResearchGate](https://www.researchgate.net/profile/Sana-Begum)

---
*"Computational drug discovery to combat infectious diseases and cancer."*
