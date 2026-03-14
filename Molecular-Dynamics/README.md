# Molecular Dynamics Projects

This directory contains projects focused on molecular dynamics simulations, trajectory analysis, and binding free energy calculations.

## Project Categories

### 1. System Preparation
- Protein structure validation and optimization
- Solvation (explicit/implicit solvent models)
- Ionization and salt bridge setup
- Box size and geometry definition

### 2. Energy Minimization
- Steepest descent minimization
- Conjugate gradient optimization
- Convergence criteria and analysis
- Minimization protocols for different systems

### 3. Equilibration Protocols
- NVT Equilibration (canonical ensemble)
  - Temperature coupling and thermostat setup
  - Velocity scaling
- NPT Equilibration (isothermal-isobaric ensemble)
  - Pressure coupling
  - Barostat configuration
- Equilibration time series analysis

### 4. Production MD Simulations
- Long timescale simulations (ns to μs)
- Trajectory generation and storage
- Timestep optimization
- Constraint algorithms (LINCS, SHAKE)
- Integration methods

### 5. Trajectory Analysis
- RMSD (Root Mean Square Deviation) calculations
- RMSF (Root Mean Square Fluctuation) analysis
- Radius of gyration (Rg) calculations
- Hydrogen bond analysis
- Protein flexibility and stability assessment

### 6. Binding Free Energy Calculations
- MM-PBSA (Molecular Mechanics Poisson-Boltzmann Surface Area)
- MM-GBSA (Molecular Mechanics Generalized Born Surface Area)
- Energy decomposition analysis
- Per-residue contribution calculations
- Entropic contributions

### 7. Advanced Analysis
- Secondary structure evolution
- Clustering analysis
- Principal Component Analysis (PCA)
- Free energy landscapes
- Transition pathway analysis

## Tools & Software
- **GROMACS** - Primary MD engine
- **AMBER** - Alternative MD engine with advanced tools
- **NAMD** - Parallel MD simulation
- **VMD** - Visualization and analysis
- **PyMOL** - Structure visualization
- **DSSP** - Secondary structure assignment

## Workflow Pipeline

```
1. Preparation
   ├── Download structure (PDB/AlphaFold)
   ├── Validate structure
   └── Add hydrogens

2. Parameterization
   ├── Assign force field (AMBER, CHARMM, GROMOS)
   ├── Generate topology
   └── Define non-bonded parameters

3. Solvation
   ├── Add solvent molecules
   ├── Add counterions
   └── Define periodic boundary conditions

4. Energy Minimization
   ├── Steepest descent (5000-10000 steps)
   └── Conjugate gradient (until convergence)

5. Equilibration
   ├── NVT equilibration (100-200 ps)
   └── NPT equilibration (200-500 ps)

6. Production Run
   ├── 1-100 ns simulation
   ├── Data collection
   └── Trajectory storage

7. Analysis
   ├── Stability analysis (RMSD, RMSF)
   ├── Binding analysis
   ├── Free energy calculation
   └── Report generation
```

## Key Parameters & Metrics

### Simulation Parameters
- Temperature: 300-310 K (physiological)
- Pressure: 1 atm (1.01325 bar)
- Timestep: 0.001-0.002 ps
- Cutoff for non-bonded interactions: 10-12 Å

### Analysis Metrics
- **RMSD**: Protein conformational stability
- **RMSF**: Per-residue flexibility
- **Rg**: Overall compactness
- **H-bonds**: Structural stability indicators
- **ΔG_binding**: Free energy of binding
- **SASA**: Solvent Accessible Surface Area

## Getting Started

### Prerequisites
```bash
# Install GROMACS
sudo apt-get install gromacs

# Install Python analysis tools
pip install mdtraj numpy pandas matplotlib seaborn scipy biopython
```

### Basic Workflow
1. See individual project folders for complete protocols
2. Refer to setup scripts for automation
3. Check trajectory analysis scripts for post-processing
4. Review visualization templates

## Output Files
- `.gro` - Structure files
- `.xtc` / `.trr` - Trajectory files
- `.edr` - Energy files
- `.log` - Log files
- Analysis results (RMSD, RMSF plots, etc.)

## References
- GROMACS Manual: http://manual.gromacs.org/
- AMBER Documentation: https://ambermd.org/
- MD Simulation Best Practices in Literature

## Troubleshooting
- Blowing up simulations: Check force field parameters, initial geometry
- Poor convergence: Increase minimization steps, review coordinates
- Unrealistic trajectories: Verify topology, check timestep, review constraints