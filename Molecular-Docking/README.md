# Molecular Docking using AutoDock

Molecular docking is a computational method used to predict the interaction between a receptor and a ligand by identifying the preferred binding mode and estimating binding affinity.

---

## Software Required

- AutoDock 4
- MGLTools
- Open Babel
- Discovery Studio Visualizer

---

## Receptor Preparation

- Open receptor structure (PDB format)
- Remove water molecules
- Remove hetero atoms/ligands
- Check for missing atoms
- Add polar hydrogens
- Assign Kollman charges
- Verify total charges
- Save receptor as `receptor.pdbqt`

---

## Ligand Preparation

- Open ligand structure
- Detect root atom
- Define torsion tree
- Set rotatable bonds
- Save ligand as `ligand.pdbqt`

---

## Grid Parameter

- Load receptor
- Select ligand map types
- Define grid box dimensions
- Save Grid Parameter File (`.gpf`)
- Run AutoGrid

---

## Docking Parameters

- Load receptor
- Load ligand
- Configure Lamarckian Genetic Algorithm
- Set population size
- Set number of runs
- Save Docking Parameter File (`.dpf`)

---

## Molecular Docking

- Run AutoDock
- Load the Docking Parameter File (`.dpf`)
- Execute docking

---

## Result Analysis

- Open Docking Log File (`.dlg`)
- Load receptor structure
- Visualize docking conformations
- Analyze interaction energy
- Build hydrogen bonds
- Save docked complex
- Export interaction image
