# README for Phylogenetic Analysis Projects

## Overview
This repository contains various projects focusing on phylogenetic analysis, a method used to infer the evolutionary relationships among biological entities such as species, genes, or populations.

## Topics Covered

### Multiple Sequence Alignment
- **Definition**: A method of aligning multiple sequences (DNA, RNA, or proteins) to identify regions of similarity.
- **Common Tools**:
  - ClustalW
  - MUSCLE
  - T-Coffee

### Tree Construction Methods
- **Phylogenetic Trees**: Representations of the evolutionary relationships.
- **Methods**:
  - **Maximum Likelihood**: Estimates the tree with the highest likelihood based on the given data.
  - **Bayesian Inference**: Uses Bayesian methods to estimate the probability of different trees.
  - **Neighbor-Joining**: Calculates pairwise distance and builds a tree based on that.

### Evolutionary Distance Calculations
- **Purpose**: Measures the genetic divergence between species.
- **Methods**:
  - Jukes-Cantor Model
  - Kimura 2-parameter Model

### Analysis Tools
#### RAxML
- **Workflow**: Upload aligned sequences, choose the substitution model, and run the analysis.
- **Key Parameters**: 
  - Model of substitution
  - Bootstrap iterations

#### PhyML
- **Workflow**: Similar to RAxML, it evaluates the likelihood of trees based on different models.
- **Key Parameters**:
  - Substitution model
  - Number of bootstrap replicates

#### MrBayes
- **Workflow**: Sets up Markov chain Monte Carlo analysis to estimate posterior distributions.
- **Key Parameters**:
  - Number of generations
  - Prior distributions

#### MEGA
- **Workflow**: User-friendly interface for analyzing DNA/protein sequences and constructing phylogenetic trees.
- **Key Parameters**:
  - Distance matrix methods
  - Tree-building algorithms

#### FigTree
- **Purpose**: Visualization of phylogenetic trees.
- **Workflow**: Import tree file and customize display options.

#### iTOL
- **Purpose**: Interactive Tree Of Life for visualizing complex trees.
- **Workflow**: Upload tree and annotate with additional data.

## Conclusion
Phylogenetic analysis involves complex processes and various tools to derive meaningful insights into evolutionary relationships. This repository serves as a guide to these processes and tools, offering complete workflows and suggested parameters for users to achieve successful analysis.