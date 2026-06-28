# 🧬 miRNA-seq: Alzheimer's Disease

**Dataset:** GSE46579 — miRNA expression profiling in AD brain (temporal cortex)  
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College

## Pipeline Overview
| Step | Tool | Output |
|------|------|--------|
| Data download | GEOquery | Count matrix |
| Pre-filtering | DESeq2 | Filtered miRNAs |
| Differential expression | DESeq2 + apeglm | DE miRNA table |
| PCA | ggplot2 + stat_ellipse | QC plot |
| Volcano plot | EnhancedVolcano | Publication figure |
| Heatmap | pheatmap | Top 40 DE miRNAs |
| AD miRNA spotlight | Custom | miR-9, miR-29, miR-107, miR-132 etc. |
| Target gene export | — | Lists for miRTarBase/TargetScan |
| GSEA export | — | Ranked .csv |

## Quick Start
```r
Rscript alzheimers_mirna_pipeline.R --output results_mirna/
```

## Key AD miRNAs Tracked
`miR-9 · miR-29a/b/c · miR-107 · miR-132 · miR-146a · miR-155 · miR-181a · miR-21 · miR-34a · miR-125b`

## Target Prediction Resources
Submit exported miRNA lists to:
- [miRTarBase](https://mirtarbase.cuhk.edu.cn/)
- [TargetScan](https://www.targetscan.org/)
- [miRDB](https://www.mirdb.org/)
- [DIANA-miRPath](https://www.microrna.gr/miRPathv3)
