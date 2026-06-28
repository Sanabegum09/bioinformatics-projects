# 🧬 RNA-seq: Alzheimer's Disease Progression

**Dataset:** GSE110226 — Human hippocampal RNA-seq: Control vs Alzheimer's Disease  
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College

## Pipeline Overview
| Step | Tool | Output |
|------|------|--------|
| Data download | GEOquery (R) | Raw count matrix |
| Pre-filtering | DESeq2 | Filtered count matrix |
| Normalisation | DESeq2 VST | Normalised expression |
| Differential expression | DESeq2 + apeglm LFC shrinkage | DE gene table |
| PCA & QC | ggplot2, pheatmap | QC plots |
| Volcano plot | EnhancedVolcano | Publication figure |
| Heatmap | pheatmap | Top 50 DE genes |
| AD hallmark genes | Custom analysis | APOE, APP, PSEN1, MAPT, TREM2 etc. |
| Pathway enrichment | clusterProfiler (KEGG + GO BP) | Enrichment plots |
| GSEA export | clusterProfiler | .rnk file for GSEA desktop |

## Quick Start
```r
# Install packages
BiocManager::install(c("GEOquery","DESeq2","EnhancedVolcano","clusterProfiler",
                       "org.Hs.eg.db","biomaRt","ggrepel","pheatmap"))
install.packages(c("ggplot2","dplyr","tibble","viridis","ggpubr","optparse"))

# Run (auto-downloads GSE110226 from GEO)
Rscript alzheimers_rnaseq_pipeline.R --output results_rnaseq/

# With local data
Rscript alzheimers_rnaseq_pipeline.R \
    --local_counts count_matrix.csv \
    --local_meta   metadata.csv \
    --output results_rnaseq/
```

## Key AD Genes Tracked
`APOE · APP · PSEN1 · PSEN2 · MAPT · TREM2 · BIN1 · CLU · CR1 · BACE1 · SNCA · GBA · TNF · IL1B · GFAP · IBA1`

## Multi-Stage Reference
For Braak-staged AD progression data, see:
- **ROSMAP** (Synapse: syn3219045) — temporal cortex, dorsolateral PFC
- **Mayo RNAseq** — TCX + CBE across AD stages
