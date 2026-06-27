# 🧬 NGS Analysis — miRNA-seq Differential Expression

Complete **miRNA-seq** differential expression analysis pipeline using **DESeq2**. Currently applied to bladder cancer miRNA profiling (ongoing research, Era's Lucknow Medical College).

## Pipeline Steps
1. Raw count matrix loading & sample QC
2. Low-count pre-filtering (≥10 reads in ≥2 samples)
3. DESeq2 normalisation + LFC shrinkage (apeglm)
4. PCA & sample distance heatmap
5. Volcano plot (EnhancedVolcano)
6. MA plot
7. Top DE miRNA heatmap (pheatmap)
8. Export of ranked lists for GSEA / target prediction

## Requirements
```r
install.packages(c("optparse","ggplot2","ggrepel","pheatmap","RColorBrewer",
                   "dplyr","tidyr","tibble"))
BiocManager::install(c("DESeq2","clusterProfiler","org.Hs.eg.db","EnhancedVolcano"))
```

## Quick Start
```bash
Rscript mirna_deseq2_analysis.R \
    --counts count_matrix.csv \
    --metadata metadata.csv \
    --output results/ \
    --condition Condition \
    --reference Normal \
    --padj 0.05 \
    --lfc 1.0
```

## Input Format
**count_matrix.csv** — rows = miRNA names, columns = sample IDs, values = raw read counts  
**metadata.csv** — rows = sample IDs, must contain a column matching `--condition`

## Output
```
results/
├── QC/
│   ├── pca_plot.pdf/.png
│   └── sample_distance_heatmap.pdf
├── DE/
│   ├── deseq2_all_results.csv
│   ├── deseq2_significant_miRNAs.csv
│   ├── upregulated_miRNAs.txt
│   ├── downregulated_miRNAs.txt
│   └── ranked_mirnas_for_gsea.csv
├── Plots/
│   ├── volcano_plot.pdf/.png
│   ├── ma_plot.pdf/.png
│   └── top_mirna_heatmap.pdf
└── session_info.txt
```

---
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College | [GitHub](https://github.com/Sanabegum09)
