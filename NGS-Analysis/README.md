# 🧬 NGS Analysis: Alzheimer's Disease Multi-Omics

**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College & Hospital  
**Email:** begumsana686@gmail.com | **GitHub:** [@Sanabegum09](https://github.com/Sanabegum09)

A comprehensive multi-type NGS analysis repository focussed on Alzheimer's disease,
demonstrating expertise across RNA-seq, miRNA-seq, and DNA-seq/WES analysis.

## Repository Structure
```
NGS-Analysis/
├── README.md                         # This file
├── Alzheimers-RNA-seq/
│   ├── alzheimers_rnaseq_pipeline.R  # Full DESeq2 pipeline (GSE110226)
│   └── README.md
├── Alzheimers-miRNA-seq/
│   ├── alzheimers_mirna_pipeline.R   # miRNA-seq DESeq2 pipeline (GSE46579)
│   └── README.md
└── Alzheimers-DNA-seq/
    ├── alzheimers_dnaseq_pipeline.sh # GATK4 WES variant calling (ADNI)
    ├── variant_annotation_analysis.R # ANNOVAR results + AD gene analysis
    └── README.md
```

## Datasets Used
| Analysis | Dataset | Source | Access |
|----------|---------|--------|--------|
| RNA-seq | GSE110226 | NCBI GEO | Public |
| miRNA-seq | GSE46579 | NCBI GEO | Public |
| DNA-seq/WES | ADNI (phs000572) | dbGaP | Requires registration |

## Key Alzheimer's Findings Targeted
- **Transcriptomic:** Dysregulation of APOE, APP, PSEN1, TREM2, MAPT, BACE1
- **miRNA:** Downregulation of miR-132, miR-107, miR-29 family; upregulation of miR-146a
- **Genomic:** APOE ε4 (rs429358), TREM2 R47H (rs75932628), PSEN1/2 pathogenic variants

## Tools & Software
### NGS
`FastQC · MultiQC · fastp · BWA-MEM2 · Samtools · GATK4 · bcftools · ANNOVAR`

### R/Bioconductor
`GEOquery · DESeq2 · EnhancedVolcano · clusterProfiler · org.Hs.eg.db · pheatmap · ggplot2`

## Common Databases Referenced
- NCBI GEO: https://www.ncbi.nlm.nih.gov/geo/
- ADNI: https://adni.loni.usc.edu/
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- miRTarBase: https://mirtarbase.cuhk.edu.cn/
- KEGG AD pathway: hsa05010
