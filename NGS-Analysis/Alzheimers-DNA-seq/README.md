# 🧬 DNA-seq / WES: Alzheimer's Disease Variant Calling

**Dataset:** ADNI (phs000572) — WES data from AD patients and controls  
**Author:** Sana Begum | Research Analyst, Era's Lucknow Medical College

## Pipeline Overview (GATK4 Best Practices)
| Step | Tool | Output |
|------|------|--------|
| Quality control | FastQC + MultiQC | QC report |
| Adapter trimming | fastp | Cleaned FASTQ |
| Alignment | BWA-MEM2 → hg38 | Sorted BAM |
| Deduplication | GATK MarkDuplicates | Dedup BAM |
| Base recalibration | GATK BQSR | Final BAM |
| Variant calling | GATK HaplotypeCaller | gVCF |
| Genotyping | GATK GenotypeGVCFs | Raw VCF |
| Filtering | GATK hard filters | Filtered VCF |
| Annotation | ANNOVAR (RefSeq, ClinVar, ExAC, CADD) | Annotated CSV |
| AD gene analysis | R + bcftools | AD variant summary |

## Quick Start
```bash
# Single sample
bash alzheimers_dnaseq_pipeline.sh SAMPLE001 R1.fastq.gz R2.fastq.gz

# Annotation analysis (R)
Rscript variant_annotation_analysis.R \
    --annovar_csv results_dnaseq/SAMPLE001/annotated/SAMPLE001_anno.hg38_multianno.csv \
    --output results_dnaseq/annotation/
```

## Key AD Variants Analysed
| rsID | Gene | Variant | Clinical Significance |
|------|------|---------|----------------------|
| rs429358 | APOE | ε4 allele (C) | Major late-onset AD risk |
| rs7412 | APOE | ε2 allele (T) | Protective |
| rs75932628 | TREM2 | R47H | Rare AD risk variant |
| rs63750847 | APP | V717I | Early-onset AD |

## Data Access
- ADNI: https://adni.loni.usc.edu/ (requires registration)
- NIAGADS NG00067: open-access AD WES
