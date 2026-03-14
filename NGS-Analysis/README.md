# NGS Analysis Projects

This directory contains projects focused on Next-Generation Sequencing (NGS) data processing, quality control, alignment, variant calling, and expression quantification.

## Project Categories

### 1. Quality Control & Data Preprocessing
- Quality Assessment: FastQC analysis, read quality metrics, GC content analysis, contamination detection
- Read Trimming & Filtering: Adapter removal, quality trimming, low-quality read removal
- Read Normalization: Digital normalization, error correction, PCR duplicate removal

### 2. Reference Genome & Indexing
- Reference genome download and preparation
- Genome indexing (BWA, Bowtie2, STAR)
- Annotation file preparation (GTF/GFF)

### 3. Read Alignment & Mapping
- DNA Sequencing Alignment: BWA, Bowtie2
- RNA Sequencing Alignment: HISAT2, STAR, TopHat
- Alignment Quality Control: Mapping statistics, coverage analysis

### 4. SAM/BAM File Processing
- SAM to BAM conversion, sorting and indexing
- Duplicate marking, base quality recalibration
- Indel realignment

### 5. Variant Calling
- SNP & Indel Detection: GATK, SAMtools, FreeBayes
- Variant Quality Metrics: QUAL score, DP, MQ filtering
- Variant Annotation: SnpEff, VEP, database annotation

### 6. RNA Expression Analysis
- Transcript Quantification: RSEM, Kallisto, Salmon, HTSeq
- Differential Expression Analysis: DESeq2, edgeR, Limma
- Quality Control: PCA plots, clustering heatmaps, batch effect detection

### 7. Genome-Wide Association Studies (GWAS)
- SNP selection and filtering, population stratification
- Association testing, Manhattan and QQ plots

### 8. Variant Effect Analysis
- Copy number variation detection, structural variant identification
- Translocation detection, loss of heterozygosity analysis

## Tools & Software
FastQC, MultiQC, Cutadapt, BWA, Bowtie2, HISAT2, STAR, Samtools, Bcftools, GATK, FreeBayes, SnpEff, VEP, RSEM, Kallisto, Salmon, HTSeq, DESeq2, edgeR, Limma, IGV

## Workflow Pipeline
Data Acquisition → Quality Control → Read Processing → Reference Preparation → Read Alignment → Post-Alignment Processing → Variant Calling → Analysis & Interpretation

## Key Parameters & Metrics
- Depth/Coverage: 20-100x typical
- Mapping Rate: >85% expected
- Q Score: Phred score typically Q30+
- QUAL Score: >20 recommended for variants

## Getting Started
Install required tools, prepare reference genomes, run quality control, trim and filter reads, align to reference, call variants or quantify expression

## References
- Li, H. (2013). Aligning sequence reads with BWA-MEM
- Dobin, A., et al. (2013). STAR: ultrafast universal RNA-seq aligner
- Love, M. I., et al. (2014). DESeq2 for RNA-seq analysis
