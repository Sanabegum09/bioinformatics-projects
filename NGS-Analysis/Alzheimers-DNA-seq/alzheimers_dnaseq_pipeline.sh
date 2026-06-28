#!/usr/bin/env bash
# ============================================================
#  DNA-seq / WES Variant Calling Pipeline
#  Alzheimer's Disease — ADNI Dataset Reference
# ============================================================
#  Author  : Sana Begum
#  Role    : Research Analyst, Era's Lucknow Medical College
#  Email   : begumsana686@gmail.com
#  GitHub  : https://github.com/Sanabegum09
#
#  Dataset Reference:
#    ADNI (Alzheimer's Disease Neuroimaging Initiative)
#    WES/WGS data — dbGaP accession: phs000572
#    Access: https://adni.loni.usc.edu/
#
#  Open-access alternative:
#    NIAGADS dataset: NG00067 (AD GWAS + WES)
#    1000 Genomes + GATK best-practices (demo)
#
#  Pipeline (GATK4 Best Practices):
#    1. Quality Control  — FastQC + MultiQC
#    2. Trimming         — Trimmomatic / fastp
#    3. Alignment        — BWA-MEM2
#    4. Post-alignment   — Samtools sort/index/flagstat
#    5. Deduplication    — GATK MarkDuplicates
#    6. Base recalibration — GATK BQSR
#    7. Variant calling  — GATK HaplotypeCaller (gVCF mode)
#    8. Joint genotyping — GATK GenomicsDBImport + GenotypeGVCFs
#    9. Variant filtering — GATK VQSR or hard filters
#   10. Annotation       — ANNOVAR / SnpEff
#   11. AD gene analysis — APOE, PSEN1, PSEN2, APP, TREM2 etc.
#
#  Key AD variants tracked:
#    APOE ε4 (rs429358, rs7412)
#    TREM2 R47H (rs75932628)
#    PSEN1, PSEN2 pathogenic mutations
#    APP duplications/mutations
#
#  Tools required (install separately):
#    FastQC, MultiQC, fastp, BWA-MEM2, Samtools,
#    GATK4, ANNOVAR, SnpEff, bcftools, vcftools
# ============================================================

set -euo pipefail

# ── CONFIGURATION ────────────────────────────────────────
SAMPLE_ID="${1:-SAMPLE001}"
FASTQ_R1="${2:-data/raw/${SAMPLE_ID}_R1.fastq.gz}"
FASTQ_R2="${3:-data/raw/${SAMPLE_ID}_R2.fastq.gz}"
REF_GENOME="reference/hg38.fa"
KNOWN_SITES_DBSNP="reference/dbsnp_146.hg38.vcf.gz"
KNOWN_SITES_MILLS="reference/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
KNOWN_SITES_1000G="reference/1000G_phase1.snps.high_confidence.hg38.vcf.gz"

OUTDIR="results_dnaseq/${SAMPLE_ID}"
THREADS=8
GATK="gatk"
JAVA_OPTS="-Xmx16g"

# AD-relevant gene regions (hg38) for targeted analysis
AD_REGIONS=(
  "chr19:44905781-44909393"   # APOE
  "chr14:73136418-73223691"   # PSEN1
  "chr1:227058274-227083858"  # PSEN2
  "chr21:25880550-26171128"   # APP
  "chr6:41158507-41163186"    # TREM2
  "chr10:104893962-105017532" # SORL1
  "chr11:85678154-85904905"   # BIN1
  "chr2:127892549-128086494"  # CLU (APOJ)
  "chr1:207577071-207666922"  # CR1
)

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
mkdir -p "${OUTDIR}"/{qc,trim,align,bqsr,gvcf,vcf,annotated,ad_genes}


# ════════════════════════════════════════════════════════════
# STEP 1: QUALITY CONTROL
# ════════════════════════════════════════════════════════════
log "STEP 1/10: Quality control (FastQC)"

fastqc \
  "${FASTQ_R1}" "${FASTQ_R2}" \
  --outdir "${OUTDIR}/qc" \
  --threads "${THREADS}" \
  --quiet

# MultiQC summary
multiqc "${OUTDIR}/qc" \
  --outdir "${OUTDIR}/qc" \
  --filename "${SAMPLE_ID}_multiqc" \
  --quiet

log "  QC complete → ${OUTDIR}/qc/"


# ════════════════════════════════════════════════════════════
# STEP 2: ADAPTER TRIMMING
# ════════════════════════════════════════════════════════════
log "STEP 2/10: Adapter trimming (fastp)"

TRIM_R1="${OUTDIR}/trim/${SAMPLE_ID}_R1_trimmed.fastq.gz"
TRIM_R2="${OUTDIR}/trim/${SAMPLE_ID}_R2_trimmed.fastq.gz"

fastp \
  --in1  "${FASTQ_R1}" \
  --in2  "${FASTQ_R2}" \
  --out1 "${TRIM_R1}" \
  --out2 "${TRIM_R2}" \
  --detect_adapter_for_pe \
  --correction \
  --thread "${THREADS}" \
  --html "${OUTDIR}/qc/${SAMPLE_ID}_fastp.html" \
  --json "${OUTDIR}/qc/${SAMPLE_ID}_fastp.json" \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --cut_front --cut_tail \
  2> "${OUTDIR}/qc/${SAMPLE_ID}_fastp.log"

log "  Trimming complete"


# ════════════════════════════════════════════════════════════
# STEP 3: ALIGNMENT (BWA-MEM2)
# ════════════════════════════════════════════════════════════
log "STEP 3/10: Alignment to hg38 (BWA-MEM2)"

# Index reference (run once)
if [[ ! -f "${REF_GENOME}.bwt.2bit.64" ]]; then
  log "  Indexing reference genome..."
  bwa-mem2 index "${REF_GENOME}"
fi

RAW_BAM="${OUTDIR}/align/${SAMPLE_ID}_raw.bam"
READ_GROUP="@RG\tID:${SAMPLE_ID}\tSM:${SAMPLE_ID}\tPL:ILLUMINA\tLB:lib1\tPU:unit1"

bwa-mem2 mem \
  -t "${THREADS}" \
  -R "${READ_GROUP}" \
  "${REF_GENOME}" \
  "${TRIM_R1}" "${TRIM_R2}" \
  | samtools sort -@ "${THREADS}" -o "${RAW_BAM}"

samtools index "${RAW_BAM}"
samtools flagstat "${RAW_BAM}" > "${OUTDIR}/align/${SAMPLE_ID}_flagstat.txt"
log "  Alignment complete → $(samtools view -c -F 4 ${RAW_BAM}) mapped reads"


# ════════════════════════════════════════════════════════════
# STEP 4: MARK DUPLICATES (GATK)
# ════════════════════════════════════════════════════════════
log "STEP 4/10: Mark duplicates (GATK MarkDuplicates)"

DEDUP_BAM="${OUTDIR}/align/${SAMPLE_ID}_dedup.bam"
METRICS="${OUTDIR}/align/${SAMPLE_ID}_dup_metrics.txt"

${GATK} --java-options "${JAVA_OPTS}" MarkDuplicates \
  --INPUT  "${RAW_BAM}" \
  --OUTPUT "${DEDUP_BAM}" \
  --METRICS_FILE "${METRICS}" \
  --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
  --CREATE_INDEX true \
  --VALIDATION_STRINGENCY SILENT \
  2> "${OUTDIR}/align/${SAMPLE_ID}_markdup.log"

DUP_RATE=$(grep -A 2 "ESTIMATED_LIBRARY_SIZE" "${METRICS}" | tail -1 | awk '{print $9}')
log "  Duplication rate: ${DUP_RATE}"


# ════════════════════════════════════════════════════════════
# STEP 5: BASE QUALITY SCORE RECALIBRATION (BQSR)
# ════════════════════════════════════════════════════════════
log "STEP 5/10: Base quality score recalibration (BQSR)"

BQSR_TABLE="${OUTDIR}/bqsr/${SAMPLE_ID}_recal.table"
FINAL_BAM="${OUTDIR}/bqsr/${SAMPLE_ID}_final.bam"

${GATK} --java-options "${JAVA_OPTS}" BaseRecalibrator \
  --input "${DEDUP_BAM}" \
  --reference "${REF_GENOME}" \
  --known-sites "${KNOWN_SITES_DBSNP}" \
  --known-sites "${KNOWN_SITES_MILLS}" \
  --known-sites "${KNOWN_SITES_1000G}" \
  --output "${BQSR_TABLE}" \
  2> "${OUTDIR}/bqsr/bqsr_calibrate.log"

${GATK} --java-options "${JAVA_OPTS}" ApplyBQSR \
  --input "${DEDUP_BAM}" \
  --reference "${REF_GENOME}" \
  --bqsr-recal-file "${BQSR_TABLE}" \
  --output "${FINAL_BAM}" \
  2> "${OUTDIR}/bqsr/bqsr_apply.log"

log "  BQSR complete → ${FINAL_BAM}"


# ════════════════════════════════════════════════════════════
# STEP 6: VARIANT CALLING (HaplotypeCaller — gVCF mode)
# ════════════════════════════════════════════════════════════
log "STEP 6/10: Variant calling (GATK HaplotypeCaller)"

GVCF="${OUTDIR}/gvcf/${SAMPLE_ID}.g.vcf.gz"

${GATK} --java-options "${JAVA_OPTS}" HaplotypeCaller \
  --input "${FINAL_BAM}" \
  --reference "${REF_GENOME}" \
  --output "${GVCF}" \
  --emit-ref-confidence GVCF \
  --dbsnp "${KNOWN_SITES_DBSNP}" \
  --native-pair-hmm-threads "${THREADS}" \
  --sample-name "${SAMPLE_ID}" \
  2> "${OUTDIR}/gvcf/haplotypecaller.log"

log "  gVCF created: ${GVCF}"


# ════════════════════════════════════════════════════════════
# STEP 7: GENOTYPING (single sample)
# ════════════════════════════════════════════════════════════
log "STEP 7/10: Genotyping"

RAW_VCF="${OUTDIR}/vcf/${SAMPLE_ID}_raw.vcf.gz"

${GATK} --java-options "${JAVA_OPTS}" GenotypeGVCFs \
  --reference "${REF_GENOME}" \
  --variant  "${GVCF}" \
  --output   "${RAW_VCF}" \
  --dbsnp    "${KNOWN_SITES_DBSNP}" \
  2> "${OUTDIR}/vcf/genotype.log"

TOTAL_VAR=$(bcftools stats "${RAW_VCF}" | grep "^SN" | grep "number of records" | awk '{print $NF}')
log "  Total variants called: ${TOTAL_VAR}"


# ════════════════════════════════════════════════════════════
# STEP 8: VARIANT FILTERING (hard filters for WES)
# ════════════════════════════════════════════════════════════
log "STEP 8/10: Variant filtering"

SNP_VCF="${OUTDIR}/vcf/${SAMPLE_ID}_snps_filtered.vcf.gz"
INDEL_VCF="${OUTDIR}/vcf/${SAMPLE_ID}_indels_filtered.vcf.gz"
FILTERED_VCF="${OUTDIR}/vcf/${SAMPLE_ID}_filtered.vcf.gz"

# SNP hard filters (GATK best practices)
${GATK} --java-options "${JAVA_OPTS}" SelectVariants \
  -R "${REF_GENOME}" -V "${RAW_VCF}" --select-type-to-include SNP \
  -O "${OUTDIR}/vcf/${SAMPLE_ID}_raw_snps.vcf.gz"
${GATK} --java-options "${JAVA_OPTS}" VariantFiltration \
  -R "${REF_GENOME}" -V "${OUTDIR}/vcf/${SAMPLE_ID}_raw_snps.vcf.gz" \
  --filter-expression "QD < 2.0"     --filter-name "QD2"  \
  --filter-expression "FS > 60.0"    --filter-name "FS60" \
  --filter-expression "MQ < 40.0"    --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  -O "${SNP_VCF}"

# INDEL hard filters
${GATK} --java-options "${JAVA_OPTS}" SelectVariants \
  -R "${REF_GENOME}" -V "${RAW_VCF}" --select-type-to-include INDEL \
  -O "${OUTDIR}/vcf/${SAMPLE_ID}_raw_indels.vcf.gz"
${GATK} --java-options "${JAVA_OPTS}" VariantFiltration \
  -R "${REF_GENOME}" -V "${OUTDIR}/vcf/${SAMPLE_ID}_raw_indels.vcf.gz" \
  --filter-expression "QD < 2.0"   --filter-name "QD2"  \
  --filter-expression "FS > 200.0" --filter-name "FS200" \
  --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
  -O "${INDEL_VCF}"

# Merge
${GATK} --java-options "${JAVA_OPTS}" MergeVcfs \
  -I "${SNP_VCF}" -I "${INDEL_VCF}" -O "${FILTERED_VCF}"

PASS_VAR=$(bcftools view -f PASS "${FILTERED_VCF}" | bcftools stats | grep "number of records" | awk '{print $NF}')
log "  PASS variants after filtering: ${PASS_VAR}"


# ════════════════════════════════════════════════════════════
# STEP 9: FUNCTIONAL ANNOTATION (ANNOVAR)
# ════════════════════════════════════════════════════════════
log "STEP 9/10: Annotation (ANNOVAR)"

# Convert VCF → ANNOVAR input format
ANNO_INPUT="${OUTDIR}/annotated/${SAMPLE_ID}.avinput"
ANNO_PREFIX="${OUTDIR}/annotated/${SAMPLE_ID}_anno"

convert2annovar.pl \
  -format vcf4old "${FILTERED_VCF}" \
  -outfile "${ANNO_INPUT}" \
  -withzyg --includeinfo 2>/dev/null

# Annotate with RefSeq, dbSNP, ClinVar, ExAC, CADD
table_annovar.pl "${ANNO_INPUT}" \
  /path/to/annovar/humandb/ \
  -buildver hg38 \
  -out "${ANNO_PREFIX}" \
  -remove \
  -protocol refGene,cytoBand,exac03,avsnp150,clinvar_20221231,dbnsfp42a,cadd16gt10 \
  -operation g,r,f,f,f,f,f \
  -nastring . \
  -csvout \
  2> "${OUTDIR}/annotated/annovar.log"

log "  Annotation complete → ${ANNO_PREFIX}.hg38_multianno.csv"


# ════════════════════════════════════════════════════════════
# STEP 10: AD-SPECIFIC GENE ANALYSIS
# ════════════════════════════════════════════════════════════
log "STEP 10/10: AD-specific variant extraction"

ANNO_CSV="${ANNO_PREFIX}.hg38_multianno.csv"
AD_VCF="${OUTDIR}/ad_genes/${SAMPLE_ID}_AD_genes.vcf.gz"
APOE_VCF="${OUTDIR}/ad_genes/${SAMPLE_ID}_APOE_status.vcf.gz"

# Extract variants in AD gene regions
REGION_ARGS=""
for REGION in "${AD_REGIONS[@]}"; do
  REGION_ARGS="${REGION_ARGS} --regions ${REGION}"
done
bcftools view ${REGION_ARGS} "${FILTERED_VCF}" -o "${AD_VCF}" -O z
bcftools index "${AD_VCF}"

# APOE genotyping (rs429358 and rs7412)
bcftools view \
  --regions "chr19:44908684,chr19:44908822" \
  "${FILTERED_VCF}" \
  -o "${APOE_VCF}" -O z 2>/dev/null

log "  AD gene variants: $(bcftools view -c 1 "${AD_VCF}" | bcftools stats | grep 'number of records' | awk '{print $NF}')"
log "  APOE variants extracted → ${APOE_VCF}"

# Summary stats with bcftools
bcftools stats "${AD_VCF}" > "${OUTDIR}/ad_genes/${SAMPLE_ID}_AD_stats.txt"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO\t[%GT]\n' "${AD_VCF}" \
  > "${OUTDIR}/ad_genes/${SAMPLE_ID}_AD_variants_table.txt"


# ════════════════════════════════════════════════════════════
# SUMMARY REPORT
# ════════════════════════════════════════════════════════════
REPORT="${OUTDIR}/${SAMPLE_ID}_pipeline_summary.txt"
cat > "${REPORT}" << EOF
============================================================
  DNA-seq Variant Calling Pipeline Summary
  Author: Sana Begum | Era's Lucknow Medical College
============================================================
Sample ID        : ${SAMPLE_ID}
Date             : $(date '+%Y-%m-%d %H:%M:%S')
Reference genome : hg38

Alignment Stats  : ${OUTDIR}/align/${SAMPLE_ID}_flagstat.txt
Duplication rate : ${DUP_RATE}
Total variants   : ${TOTAL_VAR}
PASS variants    : ${PASS_VAR}

Key output files :
  Final BAM       : ${FINAL_BAM}
  Filtered VCF    : ${FILTERED_VCF}
  Annotated CSV   : ${ANNO_CSV}
  AD Gene VCF     : ${AD_VCF}
  APOE Status     : ${APOE_VCF}

AD gene regions analysed:
  APOE, PSEN1, PSEN2, APP, TREM2, SORL1, BIN1, CLU, CR1
============================================================
EOF

log "Pipeline complete. Summary: ${REPORT}"
