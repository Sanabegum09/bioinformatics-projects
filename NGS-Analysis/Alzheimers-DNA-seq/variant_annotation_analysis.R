#!/usr/bin/env Rscript
# ============================================================
#  Variant Annotation & AD Gene Analysis (post DNA-seq)
#  Companion to: alzheimers_dnaseq_pipeline.sh
# ============================================================
#  Author : Sana Begum | begumsana686@gmail.com
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr)
  library(viridis); library(optparse)
  library(RColorBrewer)
})

option_list <- list(
  make_option("--annovar_csv", type="character", help="ANNOVAR multianno CSV"),
  make_option("--output",      type="character", default="results_dnaseq/annotation/")
)
opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)

BLUE <- "#1A5276"; RED <- "#C0392B"; GREEN <- "#1E8449"

# ── AD Gene Set ───────────────────────────────────────────
AD_GENES <- c("APOE","PSEN1","PSEN2","APP","TREM2","SORL1","BIN1",
              "CLU","CR1","ABCA7","PICALM","MS4A6A","CD33",
              "EPHA1","CD2AP","ADAM10","BACE1","SNCA","GBA")

AD_KEY_VARIANTS <- data.frame(
  rsID   = c("rs429358","rs7412","rs75932628","rs63750847"),
  Gene   = c("APOE","APOE","TREM2","APP"),
  Change = c("C>T (ε4)","T>C (ε2)","R47H","V717I"),
  Effect = c("Major AD risk","Protective","Risk, LOAD","Early-onset AD"),
  stringsAsFactors=FALSE
)

cat("\n====================================================\n")
cat("  Variant Annotation & AD Gene Analysis — Sana Begum\n")
cat("====================================================\n\n")

if (is.null(opt$annovar_csv) || !file.exists(opt$annovar_csv)) {
  cat("NOTE: ANNOVAR CSV not provided — generating demonstration output.\n")
  cat("Provide --annovar_csv path/to/annovar_multianno.csv for real analysis.\n\n")

  # Demo variant data
  set.seed(42)
  n_var <- 2000
  demo_df <- data.frame(
    Chr      = sample(paste0("chr", c(1:22,"X")), n_var, replace=TRUE),
    Start    = sample(1e6:200e6, n_var),
    Ref      = sample(c("A","T","G","C"), n_var, replace=TRUE),
    Alt      = sample(c("A","T","G","C","ATG","ATCG"), n_var, replace=TRUE),
    Func     = sample(c("exonic","intronic","splicing","UTR3","UTR5","intergenic"),
                      n_var, prob=c(.05,.50,.02,.05,.03,.35), replace=TRUE),
    ExonicFunc = sample(c("synonymous SNV","nonsynonymous SNV","stopgain",
                           "frameshift insertion","frameshift deletion","NA"),
                        n_var, prob=c(.3,.25,.02,.01,.01,.41), replace=TRUE),
    Gene     = c(sample(AD_GENES, 50, replace=TRUE),
                 paste0("GENE_", sample(1:5000, n_var-50, replace=TRUE))),
    ClinVar  = sample(c("Pathogenic","Likely pathogenic","Benign",
                        "Uncertain significance","NA"),
                      n_var, prob=c(.01,.02,.3,.07,.6), replace=TRUE),
    CADD_Phred = round(runif(n_var, 0, 40), 1),
    ExAC_ALL  = round(runif(n_var, 0, 0.5), 4),
    stringsAsFactors=FALSE
  )
  # Spike in key AD variants
  demo_df[1:nrow(AD_KEY_VARIANTS), "Gene"] <- AD_KEY_VARIANTS$Gene
  df <- demo_df
} else {
  df <- read.csv(opt$annovar_csv, stringsAsFactors=FALSE)
  cat(sprintf("Loaded: %d variants from ANNOVAR annotation\n", nrow(df)))
}

# ── 1. Variant classification summary ────────────────────
if ("Func" %in% colnames(df)) {
  func_count <- df %>% count(Func, sort=TRUE)
  p1 <- ggplot(func_count, aes(x=reorder(Func, n), y=n, fill=Func)) +
    geom_bar(stat="identity", color="white") +
    coord_flip() + scale_fill_viridis_d(option="C") +
    labs(title="Variant Distribution by Genomic Region",
         subtitle="Author: Sana Begum | ANNOVAR annotation",
         x="Region", y="Count") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", color=BLUE), legend.position="none")
  ggsave(file.path(opt$output,"variant_function_distribution.png"), p1, width=7, height=5, dpi=200)
}

# ── 2. Exonic variant breakdown ───────────────────────────
if ("ExonicFunc" %in% colnames(df)) {
  exonic <- df %>% filter(Func=="exonic", !is.na(ExonicFunc), ExonicFunc != "NA") %>%
    count(ExonicFunc, sort=TRUE)
  if (nrow(exonic) > 0) {
    p2 <- ggplot(exonic, aes(x="", y=n, fill=ExonicFunc)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y") +
      scale_fill_manual(values=colorRampPalette(brewer.pal(8,"Set2"))(nrow(exonic))) +
      labs(title="Exonic Variant Types", fill="Type",
           subtitle="Author: Sana Begum") +
      theme_void(base_size=12) +
      theme(plot.title=element_text(face="bold", color=BLUE, hjust=0.5))
    ggsave(file.path(opt$output,"exonic_variant_types.png"), p2, width=6, height=6, dpi=200)
  }
}

# ── 3. CADD score distribution ────────────────────────────
if ("CADD_Phred" %in% colnames(df)) {
  p3 <- ggplot(df %>% filter(!is.na(CADD_Phred)),
               aes(x=CADD_Phred)) +
    geom_histogram(bins=40, fill=BLUE, color="white", alpha=0.85) +
    geom_vline(xintercept=20, color=RED,    linestyle="dashed", linewidth=0.8) +
    geom_vline(xintercept=30, color=GREEN,  linestyle="dashed", linewidth=0.8) +
    annotate("text", x=21, y=Inf, label=" CADD≥20\n(predicted deleterious)",
             vjust=1.5, hjust=0, size=3.5, color=RED) +
    annotate("text", x=31, y=Inf, label=" CADD≥30\n(top 0.1%)", 
             vjust=1.5, hjust=0, size=3.5, color=GREEN) +
    labs(title="CADD Pathogenicity Score Distribution",
         subtitle="Author: Sana Begum | CADD ≥ 20: predicted deleterious",
         x="CADD Phred Score", y="Number of Variants") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", color=BLUE))
  ggsave(file.path(opt$output,"cadd_distribution.png"), p3, width=8, height=5, dpi=200)
}

# ── 4. AD gene variant summary ────────────────────────────
if ("Gene" %in% colnames(df)) {
  ad_vars <- df %>%
    filter(Gene %in% AD_GENES) %>%
    mutate(Damaging = CADD_Phred >= 20) %>%
    group_by(Gene) %>%
    summarise(
      Total_Variants  = n(),
      Damaging_CADD20 = sum(Damaging, na.rm=TRUE),
      ClinVar_Path    = sum(grepl("Pathogenic", ClinVar, ignore.case=TRUE), na.rm=TRUE),
      .groups="drop"
    ) %>%
    arrange(desc(Damaging_CADD20))

  write.csv(ad_vars, file.path(opt$output,"AD_gene_variant_summary.csv"), row.names=FALSE)

  if (nrow(ad_vars) > 0) {
    p4 <- ad_vars %>%
      pivot_longer(cols=c(Total_Variants,Damaging_CADD20,ClinVar_Path),
                   names_to="Category", values_to="Count") %>%
      ggplot(aes(x=reorder(Gene, Count), y=Count, fill=Category)) +
      geom_bar(stat="identity", position="dodge", color="white") +
      coord_flip() +
      scale_fill_manual(values=c(
        "Total_Variants"=BLUE, "Damaging_CADD20"=RED, "ClinVar_Path"="#8E44AD")) +
      labs(title="Alzheimer's Disease Gene — Variant Summary",
           subtitle="Author: Sana Begum | CADD ≥ 20 = predicted damaging",
           x="Gene", y="Number of Variants", fill="Category") +
      theme_bw(base_size=12) +
      theme(plot.title=element_text(face="bold", color=BLUE))
    ggsave(file.path(opt$output,"AD_gene_variants.png"), p4, width=9, height=6, dpi=200)
  }
}

# ── 5. Key AD variant table ───────────────────────────────
cat("\nKey AD Variant Reference:\n")
cat(sprintf("%-15s %-8s %-20s %-35s\n", "rsID","Gene","Change","Clinical Effect"))
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(AD_KEY_VARIANTS))) {
  cat(sprintf("%-15s %-8s %-20s %-35s\n",
              AD_KEY_VARIANTS$rsID[i], AD_KEY_VARIANTS$Gene[i],
              AD_KEY_VARIANTS$Change[i], AD_KEY_VARIANTS$Effect[i]))
}
cat("\n")
write.csv(AD_KEY_VARIANTS, file.path(opt$output,"key_AD_variants_reference.csv"), row.names=FALSE)

cat(sprintf("Outputs saved to: %s\n\n", opt$output))
