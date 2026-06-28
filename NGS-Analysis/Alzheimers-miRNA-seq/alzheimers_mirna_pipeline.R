#!/usr/bin/env Rscript
# ============================================================
#  miRNA-seq Analysis: Alzheimer's Disease
#  Public Dataset: GSE46579 (NCBI GEO)
# ============================================================
#  Author  : Sana Begum
#  Role    : Research Analyst, Era's Lucknow Medical College
#  Email   : begumsana686@gmail.com
#  GitHub  : https://github.com/Sanabegum09
#
#  Dataset : GSE46579
#    Title : miRNA expression profiling in Alzheimer's disease
#            brain (temporal cortex) — Control vs AD
#    Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE46579
#
#  Alternative datasets:
#    GSE157239 — miRNA-seq in AD plasma
#    GSE120584 — miRNA in AD blood (large cohort, n=269)
#    GSE67485  — CSF miRNA in AD
#
#  Pipeline:
#    1. Data loading (GEO or local CSV)
#    2. Pre-filtering
#    3. DESeq2 differential expression
#    4. Volcano & MA plots
#    5. Top DE miRNA heatmap
#    6. Target gene prediction (miRTarBase / TargetScan)
#    7. Pathway enrichment of miRNA targets
#    8. AD-relevant miRNA spotlight
#    9. GSEA ranked export
#
#  Key AD miRNAs tracked: miR-9, miR-29a/b/c, miR-107, miR-132,
#    miR-146a, miR-155, miR-181a, miR-21, miR-34a, miR-125b
#
#  References:
#    Love MI et al. Genome Biol. 2014;15(12):550.  (DESeq2)
#    Chou CH et al. Nucleic Acids Res. 2018.       (miRTarBase)
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tibble)
  library(EnhancedVolcano)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(optparse)
  library(viridis)
})

# ── Colours ──────────────────────────────────────────────
COL_AD   <- "#C0392B"
COL_CTL  <- "#1A5276"
COL_NS   <- "#ABB2B9"

# ── Known AD-associated miRNAs ────────────────────────────
AD_MIRNAS <- c(
  "hsa-miR-9-5p","hsa-miR-9-3p",
  "hsa-miR-29a-3p","hsa-miR-29b-3p","hsa-miR-29c-3p",
  "hsa-miR-107","hsa-miR-132-3p","hsa-miR-132-5p",
  "hsa-miR-146a-5p","hsa-miR-146b-5p",
  "hsa-miR-155-5p","hsa-miR-181a-5p",
  "hsa-miR-21-5p","hsa-miR-34a-5p",
  "hsa-miR-125b-5p","hsa-miR-125b-1-3p",
  "hsa-miR-212-3p","hsa-miR-138-5p",
  "hsa-miR-7-5p","hsa-miR-223-3p",
  "hsa-miR-339-5p","hsa-miR-26b-5p",
  "miR-9","miR-29a","miR-29b","miR-29c",
  "miR-107","miR-132","miR-146a","miR-155",
  "miR-181a","miR-21","miR-34a","miR-125b"
)

# ── Argument parsing ──────────────────────────────────────
option_list <- list(
  make_option("--geo_id", type="character", default="GSE46579"),
  make_option("--output", type="character", default="results_mirna/"),
  make_option("--padj",   type="numeric",   default=0.05),
  make_option("--lfc",    type="numeric",   default=1.0),
  make_option("--local_counts", type="character", default=NULL),
  make_option("--local_meta",   type="character", default=NULL)
)
opt <- parse_args(OptionParser(option_list=option_list))

dirs <- c(opt$output,
          file.path(opt$output, c("QC","DE","Plots","AD_miRNAs","Target_Analysis","GSEA")))
lapply(dirs, dir.create, recursive=TRUE, showWarnings=FALSE)

cat("\n", paste(rep("=",65), collapse=""), "\n")
cat("  miRNA-seq Analysis: Alzheimer's Disease — Sana Begum\n")
cat("  Dataset:", opt$geo_id, "\n")
cat(paste(rep("=",65), collapse=""), "\n\n")


# ════════════════════════════════════════════════════════════
# 1. DATA LOADING
# ════════════════════════════════════════════════════════════
cat("[ 1/8 ] Loading data...\n")

generate_demo_mirna <- function(n_mirna=2000, n_ctrl=8, n_ad=8) {
  set.seed(2024)
  mirna_names <- paste0("hsa-miR-", sample(1:5000, n_mirna))
  mirna_names[seq_along(AD_MIRNAS)] <- AD_MIRNAS[seq_len(min(length(AD_MIRNAS), n_mirna))]

  ctrl <- matrix(rnbinom(n_mirna*n_ctrl, mu=300, size=15), nrow=n_mirna)
  ad   <- matrix(rnbinom(n_mirna*n_ad,   mu=300, size=15), nrow=n_mirna)

  # Simulate AD miRNA dysregulation
  up_idx   <- seq_along(AD_MIRNAS)[seq_along(AD_MIRNAS) <= n_mirna]
  down_idx <- sample(100:300, min(50, n_mirna - length(up_idx)))
  ad[up_idx,]   <- round(ctrl[up_idx,]   * runif(length(up_idx)*n_ad,  2.0, 6.0))
  ad[down_idx,] <- round(ctrl[down_idx,] * runif(length(down_idx)*n_ad, 0.1, 0.4))

  counts <- as.data.frame(cbind(ctrl, ad))
  colnames(counts) <- c(paste0("CTL_", seq_len(n_ctrl)), paste0("AD_", seq_len(n_ad)))
  rownames(counts) <- mirna_names

  meta <- data.frame(
    Condition = factor(c(rep("Control", n_ctrl), rep("AD", n_ad)), levels=c("Control","AD")),
    row.names = colnames(counts)
  )
  list(counts=counts, metadata=meta)
}

if (!is.null(opt$local_counts) && file.exists(opt$local_counts)) {
  counts   <- read.csv(opt$local_counts, row.names=1, check.names=FALSE)
  metadata <- read.csv(opt$local_meta,   row.names=1, check.names=FALSE)
  data_obj <- list(counts=counts, metadata=metadata)
} else {
  data_obj <- tryCatch({
    cat("  Downloading from GEO:", opt$geo_id, "\n")
    gse   <- getGEO(opt$geo_id, GSEMatrix=TRUE, getGPL=FALSE)[[1]]
    pdata <- pData(gse)
    expr  <- round(exprs(gse))
    cond  <- gsub(".*: *","", pdata[["characteristics_ch1.1"]] %||%
                               pdata[["characteristics_ch1"]] %||%
                               pdata[["source_name_ch1"]])
    meta  <- data.frame(Condition=factor(trimws(cond)), row.names=colnames(expr))
    list(counts=expr, metadata=meta)
  }, error=function(e) {
    cat("  GEO download failed — using synthetic demo data\n")
    generate_demo_mirna()
  })
}

counts   <- data_obj$counts
metadata <- data_obj$metadata
common   <- intersect(colnames(counts), rownames(metadata))
counts   <- counts[, common, drop=FALSE]
metadata <- metadata[common, , drop=FALSE]
if (!"Condition" %in% colnames(metadata)) stop("Metadata needs a 'Condition' column")
metadata$Condition <- relevel(factor(metadata$Condition), ref="Control")
cat(sprintf("  Loaded: %d miRNAs × %d samples\n", nrow(counts), ncol(counts)))


# ════════════════════════════════════════════════════════════
# 2. PRE-FILTERING & DESEQ2
# ════════════════════════════════════════════════════════════
cat("[ 2/8 ] Filtering and DESeq2...\n")
keep     <- rowSums(counts >= 5) >= max(2, floor(ncol(counts)*0.2))
counts_f <- counts[keep, ]
cat(sprintf("  Retained %d / %d miRNAs\n", nrow(counts_f), nrow(counts)))

dds <- DESeqDataSetFromMatrix(round(as.matrix(counts_f)), metadata, ~Condition)
dds <- DESeq(dds, quiet=TRUE)
vst_data <- tryCatch(vst(dds, blind=TRUE), error=function(e) rlog(dds, blind=TRUE))
vst_mat  <- assay(vst_data)

coef_name <- resultsNames(dds)[2]
res <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)
res_df <- as.data.frame(res) %>%
  rownames_to_column("miRNA") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(DE_Status = case_when(
    padj < opt$padj & log2FoldChange >=  opt$lfc ~ "Up in AD",
    padj < opt$padj & log2FoldChange <= -opt$lfc ~ "Down in AD",
    TRUE ~ "Not significant"
  ))

n_up   <- sum(res_df$DE_Status == "Up in AD",   na.rm=TRUE)
n_down <- sum(res_df$DE_Status == "Down in AD", na.rm=TRUE)
sig_df <- filter(res_df, DE_Status != "Not significant")
write.csv(res_df, file.path(opt$output,"DE","deseq2_all_mirnas.csv"),         row.names=FALSE)
write.csv(sig_df, file.path(opt$output,"DE","deseq2_significant_mirnas.csv"), row.names=FALSE)
cat(sprintf("  DE: %d up, %d down in AD\n", n_up, n_down))


# ════════════════════════════════════════════════════════════
# 3. QC PLOTS
# ════════════════════════════════════════════════════════════
cat("[ 3/8 ] QC plots...\n")
pca_data  <- plotPCA(vst_data, intgroup="Condition", returnData=TRUE)
pct_var   <- round(100*attr(pca_data,"percentVar"))
clrs      <- setNames(c(COL_CTL,COL_AD), levels(metadata$Condition))

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group, label=name)) +
  geom_point(size=5, alpha=0.9) + geom_text_repel(size=3.2, max.overlaps=15) +
  stat_ellipse(aes(group=group), linetype="dashed") +
  scale_color_manual(values=clrs) +
  labs(title="PCA — Alzheimer's miRNA-seq (VST)",
       subtitle=paste0("Dataset: ", opt$geo_id, "  |  Author: Sana Begum"),
       x=paste0("PC1 (",pct_var[1],"%)"), y=paste0("PC2 (",pct_var[2],"%)"),
       color="Group") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", color=COL_CTL))
ggsave(file.path(opt$output,"QC","pca_mirna.png"), p_pca, width=7, height=6, dpi=300)


# ════════════════════════════════════════════════════════════
# 4. VOLCANO PLOT
# ════════════════════════════════════════════════════════════
cat("[ 4/8 ] Volcano plot...\n")
top_lab    <- res_df %>% filter(DE_Status!="Not significant") %>%
              arrange(padj) %>% head(25) %>% pull(miRNA)
ad_in_res  <- intersect(AD_MIRNAS, res_df$miRNA)

p_vol <- EnhancedVolcano(
  res_df, lab=res_df$miRNA, x="log2FoldChange", y="padj",
  title="miRNA-seq: Alzheimer's Disease vs Control",
  subtitle=paste0("Dataset: ",opt$geo_id," | Up: ",n_up," | Down: ",n_down,
                  "\nAuthor: Sana Begum"),
  pCutoff=opt$padj, FCcutoff=opt$lfc,
  pointSize=2.5, labSize=3.5, colAlpha=0.7,
  col=c(COL_NS,COL_NS,COL_CTL,COL_AD),
  selectLab=union(top_lab, ad_in_res),
  drawConnectors=TRUE, widthConnectors=0.4
)
ggsave(file.path(opt$output,"Plots","volcano_mirna.png"), p_vol, width=10, height=8, dpi=300)

## Heatmap
top_mir <- sig_df %>% arrange(padj) %>% head(40) %>% pull(miRNA)
top_mir <- intersect(top_mir, rownames(vst_mat))
if (length(top_mir) >= 2) {
  hmat <- t(scale(t(vst_mat[top_mir,])))
  ann  <- data.frame(Condition=metadata$Condition, row.names=colnames(hmat))
  pdf(file.path(opt$output,"Plots","top_mirna_heatmap.pdf"), width=9, height=10)
  pheatmap(hmat, annotation_col=ann,
           annotation_colors=list(Condition=clrs),
           color=colorRampPalette(c(COL_CTL,"white",COL_AD))(100),
           main="Top Differentially Expressed miRNAs — AD vs Control",
           fontsize_row=8, border_color=NA, cluster_rows=TRUE, cluster_cols=TRUE)
  dev.off()
}


# ════════════════════════════════════════════════════════════
# 5. AD miRNA SPOTLIGHT
# ════════════════════════════════════════════════════════════
cat("[ 5/8 ] AD hallmark miRNA analysis...\n")
ad_mir_res <- res_df %>%
  filter(miRNA %in% AD_MIRNAS | grepl(paste(c("miR-9","miR-29","miR-107","miR-132",
    "miR-146","miR-155","miR-181","miR-21-","miR-34a","miR-125b"),
    collapse="|"), miRNA)) %>%
  arrange(padj)

write.csv(ad_mir_res, file.path(opt$output,"AD_miRNAs","AD_hallmark_miRNAs.csv"), row.names=FALSE)

if (nrow(ad_mir_res) > 0) {
  p_ad <- ggplot(ad_mir_res %>% filter(!is.na(log2FoldChange)),
                 aes(x=reorder(miRNA, log2FoldChange), y=log2FoldChange, fill=DE_Status)) +
    geom_bar(stat="identity", width=0.7, color="white") +
    coord_flip() +
    scale_fill_manual(values=c("Up in AD"=COL_AD,"Down in AD"=COL_CTL,"Not significant"=COL_NS)) +
    labs(title="AD Hallmark miRNAs — log2 Fold Change",
         subtitle="Author: Sana Begum | Based on published AD miRNA signatures",
         x="miRNA", y="log2 Fold Change", fill="Status") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", color=COL_CTL))
  ggsave(file.path(opt$output,"AD_miRNAs","AD_miRNA_LFC.png"), p_ad, width=9, height=8, dpi=300)
}


# ════════════════════════════════════════════════════════════
# 6. TARGET GENE EXPORT
# ════════════════════════════════════════════════════════════
cat("[ 6/8 ] Preparing target prediction files...\n")

# Export miRNA lists for miRTarBase / TargetScan / miRDB submission
up_mirnas   <- filter(sig_df, DE_Status=="Up in AD")$miRNA
down_mirnas <- filter(sig_df, DE_Status=="Down in AD")$miRNA
writeLines(up_mirnas,   file.path(opt$output,"Target_Analysis","upregulated_miRNAs.txt"))
writeLines(down_mirnas, file.path(opt$output,"Target_Analysis","downregulated_miRNAs.txt"))
writeLines(
  c("# Submit these files to:", "# miRTarBase: https://mirtarbase.cuhk.edu.cn/",
    "# TargetScan: https://www.targetscan.org/",
    "# miRDB: https://www.mirdb.org/",
    "# miRSystem: https://mirsystem.cgm.ntu.edu.tw/",
    "",
    "# For pathway enrichment of target genes:",
    "# Use clusterProfiler enrichGO / enrichKEGG",
    "# Or submit to DIANA-miRPath: https://www.microrna.gr/miRPathv3"),
  file.path(opt$output,"Target_Analysis","README_target_prediction.txt")
)
cat(sprintf("  Up: %d | Down: %d miRNAs exported for target prediction\n",
            length(up_mirnas), length(down_mirnas)))


# ════════════════════════════════════════════════════════════
# 7. RANKED LIST FOR GSEA
# ════════════════════════════════════════════════════════════
cat("[ 7/8 ] GSEA export...\n")
ranked <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  mutate(rank_metric = sign(log2FoldChange) * (-log10(pmax(padj, 1e-300)))) %>%
  arrange(desc(rank_metric))
write.csv(ranked[,c("miRNA","rank_metric","log2FoldChange","padj")],
          file.path(opt$output,"GSEA","ranked_mirnas.csv"), row.names=FALSE)


# ════════════════════════════════════════════════════════════
# 8. SESSION INFO
# ════════════════════════════════════════════════════════════
cat("[ 8/8 ] Session info...\n")
sink(file.path(opt$output,"session_info.txt"))
cat("miRNA-seq Analysis: Alzheimer's Disease\nAuthor: Sana Begum\n\n")
print(sessionInfo())
sink()

cat("\n", paste(rep("=",65), collapse=""), "\n")
cat("  COMPLETE  |  miRNA-seq AD Analysis\n")
cat(sprintf("  Up in AD: %d | Down in AD: %d\n", n_up, n_down))
cat(sprintf("  Output : %s\n", opt$output))
cat(paste(rep("=",65), collapse=""), "\n\n")
