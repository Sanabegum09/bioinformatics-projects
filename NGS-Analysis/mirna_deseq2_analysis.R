#!/usr/bin/env Rscript
# ============================================================
#  miRNA-seq Differential Expression Analysis Pipeline
# ============================================================
#  Author  : Sana Begum
#  Role    : Research Analyst, Era's Lucknow Medical College
#  Email   : begumsana686@gmail.com
#  GitHub  : https://github.com/Sanabegum09
#
#  Description:
#  Complete differential expression analysis pipeline for
#  miRNA-seq data using DESeq2. Includes:
#    - Read count loading & QC
#    - Normalisation (VST/rlog)
#    - PCA & sample clustering
#    - DESeq2 differential expression
#    - Volcano plot & MA plot
#    - Top DE miRNA heatmap
#    - Target gene prediction prep
#    - Pathway enrichment (clusterProfiler)
#
#  Input:
#    - count_matrix.csv : raw read counts (rows = miRNAs, cols = samples)
#    - metadata.csv     : sample metadata (SampleID, Condition)
#
#  Usage (command line):
#    Rscript mirna_deseq2_analysis.R \
#        --counts count_matrix.csv \
#        --metadata metadata.csv \
#        --output results/ \
#        --condition Condition \
#        --reference Normal
#
#  References:
#    Love MI et al. Genome Biol. 2014;15(12):550.  (DESeq2)
#    Yu G et al. OMICS. 2012;16(5):284-287.        (clusterProfiler)
# ============================================================

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(EnhancedVolcano)
})

# ─── COLOUR PALETTE ─────────────────────────────────────────
BLUE   <- "#1A5276"
RED    <- "#C0392B"
GREY   <- "#AAB7B8"
ORANGE <- "#E67E22"

# ─── ARGUMENT PARSING ───────────────────────────────────────
option_list <- list(
  make_option("--counts",    type="character", help="Raw count matrix CSV"),
  make_option("--metadata",  type="character", help="Sample metadata CSV"),
  make_option("--output",    type="character", default="results/", help="Output directory"),
  make_option("--condition", type="character", default="Condition", help="Metadata column for comparison"),
  make_option("--reference", type="character", default="Normal",    help="Reference/control group name"),
  make_option("--padj",      type="numeric",   default=0.05, help="Adjusted p-value threshold"),
  make_option("--lfc",       type="numeric",   default=1.0,  help="|log2FoldChange| threshold")
)
opt <- parse_args(OptionParser(option_list=option_list))

# ─── OUTPUT DIRECTORIES ─────────────────────────────────────
dir.create(opt$output, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(opt$output, "QC"),   showWarnings=FALSE)
dir.create(file.path(opt$output, "DE"),   showWarnings=FALSE)
dir.create(file.path(opt$output, "Plots"),showWarnings=FALSE)

cat("\n", paste(rep("=",60), collapse=""), "\n")
cat("  miRNA-seq DE Analysis Pipeline — Sana Begum\n")
cat(paste(rep("=",60), collapse=""), "\n\n")


# ════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ════════════════════════════════════════════════════════════
cat("[ 1/7 ] Loading count matrix and metadata...\n")

counts   <- read.csv(opt$counts, row.names=1, check.names=FALSE)
metadata <- read.csv(opt$metadata, row.names=1, check.names=FALSE)

# Align samples
common_samples <- intersect(colnames(counts), rownames(metadata))
if (length(common_samples) == 0) stop("No common samples between counts and metadata!")
counts   <- counts[, common_samples]
metadata <- metadata[common_samples, , drop=FALSE]
metadata[[opt$condition]] <- relevel(factor(metadata[[opt$condition]]), ref=opt$reference)

cat(sprintf("  Loaded: %d miRNAs x %d samples\n", nrow(counts), ncol(counts)))
cat(sprintf("  Condition: '%s' (reference: '%s')\n", opt$condition, opt$reference))


# ════════════════════════════════════════════════════════════
# 2. PRE-FILTERING
# ════════════════════════════════════════════════════════════
cat("[ 2/7 ] Pre-filtering low-count miRNAs...\n")

# Keep miRNAs with at least 10 counts in at least 2 samples
keep <- rowSums(counts >= 10) >= 2
counts_filt <- counts[keep, ]
cat(sprintf("  Retained %d / %d miRNAs after filtering\n", nrow(counts_filt), nrow(counts)))


# ════════════════════════════════════════════════════════════
# 3. DESEQ2 ANALYSIS
# ════════════════════════════════════════════════════════════
cat("[ 3/7 ] Running DESeq2 differential expression analysis...\n")

formula <- as.formula(paste0("~ ", opt$condition))
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_filt),
  colData   = metadata,
  design    = formula
)

dds <- DESeq(dds, parallel=FALSE, quiet=TRUE)

# Variance stabilising transformation for visualisation
vst_data <- tryCatch(vst(dds, blind=TRUE), error=function(e) rlog(dds, blind=TRUE))

# Results with LFC shrinkage (apeglm)
res_name <- resultsNames(dds)[2]   # second coefficient = treatment vs reference
res <- lfcShrink(dds, coef=res_name, type="apeglm", quiet=TRUE)

res_df <- as.data.frame(res) %>%
  rownames_to_column("miRNA") %>%
  arrange(padj, desc(abs(log2FoldChange)))

# Classify miRNAs
res_df <- res_df %>%
  mutate(
    DE_Status = case_when(
      padj < opt$padj & log2FoldChange >= opt$lfc  ~ "Up-regulated",
      padj < opt$padj & log2FoldChange <= -opt$lfc ~ "Down-regulated",
      TRUE ~ "Not significant"
    )
  )

n_up   <- sum(res_df$DE_Status == "Up-regulated",   na.rm=TRUE)
n_down <- sum(res_df$DE_Status == "Down-regulated", na.rm=TRUE)
cat(sprintf("  DE results: %d up-regulated, %d down-regulated (padj<%.2f, |LFC|>%.1f)\n",
            n_up, n_down, opt$padj, opt$lfc))

# Save all results and significant subset
write.csv(res_df, file.path(opt$output, "DE", "deseq2_all_results.csv"), row.names=FALSE)
sig <- filter(res_df, DE_Status != "Not significant")
write.csv(sig, file.path(opt$output, "DE", "deseq2_significant_miRNAs.csv"), row.names=FALSE)
cat(sprintf("  Results saved: %d significant miRNAs\n", nrow(sig)))


# ════════════════════════════════════════════════════════════
# 4. QC PLOTS
# ════════════════════════════════════════════════════════════
cat("[ 4/7 ] Generating QC plots...\n")

## 4a. PCA
pca_data  <- plotPCA(vst_data, intgroup=opt$condition, returnData=TRUE)
pct_var   <- round(100 * attr(pca_data, "percentVar"))

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group, label=name)) +
  geom_point(size=4, alpha=0.85) +
  geom_text_repel(size=3.5, max.overlaps=15) +
  scale_color_manual(values=c(BLUE, RED, ORANGE, GREY)) +
  labs(
    title    = "PCA of Samples (VST-normalised)",
    subtitle = paste0("Author: Sana Begum | miRNA-seq analysis"),
    x        = paste0("PC1: ", pct_var[1], "% variance"),
    y        = paste0("PC2: ", pct_var[2], "% variance"),
    color    = opt$condition
  ) +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", colour=BLUE))
ggsave(file.path(opt$output, "QC", "pca_plot.pdf"), p_pca, width=7, height=5.5)
ggsave(file.path(opt$output, "QC", "pca_plot.png"), p_pca, width=7, height=5.5, dpi=300)

## 4b. Sample distance heatmap
vst_mat <- assay(vst_data)
sample_dists <- dist(t(vst_mat))
dist_mat     <- as.matrix(sample_dists)
colors_heat  <- colorRampPalette(rev(brewer.pal(9,"Blues")))(100)
pdf(file.path(opt$output, "QC", "sample_distance_heatmap.pdf"), width=7, height=6)
pheatmap(dist_mat, col=colors_heat,
         clustering_distance_rows=sample_dists,
         clustering_distance_cols=sample_dists,
         main="Sample Distance Matrix (VST)")
dev.off()


# ════════════════════════════════════════════════════════════
# 5. DE VISUALISATION
# ════════════════════════════════════════════════════════════
cat("[ 5/7 ] Generating DE plots...\n")

## 5a. Volcano plot
top_labels <- res_df %>%
  filter(DE_Status != "Not significant") %>%
  arrange(padj) %>%
  head(20)

p_volcano <- EnhancedVolcano(
  res_df,
  lab          = res_df$miRNA,
  x            = "log2FoldChange",
  y            = "padj",
  title        = "Volcano Plot — miRNA Differential Expression",
  subtitle     = paste0("Author: Sana Begum  |  Up: ", n_up, "  |  Down: ", n_down),
  pCutoff      = opt$padj,
  FCcutoff     = opt$lfc,
  pointSize    = 2.5,
  labSize      = 3.5,
  colAlpha     = 0.75,
  legendPosition = "right",
  selectLab    = top_labels$miRNA,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  col          = c(GREY, GREY, BLUE, RED)
)
ggsave(file.path(opt$output, "Plots", "volcano_plot.pdf"), p_volcano, width=9, height=7)
ggsave(file.path(opt$output, "Plots", "volcano_plot.png"), p_volcano, width=9, height=7, dpi=300)

## 5b. MA plot
p_ma <- ggplot(res_df, aes(x=log10(baseMean+1), y=log2FoldChange, color=DE_Status)) +
  geom_point(alpha=0.6, size=1.5) +
  geom_hline(yintercept=c(-opt$lfc, opt$lfc), linetype="dashed", color="black", linewidth=0.5) +
  scale_color_manual(values=c("Up-regulated"=RED, "Down-regulated"=BLUE, "Not significant"=GREY)) +
  labs(
    title = "MA Plot — Differential Expression",
    x     = "log10(mean normalised count + 1)",
    y     = "log2 Fold Change",
    color = "Status"
  ) +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", colour=BLUE))
ggsave(file.path(opt$output, "Plots", "ma_plot.pdf"), p_ma, width=7, height=5)
ggsave(file.path(opt$output, "Plots", "ma_plot.png"), p_ma, width=7, height=5, dpi=300)

## 5c. Top DE miRNA heatmap
top_n <- 40
top_mirnas <- sig %>% arrange(padj) %>% head(top_n) %>% pull(miRNA)
if (length(top_mirnas) >= 2) {
  heat_mat <- vst_mat[top_mirnas[top_mirnas %in% rownames(vst_mat)], ]
  heat_mat_scaled <- t(scale(t(heat_mat)))
  ann_col <- data.frame(Condition=metadata[[opt$condition]], row.names=rownames(metadata))
  ann_colors <- list(Condition=setNames(c(BLUE, RED), levels(metadata[[opt$condition]])))
  pdf(file.path(opt$output, "Plots", "top_mirna_heatmap.pdf"), width=8, height=9)
  pheatmap(
    heat_mat_scaled,
    annotation_col  = ann_col,
    annotation_colors = ann_colors,
    show_colnames   = TRUE,
    show_rownames   = TRUE,
    fontsize_row    = 8,
    cluster_rows    = TRUE,
    cluster_cols    = TRUE,
    color           = colorRampPalette(c(BLUE, "white", RED))(100),
    main            = paste0("Top ", nrow(heat_mat_scaled), " Differentially Expressed miRNAs")
  )
  dev.off()
}


# ════════════════════════════════════════════════════════════
# 6. TARGET GENE ANALYSIS PREPARATION
# ════════════════════════════════════════════════════════════
cat("[ 6/7 ] Preparing target gene analysis output...\n")

# Save up- and down-regulated miRNA lists for use with miRTarBase / TargetScan
up_mirnas   <- filter(sig, DE_Status == "Up-regulated")$miRNA
down_mirnas <- filter(sig, DE_Status == "Down-regulated")$miRNA

writeLines(up_mirnas,   file.path(opt$output, "DE", "upregulated_miRNAs.txt"))
writeLines(down_mirnas, file.path(opt$output, "DE", "downregulated_miRNAs.txt"))

# Export ranked list for GSEA
ranked <- res_df %>%
  filter(!is.na(log2FoldChange)) %>%
  select(miRNA, log2FoldChange) %>%
  arrange(desc(log2FoldChange))
write.csv(ranked, file.path(opt$output, "DE", "ranked_mirnas_for_gsea.csv"), row.names=FALSE)

cat(sprintf("  Up-regulated list: %d miRNAs\n  Down-regulated list: %d miRNAs\n",
            length(up_mirnas), length(down_mirnas)))


# ════════════════════════════════════════════════════════════
# 7. SESSION INFO & SUMMARY
# ════════════════════════════════════════════════════════════
cat("[ 7/7 ] Writing session info...\n")
sink(file.path(opt$output, "session_info.txt"))
print(sessionInfo())
sink()

cat("\n", paste(rep("=",60), collapse=""), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(paste(rep("=",60), collapse=""), "\n")
cat(sprintf("  Input miRNAs:        %d\n", nrow(counts)))
cat(sprintf("  After filtering:     %d\n", nrow(counts_filt)))
cat(sprintf("  Up-regulated:        %d\n", n_up))
cat(sprintf("  Down-regulated:      %d\n", n_down))
cat(sprintf("  Output directory:    %s\n", opt$output))
cat(paste(rep("=",60), collapse=""), "\n\n")
