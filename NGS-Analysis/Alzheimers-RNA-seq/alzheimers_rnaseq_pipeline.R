#!/usr/bin/env Rscript
# ============================================================
#  RNA-seq Analysis: Alzheimer's Disease Progression
#  Public Dataset: GSE110226 (NCBI GEO)
# ============================================================
#  Author  : Sana Begum
#  Role    : Research Analyst, Era's Lucknow Medical College
#  Email   : begumsana686@gmail.com
#  GitHub  : https://github.com/Sanabegum09
#
#  Dataset : GSE110226
#    Title : RNA-seq of human hippocampal granule cells in
#            normal aging and Alzheimer's disease
#    Groups: Control (CTL) vs Alzheimer's Disease (AD)
#    Source: NCBI GEO — https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110226
#
#  Multi-stage Alzheimer's Reference:
#    ROSMAP cohort (requires Synapse access: syn3219045)
#    Mayo Clinic RNAseq (temporal cortex + cerebellum)
#    Braak staging I–VI for disease progression analysis
#
#  Pipeline:
#    1. Download count data from GEO (GEOquery)
#    2. Pre-filtering & QC
#    3. DESeq2 normalisation + LFC shrinkage
#    4. PCA, sample QC heatmap
#    5. Differential expression analysis
#    6. Volcano plot (EnhancedVolcano)
#    7. Gene heatmap (pheatmap)
#    8. KEGG & GO pathway enrichment (clusterProfiler)
#    9. AD-specific gene analysis (APOE, APP, PSEN1, MAPT, etc.)
#   10. Export ranked gene list for GSEA
#
#  References:
#    Love MI et al. Genome Biol. 2014;15(12):550.  (DESeq2)
#    Yu G et al. OMICS. 2012;16(5):284-287.        (clusterProfiler)
#    Davis S & Meltzer PS. Bioinformatics. 2007.   (GEOquery)
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery)
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
  library(biomaRt)
  library(viridis)
  library(ggpubr)
  library(optparse)
})

# ── Colours ──────────────────────────────────────────────
COL_AD      <- "#C0392B"   # Red  = AD
COL_CTL     <- "#1A5276"   # Blue = Control
COL_MCI     <- "#E67E22"   # Orange = MCI
COL_UP      <- "#C0392B"
COL_DOWN    <- "#1A5276"
COL_NS      <- "#ABB2B9"

# ── AD Hallmark Genes ─────────────────────────────────────
AD_GENES <- c("APOE","APP","PSEN1","PSEN2","MAPT","TREM2","BIN1",
              "CLU","ABCA7","CR1","PICALM","MS4A6A","MS4A4E",
              "CD33","EPHA1","CD2AP","SORL1","ADAM10","BACE1",
              "SNCA","GBA","PARK2","LRRK2","TNF","IL1B","IL6",
              "HMOX1","SIRT1","BCL2","CASP3","CASP9",
              "SYP","PSD95","GAD1","GFAP","IBA1","AIF1")

# ── Argument parsing ─────────────────────────────────────
option_list <- list(
  make_option("--geo_id",   type="character", default="GSE110226", help="GEO accession"),
  make_option("--output",   type="character", default="results_rnaseq/", help="Output dir"),
  make_option("--padj",     type="numeric",   default=0.05),
  make_option("--lfc",      type="numeric",   default=1.0),
  make_option("--local_counts", type="character", default=NULL,
              help="Path to local count matrix CSV (if GEO download fails)"),
  make_option("--local_meta",   type="character", default=NULL,
              help="Path to local metadata CSV")
)
opt <- parse_args(OptionParser(option_list=option_list))

dir_list <- c(opt$output, file.path(opt$output, c("QC","DE","Plots","Pathways","AD_Genes","GSEA")))
lapply(dir_list, dir.create, recursive=TRUE, showWarnings=FALSE)

cat("\n", paste(rep("=",65), collapse=""), "\n")
cat("  RNA-seq Analysis: Alzheimer's Disease — Sana Begum\n")
cat("  Dataset: ", opt$geo_id, "\n")
cat(paste(rep("=",65), collapse=""), "\n\n")


# ════════════════════════════════════════════════════════════
# 1. DATA ACQUISITION
# ════════════════════════════════════════════════════════════
cat("[ 1/9 ] Acquiring RNA-seq data...\n")

load_data <- function(geo_id, local_counts=NULL, local_meta=NULL) {
  if (!is.null(local_counts) && file.exists(local_counts)) {
    cat("  Loading local count matrix:", local_counts, "\n")
    counts   <- read.csv(local_counts, row.names=1, check.names=FALSE)
    metadata <- read.csv(local_meta, row.names=1, check.names=FALSE)
    return(list(counts=counts, metadata=metadata))
  }

  cat("  Downloading from NCBI GEO:", geo_id, "...\n")
  cat("  (This may take a few minutes on first run)\n")
  tryCatch({
    gse  <- getGEO(geo_id, GSEMatrix=TRUE, getGPL=FALSE)
    gse  <- gse[[1]]
    pdata <- pData(gse)

    # Try to get count data from supplementary files
    supp_files <- getGEOSuppFiles(geo_id, makeDirectory=FALSE)
    count_file <- rownames(supp_files)[grep("count|Count|raw", rownames(supp_files), ignore.case=TRUE)][1]

    if (!is.na(count_file) && file.exists(count_file)) {
      if (grepl("\\.gz$", count_file)) {
        counts <- read.table(gzfile(count_file), header=TRUE, row.names=1, sep="\t", check.names=FALSE)
      } else {
        counts <- read.table(count_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
      }
    } else {
      # Fall back to expression matrix from GEO
      counts <- round(exprs(gse))
    }

    # Build minimal metadata
    condition <- pdata[["characteristics_ch1.1"]]
    if (is.null(condition)) condition <- pdata[["characteristics_ch1"]]
    condition <- gsub(".*: *", "", condition)
    condition <- trimws(condition)
    metadata  <- data.frame(Condition=factor(condition), row.names=colnames(counts))

    return(list(counts=counts, metadata=metadata))
  }, error=function(e) {
    cat("\n  NOTE: GEO download failed (", conditionMessage(e), ")\n")
    cat("  Generating synthetic demonstration data for pipeline testing...\n\n")
    return(generate_demo_data())
  })
}

generate_demo_data <- function(n_genes=15000, n_ctrl=8, n_ad=8) {
  ## Synthetic RNA-seq count data mimicking AD vs Control
  ## Use only for pipeline demonstration — replace with real GEO data
  set.seed(42)
  gene_names <- paste0("GENE_", seq_len(n_genes))
  # Replace first genes with known AD genes for realistic demo
  gene_names[seq_along(AD_GENES)] <- AD_GENES

  ctrl_counts  <- matrix(rnbinom(n_genes*n_ctrl, mu=500, size=20), nrow=n_genes)
  # Introduce DE signal: upregulate AD genes in AD samples
  ad_counts    <- matrix(rnbinom(n_genes*n_ad,   mu=500, size=20), nrow=n_genes)
  ad_up_idx    <- seq_along(AD_GENES)
  ad_counts[ad_up_idx,] <- round(ctrl_counts[ad_up_idx,] * runif(length(ad_up_idx)*n_ad, 2, 5))
  ad_down_idx  <- sample(100:500, 80)
  ad_counts[ad_down_idx,] <- round(ctrl_counts[ad_down_idx,] * runif(length(ad_down_idx)*n_ad, 0.1, 0.4))

  counts   <- as.data.frame(cbind(ctrl_counts, ad_counts))
  colnames(counts) <- c(paste0("CTL_", seq_len(n_ctrl)), paste0("AD_", seq_len(n_ad)))
  rownames(counts) <- gene_names

  metadata <- data.frame(
    Condition = factor(c(rep("Control", n_ctrl), rep("AD", n_ad)), levels=c("Control","AD")),
    row.names = colnames(counts)
  )
  return(list(counts=counts, metadata=metadata))
}

data_obj <- load_data(opt$geo_id, opt$local_counts, opt$local_meta)
counts   <- data_obj$counts
metadata <- data_obj$metadata

# Align
common <- intersect(colnames(counts), rownames(metadata))
counts   <- counts[, common, drop=FALSE]
metadata <- metadata[common, , drop=FALSE]

if (!"Condition" %in% colnames(metadata)) {
  stop("Metadata must have a 'Condition' column. Check your metadata CSV.")
}
metadata$Condition <- relevel(factor(metadata$Condition), ref="Control")
cat(sprintf("  Loaded: %d genes × %d samples\n  Groups: %s\n",
            nrow(counts), ncol(counts),
            paste(levels(metadata$Condition), collapse=" vs ")))


# ════════════════════════════════════════════════════════════
# 2. PRE-FILTERING
# ════════════════════════════════════════════════════════════
cat("[ 2/9 ] Pre-filtering...\n")
keep     <- rowSums(counts >= 10) >= min(3, floor(ncol(counts)*0.2))
counts_f <- counts[keep, ]
cat(sprintf("  Retained %d / %d genes (≥10 reads in ≥20%% samples)\n",
            nrow(counts_f), nrow(counts)))


# ════════════════════════════════════════════════════════════
# 3. DESEQ2 ANALYSIS
# ════════════════════════════════════════════════════════════
cat("[ 3/9 ] Running DESeq2...\n")
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_f)),
  colData   = metadata,
  design    = ~Condition
)
dds <- DESeq(dds, parallel=FALSE, quiet=TRUE)

# VST for visualisation
vst_data <- tryCatch(vst(dds, blind=TRUE), error=function(e) rlog(dds, blind=TRUE))

# LFC shrinkage
coef_name <- resultsNames(dds)[2]
res <- lfcShrink(dds, coef=coef_name, type="apeglm", quiet=TRUE)
res_df <- as.data.frame(res) %>%
  rownames_to_column("Gene") %>%
  arrange(padj, desc(abs(log2FoldChange))) %>%
  mutate(
    DE_Status = case_when(
      padj < opt$padj & log2FoldChange >=  opt$lfc ~ "Up in AD",
      padj < opt$padj & log2FoldChange <= -opt$lfc ~ "Down in AD",
      TRUE ~ "Not significant"
    ),
    neg_log10_padj = -log10(pmax(padj, 1e-300))
  )

n_up   <- sum(res_df$DE_Status == "Up in AD",   na.rm=TRUE)
n_down <- sum(res_df$DE_Status == "Down in AD", na.rm=TRUE)
cat(sprintf("  DE results: %d up, %d down in AD (padj<%.2f, |LFC|>%.1f)\n",
            n_up, n_down, opt$padj, opt$lfc))

write.csv(res_df, file.path(opt$output, "DE", "deseq2_all_genes.csv"), row.names=FALSE)
sig_df <- filter(res_df, DE_Status != "Not significant")
write.csv(sig_df, file.path(opt$output, "DE", "deseq2_significant_genes.csv"), row.names=FALSE)


# ════════════════════════════════════════════════════════════
# 4. QC PLOTS
# ════════════════════════════════════════════════════════════
cat("[ 4/9 ] QC plots...\n")
vst_mat <- assay(vst_data)

## 4a PCA
pca_data <- plotPCA(vst_data, intgroup="Condition", returnData=TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"))
cond_colors <- setNames(c(COL_CTL, COL_AD, COL_MCI),
                        levels(metadata$Condition)[seq_len(nlevels(metadata$Condition))])

p_pca <- ggplot(pca_data, aes(x=PC1, y=PC2, color=group, label=name)) +
  geom_point(size=5, alpha=0.9) +
  geom_text_repel(size=3, max.overlaps=20) +
  stat_ellipse(aes(group=group), type="norm", linetype="dashed", linewidth=0.6) +
  scale_color_manual(values=cond_colors) +
  labs(title="PCA — Alzheimer's Disease RNA-seq (VST)",
       subtitle=paste0("Dataset: ", opt$geo_id, "  |  Author: Sana Begum"),
       x=paste0("PC1 (", pct_var[1], "%)"),
       y=paste0("PC2 (", pct_var[2], "%)"),
       color="Group") +
  theme_bw(base_size=13) +
  theme(plot.title=element_text(face="bold", color=COL_CTL),
        legend.position="right")
ggsave(file.path(opt$output,"QC","pca_rnaseq.pdf"), p_pca, width=8, height=6)
ggsave(file.path(opt$output,"QC","pca_rnaseq.png"), p_pca, width=8, height=6, dpi=300)

## 4b Sample distance heatmap
sdist    <- dist(t(vst_mat))
dist_mat <- as.matrix(sdist)
ann_col  <- data.frame(Condition=metadata$Condition, row.names=rownames(metadata))
ann_col  <- ann_col[colnames(dist_mat), , drop=FALSE]
ann_colors <- list(Condition=cond_colors[levels(metadata$Condition)])

pdf(file.path(opt$output,"QC","sample_distance.pdf"), width=8, height=7)
pheatmap(dist_mat, annotation_col=ann_col, annotation_colors=ann_colors,
         col=colorRampPalette(rev(brewer.pal(9,"Blues")))(100),
         main="Sample-to-Sample Distance (VST)")
dev.off()


# ════════════════════════════════════════════════════════════
# 5. VOLCANO PLOT
# ════════════════════════════════════════════════════════════
cat("[ 5/9 ] Volcano & MA plots...\n")

top_labels <- res_df %>%
  filter(DE_Status != "Not significant") %>%
  arrange(padj) %>% head(30) %>% pull(Gene)
ad_in_res  <- intersect(AD_GENES, res_df$Gene)

p_vol <- EnhancedVolcano(
  res_df,
  lab            = res_df$Gene,
  x              = "log2FoldChange",
  y              = "padj",
  title          = "RNA-seq: Alzheimer's Disease vs Control",
  subtitle       = paste0("Dataset: ", opt$geo_id,
                          "  |  Up: ", n_up, "  |  Down: ", n_down,
                          "\nAuthor: Sana Begum"),
  pCutoff        = opt$padj,
  FCcutoff       = opt$lfc,
  pointSize      = 2.0,
  labSize        = 3.5,
  colAlpha       = 0.7,
  col            = c(COL_NS, COL_NS, COL_DOWN, COL_UP),
  selectLab      = union(top_labels, ad_in_res),
  drawConnectors = TRUE,
  widthConnectors= 0.4,
  colConnectors  = "grey60",
  legendPosition = "right"
)
ggsave(file.path(opt$output,"Plots","volcano_rnaseq.pdf"), p_vol, width=10, height=8)
ggsave(file.path(opt$output,"Plots","volcano_rnaseq.png"), p_vol, width=10, height=8, dpi=300)

## MA plot
p_ma <- ggplot(res_df, aes(x=log10(baseMean+1), y=log2FoldChange, color=DE_Status)) +
  geom_point(alpha=0.5, size=1.2) +
  geom_hline(yintercept=c(-opt$lfc, opt$lfc), linetype="dashed", color="black", linewidth=0.5) +
  scale_color_manual(values=c("Up in AD"=COL_UP, "Down in AD"=COL_DOWN,
                                "Not significant"=COL_NS)) +
  labs(title="MA Plot — AD vs Control",
       x="log10(mean normalised count + 1)", y="log2 Fold Change") +
  theme_bw(base_size=12) +
  theme(plot.title=element_text(face="bold", color=COL_CTL))
ggsave(file.path(opt$output,"Plots","ma_plot.pdf"), p_ma, width=7, height=5)
ggsave(file.path(opt$output,"Plots","ma_plot.png"), p_ma, width=7, height=5, dpi=300)


# ════════════════════════════════════════════════════════════
# 6. TOP DE GENE HEATMAP
# ════════════════════════════════════════════════════════════
cat("[ 6/9 ] Gene heatmap...\n")
top_genes <- sig_df %>% arrange(padj) %>% head(50) %>% pull(Gene)
top_genes <- intersect(top_genes, rownames(vst_mat))
if (length(top_genes) >= 2) {
  hmat   <- vst_mat[top_genes, ]
  hmat_s <- t(scale(t(hmat)))
  ann_c  <- data.frame(Condition=metadata$Condition, row.names=colnames(hmat_s))

  pdf(file.path(opt$output,"Plots","top50_heatmap.pdf"), width=10, height=12)
  pheatmap(hmat_s, annotation_col=ann_c,
           annotation_colors=ann_colors,
           show_colnames=TRUE, show_rownames=TRUE,
           fontsize_row=7, fontsize_col=8,
           color=colorRampPalette(c(COL_DOWN,"white",COL_UP))(100),
           main="Top 50 DE Genes — Alzheimer's vs Control",
           cluster_rows=TRUE, cluster_cols=TRUE,
           border_color=NA)
  dev.off()
}


# ════════════════════════════════════════════════════════════
# 7. AD-SPECIFIC GENE ANALYSIS
# ════════════════════════════════════════════════════════════
cat("[ 7/9 ] AD hallmark gene analysis...\n")
ad_res <- res_df %>% filter(Gene %in% AD_GENES) %>% arrange(padj)
write.csv(ad_res, file.path(opt$output,"AD_Genes","AD_hallmark_genes_results.csv"), row.names=FALSE)

if (nrow(ad_res) > 0) {
  p_ad <- ggplot(ad_res %>% filter(!is.na(log2FoldChange)),
                 aes(x=reorder(Gene, log2FoldChange),
                     y=log2FoldChange,
                     fill=DE_Status)) +
    geom_bar(stat="identity", color="white", width=0.7) +
    geom_hline(yintercept=0, color="black", linewidth=0.4) +
    coord_flip() +
    scale_fill_manual(values=c("Up in AD"=COL_UP, "Down in AD"=COL_DOWN,
                               "Not significant"=COL_NS)) +
    labs(title="AD Hallmark Genes — log2 Fold Change (AD vs Control)",
         subtitle="Author: Sana Begum | Dataset: GSE110226",
         x="Gene", y="log2 Fold Change", fill="DE Status") +
    theme_bw(base_size=12) +
    theme(plot.title=element_text(face="bold", color=COL_CTL),
          legend.position="bottom")
  ggsave(file.path(opt$output,"AD_Genes","AD_genes_LFC.pdf"), p_ad, width=9, height=8)
  ggsave(file.path(opt$output,"AD_Genes","AD_genes_LFC.png"), p_ad, width=9, height=8, dpi=300)

  # AD gene expression heatmap
  ad_in_vst <- intersect(ad_res$Gene, rownames(vst_mat))
  if (length(ad_in_vst) >= 2) {
    ad_hmat <- t(scale(t(vst_mat[ad_in_vst, ])))
    ann_c2  <- data.frame(Condition=metadata$Condition, row.names=colnames(ad_hmat))
    pdf(file.path(opt$output,"AD_Genes","AD_hallmark_heatmap.pdf"), width=9, height=9)
    pheatmap(ad_hmat, annotation_col=ann_c2,
             annotation_colors=ann_colors,
             color=colorRampPalette(c(COL_DOWN,"white",COL_UP))(100),
             main="AD Hallmark Gene Expression",
             fontsize_row=8, border_color=NA,
             cluster_rows=TRUE, cluster_cols=TRUE)
    dev.off()
  }
}


# ════════════════════════════════════════════════════════════
# 8. PATHWAY ENRICHMENT (KEGG + GO)
# ════════════════════════════════════════════════════════════
cat("[ 8/9 ] Pathway enrichment (KEGG + GO)...\n")

## Convert gene symbols to Entrez IDs
convert_to_entrez <- function(gene_symbols) {
  bitr(gene_symbols, fromType="SYMBOL", toType="ENTREZID",
       OrgDb=org.Hs.eg.db, drop=TRUE)$ENTREZID
}

up_entrez   <- tryCatch(convert_to_entrez(filter(sig_df, DE_Status=="Up in AD")$Gene),   error=function(e) character())
down_entrez <- tryCatch(convert_to_entrez(filter(sig_df, DE_Status=="Down in AD")$Gene), error=function(e) character())
all_entrez  <- tryCatch(convert_to_entrez(sig_df$Gene),                                   error=function(e) character())

run_enrichment <- function(gene_ids, direction, out_dir) {
  if (length(gene_ids) < 5) {
    cat(sprintf("  Skipping %s enrichment: too few genes (%d)\n", direction, length(gene_ids)))
    return(invisible(NULL))
  }

  # KEGG
  kegg <- tryCatch(
    enrichKEGG(gene=gene_ids, organism="hsa", pAdjustMethod="BH",
               pvalueCutoff=0.05, qvalueCutoff=0.2),
    error=function(e) NULL)
  if (!is.null(kegg) && nrow(kegg@result) > 0) {
    write.csv(as.data.frame(kegg), file.path(out_dir, paste0("KEGG_",direction,".csv")), row.names=FALSE)
    p <- dotplot(kegg, showCategory=15, font.size=9, title=paste("KEGG —", direction))
    ggsave(file.path(out_dir, paste0("KEGG_",direction,".png")), p, width=9, height=7, dpi=200)
  }

  # GO Biological Process
  go_bp <- tryCatch(
    enrichGO(gene=gene_ids, OrgDb=org.Hs.eg.db, ont="BP",
              pAdjustMethod="BH", pvalueCutoff=0.05, qvalueCutoff=0.2,
              readable=TRUE),
    error=function(e) NULL)
  if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
    write.csv(as.data.frame(go_bp), file.path(out_dir, paste0("GO_BP_",direction,".csv")), row.names=FALSE)
    p2 <- dotplot(go_bp, showCategory=15, font.size=9, title=paste("GO BP —", direction))
    ggsave(file.path(out_dir, paste0("GO_BP_",direction,".png")), p2, width=9, height=7, dpi=200)
  }
}

run_enrichment(up_entrez,   "Up_in_AD",   file.path(opt$output,"Pathways"))
run_enrichment(down_entrez, "Down_in_AD", file.path(opt$output,"Pathways"))


# ════════════════════════════════════════════════════════════
# 9. GSEA EXPORT
# ════════════════════════════════════════════════════════════
cat("[ 9/9 ] Exporting GSEA input & session info...\n")

ranked_list <- res_df %>%
  filter(!is.na(log2FoldChange), !is.na(padj)) %>%
  mutate(rank_metric = sign(log2FoldChange) * (-log10(pmax(padj, 1e-300)))) %>%
  arrange(desc(rank_metric)) %>%
  select(Gene, rank_metric, log2FoldChange, padj)
write.csv(ranked_list, file.path(opt$output,"GSEA","ranked_gene_list.csv"), row.names=FALSE)

# .rnk format for GSEA desktop app
write.table(ranked_list[,c("Gene","rank_metric")],
            file.path(opt$output,"GSEA","ranked_for_gsea.rnk"),
            sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

sink(file.path(opt$output,"session_info.txt"))
cat("RNA-seq Analysis: Alzheimer's Disease\nAuthor: Sana Begum\nDataset: ", opt$geo_id, "\n\n")
print(sessionInfo())
sink()

cat("\n", paste(rep("=",65), collapse=""), "\n")
cat("  ANALYSIS COMPLETE\n")
cat(paste(rep("=",65), collapse=""), "\n")
cat(sprintf("  Dataset    : %s\n", opt$geo_id))
cat(sprintf("  Genes kept : %d / %d after filtering\n", nrow(counts_f), nrow(counts)))
cat(sprintf("  Up in AD   : %d genes\n", n_up))
cat(sprintf("  Down in AD : %d genes\n", n_down))
cat(sprintf("  Outputs    : %s\n", opt$output))
cat(paste(rep("=",65), collapse=""), "\n\n")
