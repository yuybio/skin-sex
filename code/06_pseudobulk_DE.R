# Date: 2025-03-03
# Description: This script uses DaMiRseq to perform normalization, data adjustment, and feature selection on raw pseudobulk count data.
# 
# Input: 
#   - All raw pseudobulk count data
#
# Output: 
#   - Normalized data: vst_normalized_expression.txt
#   - Adjusted data: SV_adjusted_expression.txt
#   - Reduced data: reduced_expression.txt, reduced_expression_rmdup.txt, reduced_expression_rmdup_rank.txt
#   - Selected features

# script for analysis the DE of gender of NT
library(DaMiRseq)
library(limma)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggpubr)
library(data.table)
library(ggsci) 
library(extrafont)
pdf.options(family = "Arial")
library(bioplotr)
library(ggsci)

# Function to read raw count data
read_raw_count <- function(path) {
  read.delim(path, sep = ',')
}


# NT Paths to raw count data
paths <- list(
  whole_skin = "0.all_celltype/whole_skin_pseudobulk_rawCount.csv",
  KC = "0.all_celltype/KC_pseudobulk_rawCount.csv",
  KC_Channel = "0.all_celltype/KC_Channel_pseudobulk_rawCount.csv",
  FB = "0.all_celltype/FB_pseudobulk_rawCount.csv",
  VEC = "0.all_celltype/VEC_pseudobulk_rawCount.csv",
  LEC = "0.all_celltype/LEC_pseudobulk_rawCount.csv",
  Pc = "0.all_celltype/Pc-vSMC_pseudobulk_rawCount.csv",
  MEL = "0.all_celltype/MEL_pseudobulk_rawCount.csv",
  Schwann = "0.all_celltype/Schwann_pseudobulk_rawCount.csv",
  Lymphocyte = "0.all_celltype/Lymphocyte_pseudobulk_rawCount.csv",
  Mac_DC = "0.all_celltype/Mac-DC_pseudobulk_rawCount.csv",
  LC = "0.all_celltype/LC_pseudobulk_rawCount.csv",
  Mast = "0.all_celltype/Mast_pseudobulk_rawCount.csv",
  HFC = "0.all_celltype/HFC_pseudobulk_rawCount.csv",
  SGC = "0.all_celltype/SGC_pseudobulk_rawCount.csv"
)

# Read all raw count data
celltypes <- lapply(paths, read_raw_count)

sampleinfo <- read_tsv("sample_info.txt")

run_damirseq <- function(countdata, sampleinfo, group, celltype, minCounts = 10, fSample = 0.5, th.corr = 0.8, n.pred = 50) {
  if (!dir.exists(celltype)) {dir.create(celltype)}
  setwd(celltype)
  sampleinfo <- sampleinfo %>% filter(sample_group == group)
  countdata <- countdata %>% column_to_rownames("X") %>% as.data.frame()
  countdata <- countdata %>% dplyr::select(intersect(sampleinfo$sample_id, colnames(countdata)))
  keep <- rowSums(countdata) > 100#100
  countdata <- countdata[keep, ]
  
  class1 <- sampleinfo[sampleinfo$sample_id %in% colnames(countdata), ]
  class3 <- class1 %>% column_to_rownames("sample_id") %>% dplyr::rename(class = gender)
  SE <- DaMiR.makeSE(countdata, class3)
  
  data_norm <- DaMiR.normalization(SE, minCounts = minCounts, fSample = fSample, hyper = "no")
  if (!dir.exists('DaMiRseq_filtered')) {dir.create('DaMiRseq_filtered')}
  if (!dir.exists('DaMiRseq_adjust')) {dir.create('DaMiRseq_adjust')}
  write.table(assay(data_norm), file=paste0("DaMiRseq_filtered/", celltype, "_vst_normalized_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)
  data_filt <- DaMiR.sampleFilt(data_norm, th.corr = th.corr)
  sv <- DaMiR.SV(data_filt)
  data_filt$class <- as.factor(data_filt$class)
  data_adjust <- DaMiR.SVadjust(data_filt, sv)
  write.table(assay(data_adjust), file=paste0("DaMiRseq_adjust/", celltype, "_SV_adjusted_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)
  
  plot_and_save <- function(data, prefix) {
    pdf(paste0(prefix, '_all_plot.pdf'), w=10, h=10)
    DaMiR.Allplot(data, colData(data))
    dev.off()
    
    pdf(paste0(prefix, '_umap.pdf'), w=5, h=4)
    plot_umap(assay(data), group=data$class, size=3, label = TRUE, pal_group='d3')
    plot_umap(assay(data), group=data$class, size=3, pal_group='d3')
    dev.off()
    
    pdf(paste0(prefix, '_pca.pdf'), w=5, h=4)
    plot_pca(assay(data), group=data$class, size=3, label = TRUE, pal_group='d3')
    plot_pca(assay(data), group=data$class, size=3, pal_group='d3')
    plot_pca(assay(data), group=data$class, size=3, pcs = c(3L, 4L), label = TRUE, pal_group='d3')
    dev.off()
    
    pdf(paste0(prefix, '_pca_only.pdf'), w=4, h=3)
    plot_pca(assay(data), group=data$class, size=4, pal_group='d3')
    dev.off()
  }
  
  plot_and_save(data_filt, paste0('DaMiRseq_filtered/', celltype))
  plot_and_save(data_adjust, paste0('DaMiRseq_adjust/', celltype))
  
  data_clean <- DaMiR.transpose(assay(data_adjust))
  df <- colData(data_adjust)
  data_reduced <- DaMiR.FSelect(data_clean, df, th.corr = 0.4)
  write.table(data_reduced$data, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)
  
  data_reduced <- DaMiR.FReduct(data_reduced$data)
  write.table(data_reduced, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression_rmdup.txt"), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)
  
  pdf(paste0('DaMiRseq_adjust/', celltype, '_reduced_mds.pdf'), w=10, h=8)
  DaMiR.MDSplot(data_reduced, df)
  dev.off()
  
  tt <- t(data_reduced)
  pdf(paste0('DaMiRseq_adjust/', celltype, '_data_reduced_umap.pdf'), w=5, h=4)
  plot_umap(tt, group=data_adjust$class, size=4, pal_group='d3')
  plot_umap(tt, group=data_adjust$class, size=4, label=TRUE, pal_group='d3')
  dev.off()
  
  pdf(paste0('DaMiRseq_adjust/', celltype, '_data_reduced_pca.pdf'), w=5, h=4)
  plot_pca(tt, group=data_adjust$class, size=3, pal_group='d3')
  plot_pca(tt, group=data_adjust$class, size=3, label=TRUE, pal_group='d3')
  plot_pca(tt, group=data_adjust$class, size=3, pcs=c(3L, 4L), label=TRUE, pal_group='d3')
  dev.off()
  
  df.importance <- DaMiR.FSort(data_reduced, df)
  write.table(df.importance, file=paste0("DaMiRseq_adjust/", celltype, "_reduced_expression_rmdup_rank.txt"), row.names=TRUE, col.names=TRUE, sep='\t', quote=FALSE)
  selected_features4 <- DaMiR.FBest(data_reduced, ranking=df.importance, autoselect='yes')
  selected_features5 <- DaMiR.FBest(data_reduced, ranking=df.importance, n.pred=length(df.importance$RReliefF))
  
  pdf(paste0('DaMiRseq_adjust/', celltype, '_heatmap.pdf'), h=15, w=11)
  DaMiR.Clustplot(selected_features4$data, df)
  dev.off()
  
  pdf(paste0('DaMiRseq_adjust/', celltype, '_pheatmap.pdf'), h=15, w=10)
  pheatmap(t(selected_features5$data),
           clustering_distance_rows="correlation",
           clustering_distance_cols="correlation",
           show_rownames=TRUE,
           scale="row",
           annotation_col=as.data.frame(df) %>% dplyr::select(class),
           annotation_colors=list(class=c('M'="#FED439FF", 'F'="#709AE1FF")))
  dev.off()
}

setwd("0.all_celltype/")
group = 'NT'
original_wd <- getwd()

for (celltype in names(celltypes)) {
  print(celltype)
  run_damirseq(celltypes[[celltype]], sampleinfo, group, celltype, th.corr = 0.7)
  setwd(original_wd)
}


# DE analysis for all adjusted data

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(limma)
library(enrichplot)

one_vs_one_gender <- function(expr, celltype, meta_data, q_value_threshold = 0.05, logFC_threshold = log2(1.5)) {
  celltype_dir <- file.path(getwd(), celltype)
  de_dir <- file.path(celltype_dir, 'DE')
  
  if (!dir.exists(celltype_dir)) dir.create(celltype_dir)
  if (!dir.exists(de_dir)) dir.create(de_dir)
  
  meta_data_group <- meta_data %>% dplyr::filter(sample_id %in% colnames(expr)) 
  group_labels <- as.factor(meta_data_group$gender)
  
  
  group1 <- "F"
  group2 <- "M"
  
  
  subset_idx <- which(group_labels %in% c(group1, group2))
  expr_subset <- expr[, subset_idx]
  labels_subset <- factor(group_labels[subset_idx])
  
 
  design <- model.matrix(~ 0 + labels_subset)
  colnames(design) <- levels(labels_subset)
  print(levels(labels_subset))
  
  # 
  contrast <- makeContrasts(contrasts = paste(group1, "-", group2, sep = ""),
                            levels = design)
  
  # limma
  vfit <- lmFit(expr_subset, design)
  vfit <- contrasts.fit(vfit, contrast)
  efit <- eBayes(vfit)
  summary(decideTests(efit))
  
  # 
  df <- topTable(efit, coef=1, number = Inf, adjust.method = "BH")
  df$gene <- rownames(df)
  df$change = as.factor(ifelse(df$adj.P.Val < q_value_threshold & df$logFC >= logFC_threshold, group1,
                               ifelse(df$adj.P.Val < q_value_threshold & df$logFC <= -logFC_threshold, group2, 'NOT')))
  
  write.table(df, file.path(de_dir, paste0(celltype, '_', group1, '_vs_', group2, "_DE.txt")), row.names = TRUE, col.names = TRUE, sep = '\t', quote = FALSE)
  
  df %>% dplyr::filter(adj.P.Val < 0.05) %>% top_n(n = 20, wt = logFC) -> up_top10
  df %>% dplyr::filter(adj.P.Val < 0.05) %>% top_n(n = 10, wt = -logFC) -> down_top10
  gene <- c(up_top10$gene, down_top10$gene)
  df$sign <- ifelse(df$gene %in% gene, df$gene, NA)
  df <- arrange(df, df$change) # 排序
  
  p <- ggplot(data = df, aes(x = logFC, y = -log10(adj.P.Val), color = change)) + 
    geom_point() + 
    scale_color_manual(values = c("gray", 'blue', 'red'), limits = c("NOT", group2, group1)) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) + 
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 4) +
    theme_classic(base_size = 14) +
    theme(axis.title.x = element_text(size = 14), 
          axis.title.y = element_text(size = 14)) +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) + 
    geom_label_repel(aes(label = sign), 
                     fontface = "bold", 
                     color = "grey50", 
                     box.padding = unit(1, "line"), 
                     point.padding = unit(0.5, "line"), 
                     segment.color = "grey50",
                     max.overlaps = 200) +
    labs(x = "log2FC", y = "-log10(adj.P.Val)", title = "Volcano Plot") 
  
  pdf(file.path(de_dir, paste0(group1, '_', group2, '_volcano.pdf')), h = 5, w = 7)
  print(p)
  dev.off()
  
  # 
  upregulated_genes <- rownames(df[df$adj.P.Val <= q_value_threshold & df$logFC >= logFC_threshold, ])
  downregulated_genes <- rownames(df[df$adj.P.Val <= q_value_threshold & df$logFC <= -logFC_threshold, ])
  
  # 
  return(list(upregulated = upregulated_genes, downregulated = downregulated_genes))
}

# Function to read data and perform DE analysis
perform_DE_analysis <- function(file_path, celltype, meta_data) {
  expr <- read.table(file_path, sep = '\t')
  one_vs_one_gender(expr, celltype, meta_data)
}
# Read meta data
meta_data <- read_tsv("sample_info.txt")

# Define file paths and cell types
file_paths <- list(
  whole_skin = '0.all_celltype/whole_skin/DaMiRseq_adjust/whole_skin_SV_adjusted_expression.txt',
  KC = '0.all_celltype/KC/DaMiRseq_adjust/KC_SV_adjusted_expression.txt',
  KC_Channel = '0.all_celltype/KC_Channel/DaMiRseq_adjust/KC_Channel_SV_adjusted_expression.txt',
  FB = '0.all_celltype/FB/DaMiRseq_adjust/FB_SV_adjusted_expression.txt',
  VEC = '0.all_celltype/VEC/DaMiRseq_adjust/VEC_SV_adjusted_expression.txt',
  LEC = '0.all_celltype/LEC/DaMiRseq_adjust/LEC_SV_adjusted_expression.txt',
  Pc = '0.all_celltype/Pc/DaMiRseq_adjust/Pc_SV_adjusted_expression.txt',
  MEL = '0.all_celltype/MEL/DaMiRseq_adjust/MEL_SV_adjusted_expression.txt',
  Schwann = '0.all_celltype/Schwann/DaMiRseq_adjust/Schwann_SV_adjusted_expression.txt',
  Lymphocyte = '0.all_celltype/Lymphocyte/DaMiRseq_adjust/Lymphocyte_SV_adjusted_expression.txt',
  Mac_DC = '0.all_celltype/Mac_DC/DaMiRseq_adjust/Mac_DC_SV_adjusted_expression.txt',
  LC = '0.all_celltype/LC/DaMiRseq_adjust/LC_SV_adjusted_expression.txt',
  Mast = '0.all_celltype/Mast/DaMiRseq_adjust/Mast_SV_adjusted_expression.txt',
  HFC = '0.all_celltype/HFC/DaMiRseq_adjust/HFC_SV_adjusted_expression.txt',
  SGC = '0.all_celltype/SGC/DaMiRseq_adjust/SGC_SV_adjusted_expression.txt'
)

# Perform DE analysis for each cell type
setwd('0.all_celltype')
original_wd <- getwd()

for (celltype in names(file_paths)) {
  setwd(original_wd)
  perform_DE_analysis(file_paths[[celltype]], celltype, meta_data)
}

