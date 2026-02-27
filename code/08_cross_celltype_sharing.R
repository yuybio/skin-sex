library(dplyr)
library(ggplot2)
library(tidyr)
library(extrafont)
pdf.options(family = "Arial")

setwd("/data/yangyu/Project/skin_sex/DEG")
# NT cell types
###########################################
### 1.  cell types
###########################################
cell_types <- c("whole_skin", "KC", "KC_Channel", "FB", "VEC", "LEC", "Pc", "MEL", "Schwann", "Lymphocyte", "Mac_DC", "LC", "Mast", "HFC")#, "SGC", "MEC")

# Read DE data
read_de_data <- function(cell_type) {
    file_path <- paste0("/data/yangyu/Project/skin_sex/0.all_celltype/", cell_type, "/DE/", cell_type, "_F_vs_M_DE.txt")
    read.table(file_path, sep = '\t', header = TRUE)
}

de_data_list <- setNames(lapply(cell_types, read_de_data), cell_types)

###########################################
### 2. cell type  F/M up
###########################################
# Count F and M genes
count_genes <- function(de_data) {
    F_genes <- sum(de_data$change == "F")
    M_genes <- sum(de_data$change == "M")
    return(c(F_genes, M_genes))
}

results <- data.frame(cell_type = cell_types, t(sapply(de_data_list, count_genes)))
colnames(results)[2:3] <- c("F_genes", "M_genes")

# Print and save results
print(results)
write.table(results, file = "gender_gene_count.txt", sep = "\t", quote = FALSE, row.names = FALSE)

###########################################
### 3.  cell type with whole_skin F/M 
###########################################

whole_skin_DE <- de_data_list[["whole_skin"]]
de_data_excl_ws <- de_data_list[names(de_data_list) != "whole_skin"]

intersection_with_ws <- function(de_data_list, whole_skin_DE, gene_type) {
  ws_genes <- whole_skin_DE$gene[whole_skin_DE$change == gene_type]
  lapply(de_data_list, function(df) intersect(df$gene[df$change == gene_type], ws_genes))
}

F_inter_ws <- intersection_with_ws(de_data_excl_ws, whole_skin_DE, "F")
M_inter_ws <- intersection_with_ws(de_data_excl_ws, whole_skin_DE, "M")

names(F_inter_ws) <- names(de_data_excl_ws)
names(M_inter_ws) <- names(de_data_excl_ws)

F_inter_counts <- sapply(F_inter_ws, length)
M_inter_counts <- sapply(M_inter_ws, length)


F_intersection_with_whole_skin_named <- lapply(seq_along(F_inter_ws), function(i) {
  genes <- F_inter_ws[[i]]
  if (length(genes) > 0) {
    data.frame(gene = genes, cell_type = cell_types[-1][i])
  } else {
    NULL
  }
})
F_intersection_with_whole_skin_named <- Filter(Negate(is.null), F_intersection_with_whole_skin_named)

M_intersection_with_whole_skin_named <- lapply(seq_along(M_inter_ws), function(i) {
  genes <- M_inter_ws[[i]]
  if (length(genes) > 0) {
    data.frame(gene = genes, cell_type = cell_types[-1][i])
  } else {
    NULL
  }
})
M_intersection_with_whole_skin_named <- Filter(Negate(is.null), M_intersection_with_whole_skin_named)


library(openxlsx)
# Save results to Excel with each element in a separate sheet
save_to_excel <- function(intersection_list, file_name) {
  wb <- createWorkbook()
  for (i in seq_along(intersection_list)) {
    addWorksheet(wb, cell_types[-1][i])
    writeData(wb, sheet = cell_types[-1][i], intersection_list[[i]])
  }
  saveWorkbook(wb, file = file_name, overwrite = TRUE)
}

save_to_excel(F_inter_ws, "F_intersection_with_whole_skin.xlsx")
save_to_excel(M_inter_ws, "M_intersection_with_whole_skin.xlsx")


F_inter_counts <- sapply(F_inter_ws, length)
M_inter_counts <- sapply(M_inter_ws, length)

# Print intersection counts
print(F_inter_counts)
print(M_inter_counts)
###########################################
### 4.  stacked barplot（F/M）
###########################################
# Plotting with intersection percentage
plot_data <- function(intersection_counts, gene_counts, gene_type, exclude_cell_types = c()) {
  filtered_cell_types <- setdiff(cell_types[-1], exclude_cell_types)
  specified_order <- c('KC', 'KC_Channel', 'HFC', 'FB', 'Pc', 'VEC', 'LEC', 'MEL', 'Schwann', 'Lymphocyte', 'Mac_DC', 'LC', 'Mast')
  filtered_cell_types <- intersect(specified_order, filtered_cell_types)
  filtered_intersection_counts <- intersection_counts[match(filtered_cell_types, cell_types[-1])]
  filtered_gene_counts <- gene_counts[-1][match(filtered_cell_types, cell_types[-1])]
  
  plot_data <- data.frame(
    cell_type = filtered_cell_types,
    intersection_count = filtered_intersection_counts,
    non_intersection_count = filtered_gene_counts - filtered_intersection_counts
  )
  plot_data_long <- pivot_longer(plot_data, cols = c("intersection_count", "non_intersection_count"), names_to = "intersection_type", values_to = "count")
  plot_data_long$cell_type <- factor(plot_data_long$cell_type, levels = filtered_cell_types)
  
  ggplot(plot_data_long, aes(x = cell_type, y = count, fill = intersection_type)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(x = "Cell Type", y = paste(gene_type, "Gene Count"), fill = "") +
    theme_classic() +
    scale_fill_manual(values = c("intersection_count" = ifelse(gene_type == "M", "blue", "red"), "non_intersection_count" = ifelse(gene_type == "M", "lightblue", "pink"))) +
    geom_text(aes(label = ifelse(intersection_type == "intersection_count", paste0(round((count / (count + plot_data$non_intersection_count[match(cell_type, plot_data$cell_type)])) * 100), "%"), "")), 
              position = position_stack(vjust = 0.5), size = 3)+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_data(F_inter_counts, results$F_genes, "F") 
ggsave('F_intersection_with_whole_skin_counts.pdf',h=4,w=6)
plot_data(M_inter_counts, results$M_genes, "M")
ggsave('M_intersection_with_whole_skin_counts.pdf',h=4,w=6)


###########################################
### 5.  n  cell types F-/M- up
###########################################

get_gene_celltype_pairs <- function(de_data_list, gene_type){
  do.call(rbind, lapply(names(de_data_list), function(ct) {
    genes <- de_data_list[[ct]]$gene[de_data_list[[ct]]$change == gene_type]
    if (length(genes)==0) return(NULL)
    data.frame(gene=genes, cell_type=ct)
  })) %>% distinct()
}

pairs_F <- get_gene_celltype_pairs(de_data_excl_ws, "F")
pairs_M <- get_gene_celltype_pairs(de_data_excl_ws, "M")

count_exact_n <- function(df){
  df %>%
    group_by(gene) %>%
    summarise(n_ct = n_distinct(cell_type), .groups="drop")
}

count_F <- count_exact_n(pairs_F)
count_M <- count_exact_n(pairs_M)


###########################################
### 6. the count of  n  cell types /barplot
###########################################

F_exact <- table(count_F$n_ct)
M_exact <- table(count_M$n_ct)

F_exact <- as.data.frame(F_exact)
M_exact <- as.data.frame(M_exact)
colnames(F_exact) <- c("n_celltypes","gene_count")
colnames(M_exact) <- c("n_celltypes","gene_count")

F_exact$n_celltypes <- as.numeric(as.character(F_exact$n_celltypes))
M_exact$n_celltypes <- as.numeric(as.character(M_exact$n_celltypes))

gg_F <- ggplot(F_exact, aes(x=factor(n_celltypes), y=gene_count)) +
  geom_bar(stat="identity", fill="red") +
  geom_text(aes(label=gene_count), vjust=-0.3, size=3) +
  theme_classic() +
  labs(x="Number of shared cell types", y="F−SBGs Count")

gg_M <- ggplot(M_exact, aes(x=factor(n_celltypes), y=gene_count)) +
  geom_bar(stat="identity", fill="blue") +
  geom_text(aes(label=gene_count), vjust=-0.3, size=3) +
  theme_classic() +
  labs(x="Number of shared cell types", y="M−SBGs Count")

ggsave("F_share_gene_counts.pdf", gg_F, height=3, width=3)
ggsave("M_share_gene_counts.pdf", gg_M, height=3, width=3)

###########################################
### 7. save
###########################################
save_shared_gene_info <- function(pairs, out_file) {
  out <- pairs %>%
    group_by(gene) %>%
    summarise(
      n_ct = n_distinct(cell_type),
      shared_in = paste(sort(unique(cell_type)), collapse=", "),
      .groups="drop"
    ) %>%
    arrange(desc(n_ct))
  write.table(out, out_file, sep="\t", quote=FALSE, row.names=FALSE)
  return(out)
}

shared_F_info <- save_shared_gene_info(pairs_F, "shared_F_exact_genes.txt")
shared_M_info <- save_shared_gene_info(pairs_M, "shared_M_exact_genes.txt")


###########################################
### 8. presence/absence × cell type heatmap
###########################################

library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(patchwork)

# -------------------------------------
# Input:
#   shared_F_gene_info: data.frame(gene, cell_type_count, shared_in_celltypes)
#   shared_M_gene_info: same
#   celltype_order: 
# -------------------------------------


celltype_order <- c("KC", "KC_Channel", "HFC", "FB", "Pc",
                    "VEC", "LEC", "MEL", "Schwann",
                    "Lymphocyte", "Mac_DC", "LC", "Mast")

## ----------------------------
##  cell type × gene 
## ----------------------------
make_binary_mat_from_pairs <- function(pairs, celltype_order) {
  
  pairs2 <- distinct(pairs)  
  

  mat_df <- pairs2 %>%
    mutate(value = 1L) %>%
    tidyr::pivot_wider(names_from = gene,
                       values_from = value,
                       values_fill = list(value = 0L))
  
  ct_ids <- mat_df$cell_type   
  mat <- as.matrix(mat_df[, setdiff(colnames(mat_df), "cell_type") ])
  
  rownames(mat) <- ct_ids
  
  
  mat <- mat[celltype_order, , drop = FALSE]
  
  return(mat)
}


mat_F <- make_binary_mat_from_pairs(pairs_F, celltype_order)
mat_M <- make_binary_mat_from_pairs(pairs_M, celltype_order)

order_cols_by_exact_pattern <- function(mat) {
  
  
  n_ct <- colSums(mat)
  
 
  w <- 2 ^ (seq_len(nrow(mat)) - 1L)
  
 
  pattern_score <- as.vector(t(mat) %*% w)
  

  ord <- order(n_ct, pattern_score)
  
  mat_ord <- mat[, ord, drop = FALSE]
  return(mat_ord)
}

mat_F_ord <- order_cols_by_exact_pattern(mat_F)
mat_M_ord <- order_cols_by_exact_pattern(mat_M)

ht_F <- Heatmap(mat_F_ord,
                name = "F-up",
                col = c("0" = "white", "1" = "#d73027"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,  
                show_column_names = FALSE, 
                use_raster = F, raster_device = "png",
                row_title = "Cell types",border = TRUE,                        
                #rect_gp = gpar(col="black", lwd=1.2),
                column_title = "Shared F-up DEGs")
pdf("F_up_heatmap.pdf", width = 6, height = 3)
draw(ht_F)
dev.off()
ht_M <- Heatmap(mat_M_ord,
                name = "M-up",
                col = c("0" = "white", "1" = "#4575b4"),
                cluster_rows = FALSE,
                cluster_columns = FALSE,
                show_row_names = TRUE,
                show_column_names = FALSE,use_raster = F, raster_device = "png",
                row_title = "Cell types",border = TRUE,
                column_title = "Shared M-up DEGs")
pdf("M_up_heatmap.pdf", width = 6, height =3)
draw(ht_M)
dev.off()
###########################################
### 8.  F/M shared genes（ ≥3 cell types）heatmap
###########################################
# 
get_genes_in_n_cell_types_excluding_whole_skin <- function(de_data_list, gene_type, n) {
  gene_lists <- lapply(de_data_list[-1], function(de_data) de_data$gene[de_data$change == gene_type])
  gene_counts <- table(unlist(gene_lists))
  genes <- names(gene_counts[gene_counts >= n])
  return(genes)
}

# 
genes_F_3 <- get_genes_in_n_cell_types_excluding_whole_skin(de_data_list, "F", 4)
cat(genes_F_3)
df <- bitr(genes_F_3 , fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
head( df )
gene=df[,2]
egoBP <- enrichGO(gene = gene,
                  OrgDb= org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
egoBP <- setReadable(egoBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
go.all=data.frame(egoBP)
write.table(go.all,"F_3_enrich-GO_BP.txt",sep="\t",quote=F,row.names=F)




















