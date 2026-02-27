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

library(clusterProfiler)
library(org.Hs.eg.db)
library(fgsea)

file_paths <- list(
  whole_skin = '0.all_celltype/whole_skin/DE/whole_skin_F_vs_M_DE.txt',
  KC = '0.all_celltype/KC/DE/KC_F_vs_M_DE.txt',
  KC_Channel = '0.all_celltype/KC_Channel/DE/KC_Channel_F_vs_M_DE.txt',
  FB = '0.all_celltype/FB/DE/FB_F_vs_M_DE.txt',
  VEC = '0.all_celltype/VEC/DE/VEC_F_vs_M_DE.txt',
  LEC = '0.all_celltype/LEC/DE/LEC_F_vs_M_DE.txt',
  Pc = '0.all_celltype/Pc/DE/Pc_F_vs_M_DE.txt',
  MEL = '0.all_celltype/MEL/DE/MEL_F_vs_M_DE.txt',
  Schwann = '0.all_celltype/Schwann/DE/Schwann_F_vs_M_DE.txt',
  Lymphocyte = '0.all_celltype/Lymphocyte/DE/Lymphocyte_F_vs_M_DE.txt',
  Mac_DC = '0.all_celltype/Mac_DC/DE/Mac_DC_F_vs_M_DE.txt',
  LC = '0.all_celltype/LC/DE/LC_F_vs_M_DE.txt',
  Mast = '0.all_celltype/Mast/DE/Mast_F_vs_M_DE.txt',
  HFC = '0.all_celltype/HFC/DE/HFC_F_vs_M_DE.txt'
)

read_de <- function(path) {
  read.table(path, sep = '\t')
}
de_files <- lapply(file_paths, read_de)

# get_entrez
get_entrez <- function(df, orgdb = org.Hs.eg.db) {
  gene_info <- bitr(df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb=orgdb)
  dplyr::inner_join(df, gene_info, by = c("gene" = "SYMBOL"))
}

# go
do_go <- function(df, orgdb = org.Hs.eg.db) {
  df$change <- as.factor(df$change)
  cluster_res <- compareCluster(ENTREZID ~ change, data = df, 
                                fun = enrichGO, OrgDb = orgdb, ont = "BP",
                                pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                                pAdjustMethod = "BH", keyType = "ENTREZID")
  setReadable(cluster_res, OrgDb = orgdb, keyType="ENTREZID")
}
plot_go <- function(go_obj, output_file) {
  p <- dotplot(go_obj, showCategory = 20,label_format=100)
  ggsave(output_file, p, width = 10, height = 10)
}

# gsea
gene_set_paths <- list(
  BP = "/home/yangyu/Project/pipeline/database/GSEA_gmt/c5.go.bp.v2023.1.Hs.symbols.gmt",
  HALLMARK = "/home/yangyu/Project/pipeline/database/GSEA_gmt/h.all.v2023.1.Hs.symbols.gmt"
)
do_fgsea <- function(df, gene_set_path) {
  pathways <- gmtPathways(gene_set_path)
  gene_rank <- df$logFC
  names(gene_rank) <- df$gene
  gene_rank <- sort(gene_rank, decreasing = TRUE)
  fgseaRes <- fgsea(pathways = pathways, stats = gene_rank)
  fgseaRes[order(padj), ]
}
plot_fgsea_bar <- function(fgsea_df, output_file, title = "FGSEA Top NES") {
  top_n <- fgsea_df[1:20, ]
  p <- ggplot(top_n, aes(x = NES, y = reorder(pathway, NES))) +
    geom_col(aes(fill = NES > 0)) +
    scale_fill_manual(values = c("dodgerblue","firebrick")) +
    labs(title = title, x = "NES", y = "") +
    theme_minimal()
  ggsave(output_file, p, width = 7, height = 5)
}
output_dir = 'GO/'
names(file_paths)
for (celltype in names(file_paths)) {
  message("Processing ", celltype, "...")
  df <- read_de(file_paths[[celltype]])
  df_entrez <- get_entrez(df)
  
  ## GO
  go_obj <- do_go(df_entrez)
  #write.table(as.data.frame(go_obj), file = file.path(output_dir, paste0(celltype, "_GO.txt")), sep = "\t", row.names = FALSE,quote = F)
  plot_go(go_obj, file.path(output_dir, paste0(celltype, "_GO_dotplot.png")))
  plot_go(go_obj, file.path(output_dir, paste0(celltype, "_GO_dotplot.pdf")))
  
  ## GSEA
  for (gs_name in names(gene_set_paths)) {
    gsea_res <- do_fgsea(df, gene_set_paths[[gs_name]])
    gsea_res_simple <- gsea_res %>% dplyr::select(-leadingEdge)
    write.table(gsea_res_simple, file = file.path(output_dir, paste0(celltype, "_", gs_name, "_GSEA.txt")), sep = "\t", row.names = FALSE,quote = F)
    plot_fgsea_bar(gsea_res, file.path(output_dir, paste0(celltype, "_", gs_name, "_GSEA_bar.png")), title = paste(celltype, gs_name))
    plot_fgsea_bar(gsea_res, file.path(output_dir, paste0(celltype, "_", gs_name, "_GSEA_bar.pdf")), title = paste(celltype, gs_name))
  }
}




#
# whole skin level #####
de_gene <- de_files$whole_skin %>% filter(change!='NOT')%>% dplyr::select(gene,change)
df1 <- bitr(de_gene$gene , fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db )
df2 <- de_gene %>% 
  inner_join(df1,by=c("gene"="SYMBOL"))
head( df2 )
df2$change <- as.factor(df2$change)
table(df2$change)

# GO BP
cgoBP <- compareCluster(ENTREZID~change,
                        data = df2, fun = enrichGO,
                        OrgDb = org.Hs.eg.db,ont='BP',pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,qvalueCutoff = 0.05,keyType = 'ENTREZID')
cgoBP <- setReadable(cgoBP, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
cgo.bp = data.frame(cgoBP)
write.table(cgo.bp,"whole_skin/whole_skin_enrich-GO_BP.txt",sep="\t",quote=F,row.names=F)
table(cgoBP@compareClusterResult$Cluster)
pdf('whole_skin/whole_skin_de_go_BP_dotplot.pdf',h=15,w=10)
print(dotplot(cgoBP,label_format=100,showCategory=40,includeAll=F))
dev.off()

df = cgo.bp %>% filter(ID %in% c('GO:0043270', 'GO:1903522', 'GO:0007188', 'GO:0098742', 'GO:0006816',
                                 'GO:0030198', 'GO:0050900', 'GO:0030595', 'GO:0042060', 'GO:0007599'))
df_plot <- df %>%
  mutate(geneID = gsub("/", ", ", geneID)) %>%
  mutate(geneID_wrapped = str_wrap(geneID, width = 180)) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  mutate(change = change)
p<- ggplot() +
  geom_bar(data = df_plot,
           aes(x = -log10(qvalue), y = Description,fill = change),
           width= 0.5, 
           stat= 'identity') + 
  scale_fill_manual(values=c('F'='#d73027',
                             'M'='#4575b4'))+
  theme_classic() + 
  scale_x_continuous(expand = c(0,0))
p2<- p +
  theme(axis.text.y = element_blank()) + 
  geom_text(data = df_plot,
            aes(x = 0.1, 
                y= Description, 
                label= Description),
            size= 4.5,
            hjust= 0) 
p3<- p2 +
  geom_text(data = df_plot,
            aes(x = 0.1, y = Description, label = geneID_wrapped),
            size= 2,
            fontface= 'italic',
            hjust= 0,
            vjust= 3) 
p4<- p3 +
  labs(x = '-log10(qvalue)', 
       y= 'Term', 
       title= 'whole skin') +
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 14, face = 'bold'),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11),
    axis.ticks.y = element_blank())
p4
ggsave('1.whole_skin_GO.pdf',h=5,w=8)








