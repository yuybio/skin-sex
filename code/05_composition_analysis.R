library(ggpubr)
library(ggplot2)
library(dplyr)
library(data.table)
library(ggsci) 
library(ggsignif)
library(cowplot)
library(extrafont)
pdf.options(family = "Arial")

# celltype ####
meta_data <- read.csv('NT_metadata.csv',row.names = 1)
aa <- meta_data %>% dplyr::select(sample_id,gender) %>% distinct()
dat <- data.frame(table(meta_data$sample_id,meta_data$cell_type))
dat$Total <- apply(dat,1,function(x)sum(dat[dat$Var1 == x[1],3]))
dat <- dat %>% mutate(Percentage = round(Freq/Total,5) * 100)
colnames(dat) <- c('sample_id','cell_type','Freq','Total','Percentage')
dat <- dat %>% inner_join(aa,by=c('sample_id'))
dat$gender <- factor(dat$gender,levels=c('F', 'M'))
dat$cell_type <- as.factor(dat$cell_type)
dat$cell_type <- factor(dat$cell_type,levels = c('KC','KC_Channel','SGC','HFC','Sebocyte',
                                                 'MEC','FB','Pc-vSMC','VEC','LEC',
                                                 'MEL','Schwann',
                                                 'Lymphocyte','Mac-DC','LC','Mast'))
data_wt = reshape2::dcast(dat,sample_id +gender ~ cell_type, value.var = "Percentage")
head(data_wt)


plot_sig=list()
for (i in colnames(data_wt)[3:length(data_wt)]) {
  print(i)
  aa <- data_wt[,c('sample_id','gender',i)]
  colnames(aa) <- c('sample_id','gender','name')
  p <- ggboxplot(aa, x = "gender", y = "name",color = 'gender',add = "jitter") + 
    stat_compare_means(label = "p.format",
                       comparisons = rev(list(c("F","M"))),
                       method = "wilcox",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.001), 
                                          symbols = c("****", "***", "**", "*", "ns")))+
    scale_color_d3()+
    theme(legend.position = "none",panel.grid = element_blank())#+
  
  p <- p + 
    # Add labels and title
    labs(
      x = "",
      y = i
    )
  plot_sig[[i]] <- p
}

pdf('NT_celltype_composion_adj.pdf',h=5,w=13)
plot_grid(plotlist=plot_sig,nrow=2,ncol=7)
dev.off()

# barplot
aa <- meta_data %>% dplyr::select(sample_id,gender) %>% distinct()
dat <- data.frame(table(meta_data$gender,meta_data$cell_type))
dat$Total <- apply(dat,1,function(x)sum(dat[dat$Var1 == x[1],3]))
dat <- dat %>% mutate(Percentage = round(Freq/Total,5) * 100)
colnames(dat) <- c('gender','cell_type','Freq','Total','Percentage')
dat$gender <- factor(dat$gender,levels=c('F', 'M'))
dat$cell_type <- as.factor(dat$cell_type)
dat$cell_type <- factor(dat$cell_type,levels = c('KC','KC_Channel','SGC','HFC','Sebocyte',
                                                 'MEC','FB','Pc-vSMC','VEC','LEC',
                                                 'MEL','Schwann',
                                                 'Lymphocyte','Mac-DC','LC','Mast'))

ggplot(dat, aes(x = gender, y = Percentage, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c(pal_d3("category20")(20),'grey50'))+
  theme_classic()+
  labs(title = "",
       x = "Sample Group",
       y = "Percentage (%)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_flip()
ggsave('NT_celltype_barplot.pdf',h=2,w=8,dpi = 1200)





# cell subtype  ####
aa <- meta_data %>% dplyr::select(sample_id,gender) %>% distinct()
dat <- data.frame(table(meta_data$sample_id,meta_data$celltype_Granular))
dat$Total <- apply(dat,1,function(x)sum(dat[dat$Var1 == x[1],3]))
dat <- dat %>% mutate(Percentage = round(Freq/Total,5) * 100)
colnames(dat) <- c('sample_id','cell_type','Freq','Total','Percentage')
dat <- dat %>% inner_join(aa,by=c('sample_id'))
dat$gender <- factor(dat$gender,levels=c('F', 'M'))
dat$cell_type <- as.factor(dat$cell_type)
data_wt = reshape2::dcast(dat,sample_id +gender ~ cell_type, value.var = "Percentage")
head(data_wt)

plot_sig=list()
for (i in colnames(data_wt)[3:length(data_wt)]) {
  print(i)
  aa <- data_wt[,c('sample_id','gender',i)]
  colnames(aa) <- c('sample_id','gender','name')
  p <- ggboxplot(aa, x = "gender", y = "name",color = 'gender',add = "jitter") + 
    stat_compare_means(label = "p.signif",
                       comparisons = rev(list(c("F","M"))),
                       method = "wilcox",
                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1.001), 
                                          symbols = c("****", "***", "**", "*", "ns")))+
    scale_color_d3()+
    theme(legend.position = "none",panel.grid = element_blank())#+
  
  p <- p + 
    # Add labels and title
    labs(
      x = "",
      y = i
    )
  plot_sig[[i]] <- p
}

pdf('NT_subcelltype_composion_adj.pdf',h=5,w=13)
plot_grid(plotlist=plot_sig,nrow=2,ncol=7)
dev.off()

pdf('NT_subcelltype_composion_adj_all.pdf',h=20,w=13)
plot_grid(plotlist=plot_sig,nrow=10,ncol=7)
dev.off()
# barplot
aa <- meta_data %>% dplyr::select(sample_id,gender) %>% distinct()
dat <- data.frame(table(meta_data$gender,meta_data$celltype_Granular))
dat$Total <- apply(dat,1,function(x)sum(dat[dat$Var1 == x[1],3]))
dat <- dat %>% mutate(Percentage = round(Freq/Total,5) * 100)
colnames(dat) <- c('gender','cell_type','Freq','Total','Percentage')
dat$gender <- factor(dat$gender,levels=c('F', 'M'))
dat$cell_type <- as.factor(dat$cell_type)


