###################################

### Differentially expressed genes analysis of single cell 

### author honghui Zhao et al

###################################
rm(list = ls())

library(Seurat)
library(tidyverse)
library(scater)
library(patchwork)
library(ggrepel)
options(Seurat.object.assay.version = "v3")

setwd("C:/Users/ZHH/Desktop/program/Human LUAD/imm")
sce <- readRDS('imm_anno.rds')
log2FC = 0.5
padj = 0.05 

### DEGs of tumor and normal all -----------
degsall <- FindMarkers(sce, 
                       min.pct = 0, 
                       logfc.threshold = 0,
                       group.by = "tissue_type",
                       ident.1 ="Tumor",
                       ident.2="Normal")

  write.csv(degsall,file = "degsall_imm.csv")

  degsall <- degsall %>%
                     mutate(Difference = pct.1 - pct.2) %>% 
                     rownames_to_column("gene")


  degsall$threshold="ns";
  degsall[which(degsall$avg_log2FC  > log2FC & degsall$p_val_adj <padj),]$threshold="up";
  degsall[which(degsall$avg_log2FC  < (-log2FC) & degsall$p_val_adj < padj),]$threshold="down";
  degsall$threshold=factor(degsall$threshold, levels=c('down','ns','up'))
  table(degsall$threshold)


### genetype 
genename <- read.table("C:/Users/ZHH/Desktop/program/Human LUAD/gene_length_Table.txt",header = T)
lncname <- subset(genename,genetype == 'lincRNA')

## gene of up
up <- subset(degsall, avg_log2FC >= 1  & p_val_adj <= 0.05)
names(up)[1] <- "genename"
uprna <- left_join(up,genename,by = 'genename')
write.csv(uprna,file = 'uprna.csv')
up_down <- subset(degsall,abs(avg_log2FC)>= 0.8 & p_val_adj <= 0.05)
### gene of down
down <- subset(degsall, avg_log2FC <= -1  & p_val_adj <= 0.05)
names(up_down)[1] <- "genename"
up_down <- left_join(up_down,genename,by = 'genename')
write.csv(up_down,file = 'up_down_imm.csv')

### gene of EnhancedVolcano
ggplot(degsall, aes(x=Difference, y=avg_log2FC, color = threshold)) + 
  geom_point(size=0.5) + 
  scale_color_manual(values=c( "blue","grey","red") ) + 
  geom_label_repel(data=subset(degsall, avg_log2FC >= 1 & Difference >= 0.2 & p_val_adj <= 0.05), 
                   aes(label=gene),  #添加label
                   color="black", #设置label中标签的颜色
                   segment.colour = "black",#设置label框的颜色
                   label.padding = 0.1, 
                   #max.overlaps = 200,
                   segment.size = 0.3,  #框的大小
                   size=4)+
  geom_label_repel(data=subset(degsall, avg_log2FC <= -1 & Difference <= -0.2 & p_val_adj <= 0.05), 
                   aes(label=gene), label.padding = 0.1, 
                   color="black",
                   segment.colour = "black",
                   segment.size = 0.3, size=4)+
  geom_vline(xintercept = 0.0,linetype=2)+
  geom_hline(yintercept = 0,linetype=2)+
  theme_classic()
ggsave("火山图.png",width = 50, height = 16, units = "cm")

