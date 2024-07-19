library(Seurat)
library(tidyverse)
setwd('~/KRAS/Basis')
unique(sce$celltype)
celltype <- unique(sce$celltype)

DEGslist <- list()
for (i in 1:length(celltype)) {
  
  sce_data <- subset(sce, idents = celltype[i])
  DEGs <- FindMarkers(sce_data,
                      ident.1 = "G12D", 
                      ident.2 = "Non_G12D",
                      logfc.threshold = 0,
                      group.by="mutation",
                      verbose=T)
  DEGs$celltype <- celltype[i]
  assign(paste0("diff_", celltype[i]),DEGs)
}

diff_cell<-rbind(diff_Monocyte,`diff_T and NK`,diff_MARCO,diff_Epithelial,
                 `diff_B Cell`,diff_Endothlial,diff_Mast, diff_Fibroblast,diff_Plasma)
diff_cell$GENE <- rownames(diff_cell)
write.csv(diff_cell,file = 'G12D vs other mutation.csv')
head(diff_cell)

#显著性
diff_cell$label <- ifelse(diff_cell$p_val_adj<0.01,"adjust P-val<0.01","adjust P-val>=0.01")
for (i in 1:length(celltype)) {
  current_celltype <- celltype[i]
  top <- diff_cell %>%
    filter(celltype == current_celltype) %>%
    distinct(GENE, .keep_all = TRUE) %>%
    top_n(20, abs(avg_log2FC))
  
  assign(paste0("top_", current_celltype), top)
}
top10 <- rbind(top_Monocyte,`top_T and NK`,top_MARCO,top_Epithelial,`top_B Cell`,
               top_Endothlial,top_Mast,top_Fibroblast, top_Plasma)
diff_cell$size <- case_when(!(diff_cell$GENE %in% top10$GENE)~ 1,
                            diff_cell$GENE %in% top10$GENE ~ 2)

#提取非Top10的基因表格；
dt <- filter(diff_cell,size==1)
head(dt)


p <- ggplot()+
  geom_jitter(data = dt,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)

celltype
#根据图p中log2FC区间确定背景柱长度：
dfbar<-data.frame(x=c('Monocyte','T and NK','MARCO','Epithelial','B Cell',
                      'Endothlial','Mast','Fibroblast', 'Plasma'),
                  y=c(2,2,2,2,2,2,2,2,2))
dfbar1<-data.frame(x=c('Monocyte','T and NK','MARCO','Epithelial','B Cell',
                       'Endothlial','Mast','Fibroblast', 'Plasma'),
                   y=c(-2,-2,-2,-2,-2,-2,-2,-2,-2))
#绘制背景柱：
p1 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)



p2 <- ggplot()+
  geom_col(data = dfbar,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = dfbar1,
           mapping = aes(x = x,y = y),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = dt,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 0.85,
              width =0.4)+
  geom_jitter(data = top10,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)


#添加X轴的cluster色块标签：
dfcol<-data.frame(x=c('Monocyte','T and NK','MARCO','Epithelial','B Cell',
                      'Endothlial','Mast','Fibroblast', 'Plasma'),
                  y=0,
                  label=c('Monocyte','T and NK','MARCO','Epithelial','B Cell',
                          'Endothlial','Mast','Fibroblast', 'Plasma'))
mycol <- c("#E64B357F","#00A0877F","#F39B7F7F","#00FFFF", "#C1FFC1", "#E066FF", "#00FA9A", "#8DEEEE", "#FF7F50")
p3 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)
p3


library(ggrepel)
#给每个Cluster差异表达前Top10基因加上标签：
p4 <- p2 + geom_tile(data = dfcol,
                     aes(x=x,y=y),
                     height=0.8,
                     color = "black",
                     fill = mycol,
                     alpha = 0.6,
                     show.legend = F)+
  geom_text_repel(
    data=top10,
    aes(x=celltype,y=avg_log2FC,label=GENE),
    size =3,
    arrow = arrow(length = unit(0.008, "npc"),
                  type = "open", ends = "last"),
    max.overlaps= 100
  )
p4




p5 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","black"))
p5



p6 <- p5+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")
p6



p7 <- p4 +
  scale_color_manual(name=NULL,
                     values = c("red","grey"))+
  labs(x="celltype",y="avg_log2FC")+
  geom_text(data=dfcol,
            aes(x=x,y=y,label=label),
            size =4,
            color ="white")+
  theme_minimal()+
  theme(
    axis.title = element_text(size = 12,
                              color = "black",
                              face = "bold"),
    axis.line.y = element_line(color = "black",
                               size = 0.8),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.direction = "vertical",
    legend.justification = c(1,0),
    legend.text = element_text(size = 12)
  )
p7

ggsvae('Deg.pdf',p7,width = 8,height = 8)

save.image('Degs of G12D vs other mutation')
  