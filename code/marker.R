library(cowplot)
library(paletteer)
pal <- paletteer_d("ggsci::nrc_npg")[c(1,2,3,4,5,6,7,8,10,11,12,13,14,15,9)]
pdf("epi_celltype.pdf",width=8,height=8)
 DimPlot(sce,label = T ,pt.size = 1,cols = pal)+labs(x = "UMAP1", y = "UMAP2",title = "Celltype") + NoLegend() +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
library(RColorBrewer)
cell_type_cols <- c(brewer.pal(9, "Set1"), 
                    "#FF34B3","#BC8F8F","#20B2AA","#00F5FF","#FFA500","#ADFF2F",
                    "#FF6A6A","#7FFFD4", "#AB82FF","#90EE90","#00CD00","#008B8B",
                    "#6495ED","#FFC1C1","#CD5C5C","#8B008B","#FF3030", "#7CFC00")  

pdf("epi_mutation.pdf",width=8,height=8)
DimPlot(sce, pt.size = 1,group.by = 'mutation')+labs(x = "UMAP1", y = "UMAP2",title = "mutation") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
pdf("epi_CNV.pdf",width=8,height=8)
DimPlot(sce, pt.size = 1,group.by = 'CNV')+labs(x = "UMAP1", y = "UMAP2",title = "CNV") +
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
dev.off()
plot_grid(plot6,plot7,plot8)
