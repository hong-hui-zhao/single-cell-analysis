library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(Seurat)
library(tidyverse)
setwd('~/KRAS/Monocyte/Basis')
sce@meta.data$celltype
anno_col <- c("#00FFFF", "#C1FFC1", "#FF3030", "#8B008B", "#FF7F00", "#FFC125", "#E066FF", "#C6E2FF")#, "#B3EE3A")#, "#FFFF00")#", "#00FF00", "#FF8C00", "#4876FF", "#8B8B00")
names(anno_col) <- c('Mon5', 'Mon2', 'Mon1', 'Mon6', 'Mon4', 'Mon7', 'Mon3', 'Mon8')
markers_genes <- FindAllMarkers(sce, logfc.threshold = 0.2, test.use = "wilcox", assay = "RNA",
                                min.pct = 0.1, min.diff.pct = 0.1, only.pos = TRUE, max.cells.per.ident = 50)

top_genes <- markers_genes %>%
  group_by(cluster) %>%
  top_n(5, wt = avg_log2FC)
top_genes <- top_genes$gene
# get cells mean gene expression
sce[['origCounts']] <- NULL

mean_gene_exp <- AverageExpression(sce,
                                   features = top_genes,
                                   group.by = 'celltype',
                                   slot = 'data') %>%
  data.frame() %>%
  as.matrix()

# add colnames
colnames(mean_gene_exp) <-  c('Mon5', 'Mon2', 'Mon1', 'Mon6', 'Mon4', 'Mon7', 'Mon3', 'Mon8')

# Z-score
htdf <- t(scale(t(mean_gene_exp),scale = T,center = T))

# color
col_fun = colorRamp2(c(-2, 0, 2), c("#0099CC", "white", "#CC0033"))

# top annotation
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = anno_col))
pdf('gene_heatmap.pdf',width = 8,height = 8)
Heatmap(htdf,
        name = "Z-score",
        cluster_columns = F,cluster_rows = F,
        row_names_gp = gpar(fontface = 'italic',fontsize = 8),
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 2),
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = column_ha,
        col = col_fun,
        heatmap_legend_param = list(
          title = "Z-score", at = c(-2,1, 0,-1, -2)))
dev.off()
