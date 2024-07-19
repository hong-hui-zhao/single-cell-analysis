setwd('/home/zhh/KRAS/TCell')
###################################
path = "/home/zhh/KRAS/TCell"

library(reticulate)
library(cowplot)
library(viridis)
library(magrittr)
library(reshape2)
library(readr)
library(stringr)
library(Seurat)
library(tidyverse)
library(scater)
library(patchwork)
library(readxl)

sce <- readRDS('imm.RDS')
### View the number of clusters ------------------
p1 <- DimPlot(sce,label = T,group.by = "RNA_snn_res.0.3") + NoLegend()
ggsave("RNA_snn_res.0.3.pdf",p1, path = path,width = 8, height = 8)
sce <- SetIdent(sce,value = "RNA_snn_res.0.3")

###,Main cell type annotation --------------------
markers_genes <- FindAllMarkers(sce, logfc.threshold = 0.2, test.use = "wilcox",
                                assay = "RNA",min.pct = 0.1, min.diff.pct = 0.1, 
                                only.pos = TRUE, max.cells.per.ident = 50)
top_genes <- markers_genes %>%
             group_by(cluster) %>%
             top_n(10, wt = avg_log2FC)

write_csv(markers_genes,file = 'markers_genes.csv')

mainmarkers <-c("EPCAM","KRT18","KRT19","CDH1",#epi
                "ACTA2","PDGFRA","VWF","PECAM1",#str
                "JCHAIN","MS4A1","CD68","KIT","PTPRC")

epimarkers <- c('ABCA3', 'SFTPC','SFTPD',"PGC",#AT2
                 'AGER', 'PDPN',#AT1
                 'KRT5','NGFR','SCGB1A1','MUC5B',#Club
                 'TMEM190','FOXJ1',#Ciliated
                 'CHGA','CALCA',#Neuroendocrine'
                 'KRT14', 'TP63', 'DAPL1','CXCL1',#Basal
                 'MUC5AC', 'SPDEF',
                 'PRR4', 'LPO', 'LTF',
                 'CFTR', 'FOXI1', 'ASCL3',
                 'DCLK1', 'ASCL2',
                 'CLTC','PHB','SCGB3A2','CEACAM6'
                 )

immmarkers <- c('CD79A', 'CD24', 'MS4A1', 'CD19',#B
               'CD27', "JCHAIN",'SLAMF7',#Plasma
               'CD3E', 'CD8A', 'GZMK', 'DUSP2',#CD8+ Mem/Eff T Cell
               'CD8', 'GZMH', 'GZMB',#CD8+ Naive T Cell
               'COTL1', 'LDHB',#CD4+ Mem/Eff Cell
               'CD4', 'CCR7', 'LEF1',#CD4+ Naive T Cell
               'KLRD1', 'NKG7', 'TYROBP',#Natural Killer Cell
               'FCER1G',#Natural Killer T Cell
               'TNFRSF4','LTB',#Treg
               'MS4A2', 'CPA3', 'TPSAB1',#Mast
               'MARCO', 'MSR1', 'MRC1',#Macrophage
               'MHCII', 'CD1C', 'PLD4',#Myeloid Dendritic Cell 2
               'CLEC9A', 'LAMP3',#Myeloid Dendritic Cell 1
               'SPP1','APOE','CCL2',"FABP4",
               "CXCL13","RGS1","GAPDH",
               'CTLA4','LDLRAD4','VAV3','IKZF2','TNFRSF18',	
               'IL2R',"HPGD","SFTPB"
)
#sce <- sce[,sce@meta.data$RNA_snn_res.0.5 %in% c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)]

ma <- c('CXCL10','CXCL9','CD68',#M1
        'CD206','CD163','CCL18','MRC1',#M2
        'HLA-DR','TREM2','MARCO', 'MSR1'
)

t <- c(
  'CD3E', 'CD8A', 'GZMK', 'DUSP2',#CD8+ Mem/Eff T Cell
  'CD8', 'GZMH', 'GZMB',#CD8+ Naive T Cell
  'COTL1', 'LDHB',#CD4+ Mem/Eff Cell
  'CD4', 'CCR7', 'LEF1',#CD4+ Naive T Cell
  'KLRD1', 'NKG7', 'TYROBP',#Natural Killer Cell
  'FCER1G',#Natural Killer T Cell
  'TNFRSF4','LTB',#Treg
  'SPP1','APOE','CCL2',"FABP4",
  "CXCL13","RGS1","GAPDH",
  'CTLA4','LDLRAD4','VAV3','IKZF2','TNFRSF18',	
  'IL2R',"HPGD","SFTPB"
)
DotPlot(sce, features =t, group.by = "RNA_snn_res.0.3") + 
        coord_flip() + 
        scale_color_viridis()

ggsave("DotPlot_mainmarkers.pdf", width = 8, height = 8)
sce@meta.data$seurat_clusters

library(readxl)
annotation_curated_main <- read_excel("curated_annotation_main.xlsx")
new_ids_main <- annotation_curated_main$celltype
names(new_ids_main) <- levels(sce)
sce <- RenameIdents(sce, new_ids_main)
sce@meta.data$celltype <- Idents(sce)
Idents(sce) <- sce@meta.data$celltype
sce <-subset(sce,celltype != 'NO')
### figure for all --------------------------
cols<-c('Normal' = "#88c4e8",'Tumor'= '#FD7F76')
p1 <- DimPlot(sce, group.by = "celltype",pt.size = 0.1)
ggsave("CNV_type.pdf",p1, width = 8, height = 5)


p2 <- DimPlot(sce, label = T,group.by = "main_cell_type", pt.size = 0.08,reduction = 'umap') + NoLegend()
ggsave("imm_cell_type.pdf",p2, width = 8, height = 5)


ggsave("Fig2.pdf",p1+p2,width = 12.5, height = 5)

cell_types <- FetchData(sce, vars = c("sample_id", "main_cell_type", "tissue_type")) %>% 
  mutate(main_cell_type = factor(main_cell_type, levels = c("Stromal", "Immune", "Epithelial"))) %>% 
  mutate(sample_id = factor(sample_id, levels = rev(c("p031t","p033t", "p031n", "p033n"))))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = sample_id, fill = main_cell_type, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = c("#456791","#789438","#369854")) +
  coord_flip()
#ggsave("all_barplot.pdf", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 15, height = 30, units = "cm")

### Cell cycle scoring ----------------------------

### add cell cycle, cc.genes loaded with Seurat

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(sce) {
  sce <- CellCycleScoring(sce, s.genes, g2m.genes)
  sce@meta.data$CC.Diff <- sce@meta.data$S.Score - sce@meta.data$G2M.Score
  return(sce)
}

sce <- score_cc(sce)

FeatureScatter(sce, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)
ggsave("cell_cycle.pdf", width = 5, height = 5)

### Subset, rescale and save RDS files ------------------------------

epi <- subset(sce, idents = "Epithelial")
imm <- subset(sce, idents = "T and NK")
MARCO <- subset(sce, idents = "MARCO")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
MARCO <- ScaleData(MARCO)

saveRDS(epi, file = "epi.RDS")
saveRDS(imm, file = "T and NK.RDS")
saveRDS(MARCO, file = "MARCO.RDS")

### the sceond step analysis of single cell was end --------------------------