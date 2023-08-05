#####################

### author honghui Zhao et al

### the sceond step analysis of single cell

#####################

library(Seurat)
library(tidyverse)
library(scater)
library(patchwork)
library(viridis)
options(Seurat.object.assay.version = "v3")

setwd("C:/Users/ZHH/Desktop/program/Human LUAD")

### Load data after processing -------------------
# sce <- readRDS("sce_after_processing.rds")

### View the number of clusters ------------------
head(Idents(sce), 5)
DimPlot(sce,reduction = "umap")

### cell subsets coincide
DimPlot(sce,reduction = "umap",group.by = "RNA_snn_res.0.3")
#ggsave('umap_0.3.pdf', path = "C:/Users/ZHH/Desktop/results", width = 30, height = 16, units = "cm")

DimPlot(sce, reduction = "umap", group.by = "orig.ident", pt.size = .1) 
#ggsave("umap_hamony_orig.ident.png", path = "C:/Users/ZHH/Desktop/results", width = 30, height = 16, units = "cm")

DimPlot(sce,label = T,group.by = "RNA_snn_res.0.3") + NoLegend()
#ggsave("RNA_snn_res.0.9.png", path = "C:/Users/ZHH/Desktop/results/piture", width = 30, height = 16, units = "cm")
sce <- SetIdent(sce,value = "RNA_snn_res.0.3")

###  Main cell type annotation --------------------

mainmarkers <-c("PECAM1", "VWF", "ACTA2",#stromal
                "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", #imm
                "EPCAM", "CDH1", "KRT7", "KRT19"#epi
                )

DotPlot(sce, features =mainmarkers, group.by = "RNA_snn_res.0.5") + 
  coord_flip() + 
  scale_color_viridis()
#ggsave("DotPlot_mainmarkers.png", path = "C:/Users/ZHH/Desktop/results/output", width = 30, height = 8, units = "cm")

Idents(sce) <- sce$RNA_snn_res.0.5

library(readxl)
annotation_curated_main <- read_excel("maincell.xlsx")
new_ids_main <- annotation_curated_main$cellmark
names(new_ids_main) <- levels(sce)
sce <- RenameIdents(sce, new_ids_main)
sce@meta.data$main_cell_type <- Idents(sce)
Idents(sce) <- sce@meta.data$main_cell_type

### figure for all --------------------------

p1 <- DimPlot(sce, group.by = "tissue_type",  pt.size = 0.1)
#ggsave("tissue_type_all.png", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 15, height = 15, units = "cm")

p2 <- DimPlot(sce, group.by = "patient_id", pt.size = 0.1)
#ggsave("patient.id_all.png", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 15, height = 15, units = "cm")

p3 <- DimPlot(sce, group.by = "main_cell_type", pt.size = 0.1)
#ggsave("main_cell_type.png", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 15, height = 15, units = "cm")

p1+p2+p3
ggsave("Fig1.png", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 45, height = 15, units = "cm")

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
#ggsave("cell_cycle.png", path = "C:/Users/ZHH/Desktop/results/output/figure of all", width = 30, height = 15, units = "cm")

### Subset, rescale and save RDS files ------------------------------

epi <- subset(sce, idents = "Epithelial")
imm <- subset(sce, idents = "Immune")
str <- subset(sce, idents = "Stromal")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

saveRDS(epi, file = "epi.RDS")
saveRDS(imm, file = "imm.RDS")
saveRDS(str, file = "str.RDS")

### the sceond step analysis of single cell was end --------------------------