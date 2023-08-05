#####################preprocessing

### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)

nFeature_lower <- 500
nFeature_upper <- 10000
nCount_lower <- 1000
nCount_upper <- 100000
pMT_lower <- 0
pMT_upper <- 30
pHB_lower <- 0
pHB_upper <- 5

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


# Data loading and QC
setwd("C:/Users/ZHH/Desktop/公共数据/cellranger~")

### sample list
samples <- read_excel("patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id

### import cellranger files from different data sets
for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0( samples[i])))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i], min.cells = 3))
}

### merge data sets
seu_obj <- merge(seu_obj1, y = c(seu_obj2, seu_obj3, seu_obj4, seu_obj5, seu_obj6, seu_obj7, seu_obj8, seu_obj9, seu_obj10, seu_obj11, seu_obj12, seu_obj13, seu_obj14, seu_obj15, seu_obj16, seu_obj17, seu_obj18, seu_obj19, seu_obj20), add.cell.ids = samples, project = "lung")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^MT-", col.name = "pMT")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^HBA|^HBB", col.name = "pHB")
seu_obj <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL", col.name = "pRP")

qcparams <- c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")
for (i in seq_along(qcparams)){
  print(VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0))
}
for (i in seq_along(qcparams)){
  print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
}

VlnPlot(seu_obj, features = c("nFeature_RNA", "nCount_RNA", "pMT"), pt.size = 0, group.by = "orig.ident", ncol = 1, log = T)
ggsave2("SuppFig1B.pdf", path = "../results", width = 30, height = 20, units = "cm")

### clear environment
remove(seu_obj1)
remove(seu_obj2)
remove(seu_obj3)
remove(seu_obj4)
remove(seu_obj5)
remove(seu_obj6)
remove(seu_obj7)
remove(seu_obj8)
remove(seu_obj9)
remove(seu_obj10)
remove(seu_obj11)
remove(seu_obj12)
remove(seu_obj13)
remove(seu_obj14)
remove(seu_obj15)
remove(seu_obj16)
remove(seu_obj17)
remove(seu_obj18)
remove(seu_obj19)
remove(seu_obj20)

remove(scs_data1)
remove(scs_data2)
remove(scs_data3)
remove(scs_data4)
remove(scs_data5)
remove(scs_data6)
remove(scs_data7)
remove(scs_data8)
remove(scs_data9)
remove(scs_data10)
remove(scs_data11)
remove(scs_data12)
remove(scs_data13)
remove(scs_data14)
remove(scs_data15)
remove(scs_data16)
remove(scs_data17)
remove(scs_data18)
remove(scs_data19)
remove(scs_data20)


# Data Filtering 

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "pMT", "pHB", "pRP")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pMT))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pHB))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = pRP))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pHB, color = nFeature_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), pRP, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pHB, color = nCount_RNA))) + 
      geom_hline(yintercept = pHB_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pHB_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), pRP, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nCount_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(pRP, pMT, color = nFeature_RNA))) + 
      geom_hline(yintercept = pMT_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = pMT_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}


## Before filtering

seu_obj_unfiltered <- seu_obj

#saveRDS(seu_obj_unfiltered, "seurat_obj_unfiltered.RDS")

qc_std_plot(seu_obj_unfiltered)

ggsave2("SuppFig1A.png", path = "../results", width = 30, height = 30, units = "cm")


## After filtering

#seu_obj_unfiltered <- readRDS("all_unfiltered.RDS")

seu_obj <- subset(seu_obj_unfiltered, subset = nFeature_RNA > nFeature_lower & nFeature_RNA < nFeature_upper & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper & pMT < pMT_upper & pHB < pHB_upper)

qc_std_plot(seu_obj)

seu_obj_unfiltered
seu_obj

remove(seu_obj_unfiltered)
# Data normalization

seu_obj <- SCTransform(seu_obj, verbose = T, vars.to.regress = c("nCount_RNA", "pMT"), conserve.memory = T)

#saveRDS(seu_obj, "seurat_objSCTransform.RDS")


# Dimensionality reduction

seu_obj <- RunPCA(seu_obj)
ElbowPlot(seu_obj, ndims = 50)

#for (i in c(15, 20)) {
#  umaptest <- RunUMAP(seu_obj, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "orig.ident") + labs(title = paste0(i, " dimensions")))
#  print(FeaturePlot(umaptest, features = c("EPCAM", "PTPRC"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("MARCO", "KIT"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("FOXJ1", "AGER"), sort.cell = T))
#  print(FeaturePlot(umaptest, features = c("JCHAIN", "VWF"), sort.cell = T))
#  remove(umaptest)
#}

seu_obj <- RunUMAP(seu_obj, dims = 1:15, verbose = T)
seu_obj <- FindNeighbors(seu_obj, dims = 1:15)

for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  seu_obj <- FindClusters(seu_obj, resolution = i)
  print(DimPlot(seu_obj, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

for (i in c("nFeature_RNA", "nCount_RNA", "pMT", "pHB", "pRP")) {
  print(FeaturePlot(seu_obj, features = i, coord.fixed = T, sort.cell = T))
}

#DimPlot(seu_obj, group.by = "orig.ident")
#DimPlot(seu_obj, group.by = "SCT_snn_res.0.5", label = T)


# Main cell type annotation

mainmarkers <- c("PECAM1", "VWF", "ACTA2", "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", "EPCAM", "CDH1", "KRT7", "KRT19")

for (i in seq_along(mainmarkers)) {
  FeaturePlot(seu_obj, features = mainmarkers[i], coord.fixed = T, order = T, cols = viridis(10))
  #ggsave2(paste0("FeaturePlot_mainmarkers_", mainmarkers[i], ".png"), path = "output/annotation", width = 10, height = 10, units = "cm")
}

DotPlot(seu_obj, features = mainmarkers, group.by = "SCT_snn_res.0.2") + 
  coord_flip() + 
  scale_color_viridis()
#ggsave2("DotPlot_mainmarkers.png", path = "output/annotation", width = 30, height = 8, units = "cm")

DimPlot(seu_obj, group.by = "SCT_snn_res.0.2", label = T, label.size = 5)
#ggsave2("DimPlot_all_clusters.png", path = "output/annotation", width = 20, height = 20, units = "cm")

Idents(seu_obj) <- seu_obj$SCT_snn_res.0.2
annotation_curated_main <- read_excel("curated_annotation_main.xlsx")
new_ids_main <- annotation_curated_main$main_cell_type
names(new_ids_main) <- levels(seu_obj)
seu_obj <- RenameIdents(seu_obj, new_ids_main)
seu_obj@meta.data$main_cell_type <- Idents(seu_obj)


# Add metadata

metatable <- read_excel("patients_metadata.xlsx")

metadata <- FetchData(seu_obj, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
rownames(metadata) <- metadata$cell_id

seu_obj <- AddMetaData(seu_obj, metadata = metadata)

# Cell cycle scoring

### add cell cycle, cc.genes loaded with Seurat

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

score_cc <- function(seu_obj) {
  seu_obj <- CellCycleScoring(seu_obj, s.genes, g2m.genes)
  seu_obj@meta.data$CC.Diff <- seu_obj@meta.data$S.Score - seu_obj@meta.data$G2M.Score
  return(seu_obj)
}

seu_obj <- score_cc(seu_obj)

FeatureScatter(seu_obj, "G2M.Score", "S.Score", group.by = "Phase", pt.size = .1) +
  coord_fixed(ratio = 1)


# Subset, rescale and save RDS files

###subset and rescale

#saveRDS(seu_obj, file = "seurat_objall.RDS")

Idents(seu_obj) <- seu_obj@meta.data$main_cell_type

epi <- subset(seu_obj, idents = "Epithelial")
imm <- subset(seu_obj, idents = "Immune")
str <- subset(seu_obj, idents = "Stromal")

epi <- ScaleData(epi)
imm <- ScaleData(imm)
str <- ScaleData(str)

epi
imm
str

#saveRDS(epi, file = "seurat_objectsepi.RDS")
#saveRDS(imm, file = "seurat_objectsimm.RDS")
#saveRDS(str, file = "seurat_objectsstr.RDS")


# Plots for figure 1

###r plots for figure 1

DimPlot(seu_obj, group.by = "tissue_type", cols = use_colors, pt.size = 0.1)
ggsave2("Fig1B.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(seu_obj, group.by = "patient_id", cols = use_colors, pt.size = 0.1)
ggsave2("Fig1C.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(seu_obj, group.by = "main_cell_type", cols = use_colors, pt.size = 0.1)
ggsave2("Fig1D_umap.png", path = "../results", width = 15, height = 15, units = "cm")

cell_types <- FetchData(seu_obj, vars = c("sample_id", "main_cell_type", "tissue_type")) %>% 
  mutate(main_cell_type = factor(main_cell_type, levels = c("Stromal", "Immune", "Epithelial"))) %>% 
  mutate(sample_id = factor(sample_id, levels = rev(c("p018t", "p019t", "p023t", "p024t", "p027t", "p028t", "p030t", "p031t", "p032t", "p033t", "p034t", "p018n", "p019n", "p027n", "p028n", "p029n", "p030n", "p031n", "p032n", "p033n", "p034n"))))

ggplot(data = cell_types) + 
  geom_bar(mapping = aes(x = sample_id, fill = main_cell_type, ), position = "fill", width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("Fig1D_barplot.pdf", path = "../results", width = 15, height = 30, units = "cm")

remove(seu_obj)


##########################################cell type annotation


### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(progeny)
library(scales)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467")


# Load data

#epi <- readRDS("seurat_objects/epi.RDS")
#imm <- readRDS("seurat_objects/imm.RDS")
#str <- readRDS("seurat_objects/str.RDS")


# Rerun PCA, reclustering

### epithelial subclustering

epi <- RunPCA(epi)
ElbowPlot(epi,  ndims = 50)

for (i in c(10, 15, 20, 25)){
  umaptest <- RunUMAP(epi, dims = 1:i, verbose = F)
  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, "dimensions")))
  remove(umaptest)
}

epi <- RunUMAP(epi, dims = 1:20)
epi <- FindNeighbors(epi, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  epi <- FindClusters(epi, resolution = i)
  print(DimPlot(epi, reduction = "umap", label = T) + labs(title = paste0("resolution: ", i)))
}

Idents(epi) <- epi@meta.data$SCT_snn_res.1


### immune subclustering

imm <- RunPCA(imm)
ElbowPlot(imm,  ndims = 50)

#for (i in c(10, 15, 20, 25)){
#  umaptest <- RunUMAP(imm, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  remove(umaptest)
#}

imm <- RunUMAP(imm, dims = 1:20)
imm <- FindNeighbors(imm, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  imm <- FindClusters(imm, resolution = i)
  print(DimPlot(imm, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

Idents(imm) <- imm@meta.data$SCT_snn_res.0.5


### stromal sublustering

str <- RunPCA(str)
ElbowPlot(str, ndims = 50)

#for (i in c(5, 10, 15, 20, 25, 30)){
#  umaptest <- RunUMAP(str, dims = 1:i, verbose = F)
#  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", split.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  print(DimPlot(umaptest, reduction = "umap", group.by = "tissue_type") + labs(title = paste0(i, " dimensions")))
#  remove(umaptest)
#}

str <- RunUMAP(str, dims = 1:20)
str <- FindNeighbors(str, dims = 1:20)
for (i in c(0.2, 0.3, 0.4, 0.5, 1, 2)) {
  str <- FindClusters(str, resolution = i)
  print(DimPlot(str, reduction = "umap") + labs(title = paste0("resolution: ", i)))
}

Idents(str) <- str@meta.data$SCT_snn_res.1



# Define normal and tumor cell clusters

DimPlot(epi, group.by = "SCT_snn_res.1", label = T, repel = T, split.by = "tissue_type")
ggsave2("SuppFig3A.png", path = "../results", width = 30, height = 15, units = "cm")

### compare proportion of cells in a cluster to all epithelial cells for tumor and normal separately, clusters overrepresented in normal samples are supposed to be cell of normal lung parenchyma, all other clusters are supposed to be tumor cells

epi_clusters <- FetchData(epi, vars = c("SCT_snn_res.1", "tissue_type"))

count_tumor <- epi_clusters %>% filter(tissue_type == "Tumor") %>% count() %>% as.numeric()
count_normal <- epi_clusters %>% filter(tissue_type == "Normal") %>% count() %>% as.numeric()

epi_counts <- epi_clusters %>% group_by(tissue_type) %>% count(SCT_snn_res.1)
proportion_tumor <- epi_counts %>% filter(tissue_type == "Tumor") %>% mutate(proportion = n/count_tumor)
proportion_normal <- epi_counts %>% filter(tissue_type == "Normal") %>% mutate(proportion = n/count_normal)

proportion_epi <- full_join(proportion_normal, proportion_tumor, by = "SCT_snn_res.1") %>% 
  mutate(proportion.x = ifelse(is.na(proportion.x), 0,  proportion.x)) %>%  
  mutate(proportion.y = ifelse(is.na(proportion.y), 0,  proportion.y)) %>%
  mutate(tissue_type.x = "Normal") %>%
  mutate(tissue_type.y = "Tumor") %>%
  mutate(cluster_type = ifelse(proportion.x > proportion.y, "Normal", "Tumor"))

cluster_type_data <- left_join(x = epi_clusters, y = proportion_epi, by = "SCT_snn_res.1")
rownames(cluster_type_data) <- rownames(epi_clusters)

epi <- AddMetaData(epi, select(cluster_type_data, cluster_type))


### Bar plot for figure 2

n1 <- select(proportion_epi, c(tissue_type.x, SCT_snn_res.1, proportion.x)) %>%
  mutate(tissue_type = tissue_type.x) %>% 
  mutate(proportion = proportion.x) %>%
  mutate(tissue_type.x = NULL) %>%
  mutate(proportion.x = NULL)
t1 <- select(proportion_epi, c(tissue_type.y, SCT_snn_res.1, proportion.y)) %>%
  mutate(tissue_type = tissue_type.y) %>% 
  mutate(proportion = proportion.y) %>%
  mutate(tissue_type.y = NULL) %>%
  mutate(proportion.y = NULL)

proportion_epi2 <- rbind(n1, t1)

ggplot(proportion_epi2, aes(fill = tissue_type, y = proportion, x = SCT_snn_res.1)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values = use_colors)

ggsave2("SuppFig3B.pdf", path = "../results", width = 40, height = 20, units = "cm")



# Cell type scoring

myeloid_markers <- c("S100A12", "FCN1", "S100A8", "S100A9", "CD14", "CTSS","VCAN", "LYZ", 
                     "MARCO", "FCGR1A", "C1QA", "APOC1", "LGMN", "CTSB", "FCGR3A", 
                     "MAFB", "MAF", "CX3CR1", "ITGAM", "CSF1R",
                     "FABP4", "MCEMP1", 
                     "IL1B", "CXCL8", 
                     "APOE", "CD163", "C1QB", "C1QC", 
                     "FCER1A", "CD1C", "CLEC9A", 
                     "LILRA4", "CLEC4C", "JCHAIN", "IL3RA", "NRP1", 
                     "CLEC10A", "PTCRA", "CCR7", "LAMP3", 
                     "ITGAX", "CD68", "MKI67", "CDK1", "EPCAM")

tcell_nk_markers <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "NKG7", "GNLY", "NCAM1", "KLRD1", "KLRB1", "CD69", "KLRG1", "MKI67", "CDK1", "EPCAM")

bcell_plasma_mast_markers <- c("MS4A1", "CD19", "CD79A", "JCHAIN", "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGKC", "IGLC2", "IGLC3", "CPA3", "KIT", "MS4A2", "GATA2",  "MKI67", "CDK1", "EPCAM")

DotPlot(imm, features = myeloid_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_myeloid_markers.png", path = "output", width = 20, height = 20, units = "cm")

DotPlot(imm, features = tcell_nk_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_T_NK_markers.png", path = "output", width = 20, height = 20, units = "cm")

DotPlot(imm, features = bcell_plasma_mast_markers, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DotPlot_B_Plasma_markers.png", path = "output", width = 20, height = 20, units = "cm")

DimPlot(imm, group.by = "SCT_snn_res.0.5", label = T, split.by = "tissue_type")
ggsave2("DimPlot_immune_clusters.png", path = "output", width = 30, height = 15, units = "cm")

DimPlot(imm, group.by = "patient_id", split.by = "tissue_type", cols = use_colors)
ggsave2("DimPlot_immune_patients.png", path = "output", width = 30, height = 15, units = "cm")

  
  
## Habermann et al.
#https://www.biorxiv.org/content/10.1101/753806v1
  
habermann_epi <- c("ABCA3", "SFTPB", "SFTPC", #AT2
                   "AGER", "PDPN",  #AT1
                   "KRT5",  "NGFR", #Basal
                   "KRT17", "SCGB1A1", "MUC5B", #Goblet
                   "FOXJ1", "TMEM190", "CAPS",#Ciliated
                   "CHGA", "CALCA", "ASCL1" #Neuroendocrine
                   )

habermann_imm <- c("CD3E", "CD4", "FOXP3", "IL7R", "IL2RA", "CD40LG", "CD8A", "CCL5", "NCR1", "KLRB1", "NKG7", "LYZ", "CD68", "ITGAX", "MARCO", "FCGR1A", "FCGR3A", "C1QA", "APOC1", "S100A12", "FCN1", "S100A9", "CD14", "FCER1A", "CD1C", "CD16", "CLEC9A", "LILRA4", "CLEC4C", "JCHAIN", "IGHG1", "IGLL5", "MS4A1", "CD19", "CD79A", "CPA3", "KIT", "MKI67", "CDK1", "EPCAM")

habermann_oth <- c("VWF", "PECAM1", "CCL21", "PROX1", "ACTA2", "MYH11", "PDGFRB", "WT1", "UPK3B", "LUM", "PDGFRA", "MYLK", "HAS1", "PLIN2", "FAP", "PTPRC", "EPCAM")


### Epithelial genes according to Habermann et al.

### manual annotation of normal cell types & detection of immune cell contaminated clusters

DotPlot(subset(epi, subset = cluster_type == "Normal"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2("DimPlot_epi_normal.png", path = "output", width = 30, height = 15, units = "cm")
### manual detection of immune cell contaminated clusters

DotPlot(subset(epi, subset = cluster_type  == "Tumor"), features = habermann_epi, group.by = "SCT_snn_res.1") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave2("DimPlot_epi_tumor.png", path = "output", width = 30, height = 15, units = "cm")


### Immune genes according to Habermann et al.

DotPlot(imm, features = habermann_imm, group.by = "SCT_snn_res.0.5") + coord_flip()
ggsave2("DimPlot_imm_habermann.png", path = "output", width = 30, height = 15, units = "cm")
#for (i in seq_along(habermann_imm)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(imm, features = habermann_imm[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(imm, features = habermann_imm[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist))
#}


### Stromal genes according to Habermann et al.

DotPlot(str, features = habermann_oth, group.by = "SCT_snn_res.1") + 
  coord_flip()
ggsave2("DimPlot_str_habermannoth.png", path = "output", width = 30, height = 15, units = "cm")
#for (i in seq_along(habermann_oth)) {
#  plotlist <- list()
#  plotlist[1] <- FeaturePlot(str, features = habermann_oth[i], sort.cell = T, combine = F)
#  plotlist[2] <- VlnPlot(str, features = habermann_oth[i], pt.size = 0, combine = F)
#  print(CombinePlots(plots = plotlist, ncol = 3))
#}


## Travaglini et al.
#https://www.biorxiv.org/content/10.1101/742320v1

#load marker gene lists

sheets <- paste0("Cluster ", c(1:58))
sheets <- sheets[-43]

signaturelist <- list()

for (i in seq_along(sheets)) {
  a <- read_excel("media-3.xlsx", sheet = sheets[[i]])
  a <- filter(a, a$...2 > 0.7 & a$...4 < 0.3)
  signaturelist <- c(signaturelist,a[1])
  remove(a)
}

#generate list with names of module scores in seurat object

names_of_modulescores <- c()
for (i in seq_along(signaturelist)){
  names_of_modulescores <- c(names_of_modulescores, paste0("T_", names(signaturelist[i]), i))
}

names_of_modulescores <- gsub(names_of_modulescores, pattern = " ", replacement = ".", fixed = TRUE)
names_of_modulescores <- gsub(names_of_modulescores, pattern = "+", replacement = ".", fixed = TRUE)
names_of_modulescores <- gsub(names_of_modulescores, pattern = "/", replacement = ".", fixed = TRUE)


#names_of_modulescores_unfiltered <- c()
#for (i in seq_along(signaturelist)){
#  names_of_modulescores_unfiltered <- c(names_of_modulescores_unfiltered, paste0(names(signaturelist[i]), "_unfiltered", i))
#}
#names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = " ", replacement = ".", fixed = TRUE)
#names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = "+", replacement = ".", fixed = TRUE)
#names_of_modulescores_unfiltered <- gsub(names_of_modulescores_unfiltered, pattern = "/", replacement = ".", fixed = TRUE)
#signature_list_updated <- list()
#for (i in seq_along(sheets)) {
#  signature_list_updated[[i]] <- checkGeneSymbols(signature_list[[i]])
#}

#calculate module scores for different subsets

epi <- AddModuleScore(object = epi, features = signaturelist, name = paste0("T_", names(signaturelist)))
imm <- AddModuleScore(object = imm, features = signaturelist, name = paste0("T_", names(signaturelist)))
str <- AddModuleScore(object = str, features = signaturelist, nbin = 12 , name = paste0("T_", names(signaturelist)))


## Vieira Braga et al.
#https://www.nature.com/articles/s41591-019-0468-5

#load marker gene lists

teichmann_signatures_epi <- read.csv("../data/data/Fig1_DE_Lung_atlas_epithelial.csv")
teichmann_signatures_epi$gene <- as.character(teichmann_signatures_epi$gene)
teichmann_epi <- levels(teichmann_signatures_epi$cluster)

teichmann_signatures_imm <- read.csv("../data/data/Fig2_DE_Lung_atlas_immune.csv")
teichmann_signatures_imm$gene <- as.character(teichmann_signatures_imm$gene)
teichmann_imm <- levels(teichmann_signatures_imm$cluster)

signaturelist2 <- list()

for (i in seq_along(teichmann_epi)) {
  signaturelist2 <- c(signaturelist2, teichmann_signatures_epi %>% filter(cluster == teichmann_epi[i], pct.2 < 0.3, avg_logFC > 0.7) %>% select(gene))
}

for (i in seq_along(teichmann_imm)) {
  signaturelist2 <- c(signaturelist2, teichmann_signatures_imm %>% filter(cluster == teichmann_imm[i], pct.2 < 0.3, avg_logFC > 0.7) %>% select(gene))
}

names(signaturelist2) <- gsub(c(teichmann_epi, teichmann_imm), pattern = "_", replacement = " ")

#generate list with names of module scores in seurat object

names_of_modulescores2 <- c()
for (i in seq_along(signaturelist2)){
  names_of_modulescores2 <- c(names_of_modulescores2, paste0("VB_", names(signaturelist2[i]), i))
}
names_of_modulescores2 <- gsub(names_of_modulescores2, pattern = " ", replacement = ".", fixed = TRUE)

#calculate module scores for different subsets

epi <- AddModuleScore(epi, features = signaturelist2, name = paste0("VB_", names(signaturelist2)))
imm <- AddModuleScore(imm, features = signaturelist2, name = paste0("VB_", names(signaturelist2)))
str <- AddModuleScore(str, features = signaturelist2, nbin = 12, name = paste0("VB_", names(signaturelist2)))


# Curated cell type annotation and subsetting

###epithelial

annotation_curated_epi <- read_excel("curated_annotation_epi.xlsx")
epi_anno <- epi
new_ids_epi <- annotation_curated_epi$cell_type_epi
names(new_ids_epi) <- levels(epi_anno)
epi_anno <- RenameIdents(epi_anno, new_ids_epi)
epi_anno@meta.data$cell_type_epi <- Idents(epi_anno)

epi_anno <- subset(epi_anno, subset = cell_type_epi != "Immune_contamination")
epi_anno <- ScaleData(epi_anno)


###immune

annotation_curated_imm <- read_excel("../data/curated_annotation/curated_annotation_imm.xlsx")
imm_anno <- imm
new_ids_imm <- annotation_curated_imm$cell_type_imm
names(new_ids_imm) <- levels(imm_anno)
imm_anno <- RenameIdents(imm_anno, new_ids_imm)
imm_anno@meta.data$cell_type_imm <- Idents(imm_anno)

imm_anno <- subset(imm_anno, subset = cell_type_imm != "Epithelial_contamination")
imm_anno <- ScaleData(imm_anno)


###stromal

annotation_curated_str <- read_excel("../data/curated_annotation/curated_annotation_str.xlsx")
str_anno <- str
new_ids_str <- annotation_curated_str$cell_type_str
names(new_ids_str) <- levels(str_anno)
str_anno <- RenameIdents(str_anno, new_ids_str)
str_anno@meta.data$cell_type_str <- Idents(str_anno)

str_anno <- subset(str_anno, subset = cell_type_str != "Immune/Epithelial contamination")
str_anno <- ScaleData(str_anno)



###Travaglini et al.

names_of_modulescores_original <- names(signaturelist)


#Epithelial

epi_type <- FetchData(epi_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(epi_type)[i] <- names_of_modulescores_original[i]
}

epi_type %>%
  merge(FetchData(epi_anno, vars = c("cluster_type", "cell_type_epi")), by = 0) %>%
  filter(cluster_type == "Normal") %>%
  group_by(cell_type_epi) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_epi, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_epi, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_epithelial.pdf", path = "../results", width = 30, height = 15, units = "cm")


#Immune

imm_type <- FetchData(imm_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(imm_type)[i] <- names_of_modulescores_original[i]
}

imm_type %>%
  merge(FetchData(imm_anno, vars = "cell_type_imm"), by = 0) %>%
  group_by(cell_type_imm) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_imm, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_imm, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_immune.pdf", path = "../results", width = 30, height = 30, units = "cm")


#Stromal

str_type <- FetchData(str_anno, vars = c(names_of_modulescores))

for(i in seq_along(names_of_modulescores_original)) {
  colnames(str_type)[i] <- names_of_modulescores_original[i]
}

str_type %>%
  merge(FetchData(str_anno, vars = "cell_type_str"), by = 0) %>%
  group_by(cell_type_str) %>%
  pivot_longer(cols = names_of_modulescores_original, names_to = "T_cell_type") %>%
  group_by(cell_type_str, T_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = T_cell_type, y = cell_type_str, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5B_stromal.pdf", path = "../results", width = 30, height = 20, units = "cm")



###Vieira Braga et al.

names_of_modulescores_original2 <- names(signaturelist2)


#Epithelial

epi_type <- FetchData(epi_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(epi_type)[i] <- names_of_modulescores_original2[i]
}

epi_type %>%
  merge(FetchData(epi_anno, vars = c("cluster_type", "cell_type_epi")), by = 0) %>%
  filter(cluster_type == "Normal") %>%
  group_by(cell_type_epi) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_epi, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_epi, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_epithelial.pdf", path = "../results", width = 20, height = 10, units = "cm")


#Immune

imm_type <- FetchData(imm_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(imm_type)[i] <- names_of_modulescores_original2[i]
}

imm_type %>%
  merge(FetchData(imm_anno, vars = "cell_type_imm"), by = 0) %>%
  group_by(cell_type_imm) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_imm, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_imm, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_immune.pdf", path = "../results", width = 20, height = 30, units = "cm")


#Stromal

str_type <- FetchData(str_anno, vars = c(names_of_modulescores2))

for(i in seq_along(names_of_modulescores_original2)) {
  colnames(str_type)[i] <- names_of_modulescores_original2[i]
}

str_type %>%
  merge(FetchData(str_anno, vars = "cell_type_str"), by = 0) %>%
  group_by(cell_type_str) %>%
  pivot_longer(cols = names_of_modulescores_original2, names_to = "VB_cell_type") %>%
  group_by(cell_type_str, VB_cell_type) %>%
  summarise(mean = mean(value)) %>%
  mutate(mean = rescale(mean)) %>%
  ggplot() +
  geom_tile(aes(x = VB_cell_type, y = cell_type_str, fill = mean))+
  scale_fill_gradientn(colors = c("blue", "white", "red"),
                       breaks = c(0, 1),
                       labels = c("0", "1")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(hjust = 1))
ggsave2("SuppFig5A_stromal.pdf", path = "../results", width = 20, height = 20, units = "cm")




#hallmark signatures

broad_pws <- read_lines("../data/data/h.all.v6.2.symbols.gmt") %>%
  lapply(str_split, "//t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

epi_anno <- AddModuleScore(object = epi_anno, features = broad_pws, name = names(broad_pws))
imm_anno <- AddModuleScore(object = imm_anno, features = broad_pws, name = names(broad_pws))
str_anno <- AddModuleScore(object = str_anno, features = broad_pws, name = names(broad_pws), nbin = 12)

#progeny signatures

epi_anno <- progeny(epi_anno, scale = F, organism="Human", top=500, perm=1, return_assay=T)
epi_anno <- ScaleData(epi_anno, assay = "progeny")

### obviously too many cells for progeny function, progeny scores added after immune cells split into myeloid and lymphoid subset
#imm_anno <- progeny(imm_anno, scale = F, organism="Human", top=500, perm=1, return_assay=T)
#imm_anno <- ScaleData(imm_anno, assay = "progeny")

str_anno <- progeny(str_anno, scale = F, organism="Human", top=500, perm=1, return_assay=T)
str_anno <- ScaleData(str_anno, assay = "progeny")

#saveRDS(epi_anno, file = "seurat_objects/epi_anno.RDS")
#saveRDS(imm_anno, file = "seurat_objects/imm_anno.RDS")
#saveRDS(str_anno, file = "seurat_objects/str_anno.RDS")
remove(epi)
remove(imm)
remove(str)

#########################################epithelial analyses


### load libraries
library(ggplot2)
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(tibble)
library(gplots)
library(RColorBrewer)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467",
  AT1 = "#2B8CBE",
  AT2 = "#045A8D",
  Club = "#006D2C",
  DifferentiatingCiliated = "#31A354",
  Ciliated = "#74C476",
  Neuroendocrine = "#8856A7",
  lepidic = "#0B775E",
  acinar = "#74A089",
  `mucinuous (papillary)` = "#E2D200",
  `(micro)papillary` = "#CEAB07",
  solid = "#B40F20",
  sarcomatoid = "#5B1A18",
  CNN = "chartreuse4",
  CNA = "orange")


#load data

#epi_anno <- readRDS("seurat_objects/epi_anno.RDS")

epi_anno@meta.data$cell_type_epi <- factor(epi_anno@meta.data$cell_type_epi, levels = c("AT2",
                                                                                        "AT1",
                                                                                        "Club",
                                                                                        "Ciliated",
                                                                                        "Neuroendocrine",
                                                                                        "Tumor"))

#add inferCNV clone scores (for inferCNV, use seurat object epi_anno and the following code: https://github.com/bischofp/single_cell_lung_adenocarcinoma, computation time ~12h)

scna_scores <- read.delim("../data/inferCNV_output/infercnv_clone_scores_nsclc.tsv") %>% filter(tissue_type == "Tumor") %>% filter(!is.na(cna_clone))
rownames(scna_scores) <- scna_scores$cell_id
scna_scores <- scna_scores %>% select(cna_clone) %>% mutate(cna_clone = as.character(cna_clone))
epi_anno <- AddMetaData(epi_anno, scna_scores)

scna_data <- FetchData(epi_anno, c("tissue_type", "cna_clone"))
scna_data <- scna_data %>% mutate(cna_clone = ifelse(is.na(cna_clone), "CNN", cna_clone))
epi_anno <- AddMetaData(epi_anno, scna_data)



#some plots

DotPlot(subset(epi_anno, subset = cluster_type == "Normal"), features = c("ABCA3", "SFTPC", "AGER", "PDPN",  "KRT5", "TRP63", "NGFR", "SCGB1A1", "MUC5B", "FOXJ1", "TMEM190", "CHGA", "CALCA"), group.by = "cell_type_epi") + 
  coord_flip() + 
  scale_color_viridis()
#ggsave2("DotPlot_markergenes_epi_cell_type.emf", path = "../results", width = 11, height = 8, units = "cm")
ggsave2("Fig2B.png", path = "../results", width = 11, height = 8, units = "cm")


DimPlot(epi_anno, group.by = "cell_type_epi", pt.size = 0.5)
ggsave2("Fig2A_celltype.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
ggsave2("Fig2A_patients.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("Fig2A_tissuetype.png", path = "../results", width = 15, height = 15, units = "cm")

epi_cell_counts <- FetchData(epi_anno, vars = c("tissue_type", "cell_type_epi", "cna_clone", "cluster_type")) %>%
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cell_type_epi)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  scale_y_reverse() +
  coord_flip()
ggsave2("Fig2A_barplot.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cna_clone)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  scale_y_reverse() +
  coord_flip()
ggsave2("SuppFig4_barplot_cna_clone.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = epi_cell_counts, aes(x = tissue_type, fill = cluster_type)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("cyan3", "darkorange2")) +
  coord_flip()
ggsave2("SuppFig4_barplot_cluster_type.pdf", path = "../results", width = 20, height = 5, units = "cm")

DimPlot(epi_anno, group.by = "cna_clone", cols = use_colors, pt.size = 0.5)
ggsave2("SuppFig4_umap_cna_clone.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(epi_anno, group.by = "cluster_type", cols = c("cyan3", "darkorange2"), pt.size = 0.5)
ggsave2("SuppFig4_umap_cluster_type.png", path = "../results", width = 15, height = 15, units = "cm")



###subset tumor epithelial cells

epi_tumor <- SubsetData(SubsetData(epi_anno, subset.name = "cluster_type", accept.value = "Tumor"), subset.name = "tissue_type", accept.value = "Tumor")

epi_tumor <- ScaleData(epi_tumor)


###tumor cell marker genes

Idents(epi_tumor) <- epi_tumor@meta.data$patient_id

markers <- FindAllMarkers(epi_tumor, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_TC_markers <- markers %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

DoHeatmap(epi_tumor, features = top_TC_markers$gene, group.by = "patient_id", draw.lines = F, group.colors = use_colors) +
  scale_fill_viridis()
#ggsave2("HeatMap_Tumor.pdf", path = "output/fig2", width = 30, height = 30, units = "cm")
ggsave2("Fig2C.png", path = "../results", width = 30, height = 30, units = "cm")


###mitotic activity

mitotic_activity <- FetchData(epi_tumor, c("tissue_type", "cell_type_epi", "Phase", "sample_id")) %>%
  mutate(sample_id = factor(sample_id, levels = c("p034t", "p033t", "p032t", "p031t", "p030t", "p027t", "p024t", "p023t", "p019t", "p018t")))

ggplot(mitotic_activity, aes(x = sample_id, fill = Phase)) +
  geom_bar(position = "fill", width = 0.75) +
  scale_fill_manual(values = use_colors) +
  coord_flip()
ggsave2("SuppFig3D.pdf", path = "../results", width = 12, height = 10, units = "cm")



###progeny pathway scores

#clustered heatmap progeny scores

progeny_scores <- as.data.frame(t(GetAssayData(epi_tumor, assay = "progeny", slot = "scale.data")))
progeny_scores$cell_id <- rownames(progeny_scores)
progeny_scores <- gather(progeny_scores, Pathway, Activity, -cell_id)

cells_clusters <- FetchData(epi_tumor, c("sample_id", "cluster_type")) %>% filter(str_detect(sample_id, "t"))
cells_clusters$cell_id <- rownames(cells_clusters)

progeny_scores <- inner_join(progeny_scores, cells_clusters)

summarized_progeny_scores <- progeny_scores %>% 
  group_by(Pathway, sample_id) %>% 
  summarise(avg = mean(Activity), std = sd(Activity)) %>%
  pivot_wider(id_cols = Pathway, names_from = sample_id, values_from = avg) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

pdf("../results/Fig2E.pdf", width = 6, height = 8)
heatmap.2(summarized_progeny_scores, trace = "none", density.info = "none", col = bluered(100))
dev.off()


###correlation mutational status ~ progeny scores

summarized_progeny_scores_mutations <- summarized_progeny_scores %>%
  t() %>%
  as.data.frame()

summarized_progeny_scores_mutations$Sample <- rownames(summarized_progeny_scores_mutations)
summarized_progeny_scores_mutations <- summarized_progeny_scores_mutations %>%
  mutate(KRAS_status = ifelse(Sample %in% c("p018t", "p023t", "p030t", "p031t", "p032t", "p033t"), "mutated", "wildtype")) %>%
  mutate(TP53_status = ifelse(Sample %in% c("p023t", "p027t"), "mutated", "wildtype")) %>%
  mutate(PIK3CA_status = ifelse(Sample %in% c("p031t"), "mutated", "wildtype")) %>%
  mutate(KRAS_status = factor(TP53_status, levels = c("wildtype", "mutated"))) %>%
  mutate(TP53_status = factor(TP53_status, levels = c("wildtype", "mutated")))

ggplot(summarized_progeny_scores_mutations, aes(x = KRAS_status, y = MAPK)) +
  geom_boxplot() +
  ggtitle(paste0(t.test(formula = MAPK~KRAS_status, data = summarized_progeny_scores_mutations, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig3C_KRAS~MAPK.pdf", path = "../results", width = 8, height = 8, units = "cm")

ggplot(summarized_progeny_scores_mutations, aes(x = KRAS_status, y = EGFR)) +
  geom_boxplot() +
  ggtitle(paste0(t.test(formula = EGFR~KRAS_status, data = summarized_progeny_scores_mutations, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig3C_KRAS~EGFR.pdf", path = "../results", width = 8, height = 8, units = "cm")

ggplot(summarized_progeny_scores_mutations, aes(x = TP53_status, y = p53)) +
  geom_boxplot() +
  ggtitle(paste0(t.test(formula = p53~TP53_status, data = summarized_progeny_scores_mutations, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig3C_p53~TP53.pdf", path = "../results", width = 8, height = 8, units = "cm")


### rerun PCA

epi_pca <- epi_tumor

epi_pca <- RunPCA(epi_pca)

DimPlot(epi_pca, reduction = "pca", group.by = "patient_id", dims = c(1,2))
DimPlot(epi_pca, reduction = "pca", group.by = "patient_id", dims = c(3,4))

DimPlot(epi_pca, reduction = "pca", group.by = "histo_subtype")

DimHeatmap(epi_pca, dims = 1, cells = 1000, balanced = T, fast = F, nfeatures = 60) +
  scale_fill_viridis()
#ggsave2("DimHeatmap_epitumor_PC1.pdf", path = "output/fig2", width = 10, height = 20, units = "cm")
ggsave2("SuppFig3F.png", path = "../results", width = 10, height = 20, units = "cm")


### low-dimensional UMAPs

DimPlot(epi_pca, group.by = "patient_id")
DimPlot(epi_pca, group.by = "histo_subtype")

for (i in c(4, 6, 8)){
  umaptest <- RunUMAP(epi_pca, dims = 1:i, verbose = F)
  print(DimPlot(umaptest, reduction = "umap", group.by = "patient_id", cols = use_colors) + labs(title = paste0(i, " dimensions")) + coord_fixed())
  ggsave2(paste0("SuppFig3E_", i, "dimensions_patient.png"), path = "../results", width = 15, height = 10, units = "cm")
  print(DimPlot(umaptest, reduction = "umap", group.by = "histo_subtype", cols = use_colors) + labs(title = paste0(i, " dimensions")) + coord_fixed())
  ggsave2(paste0("SuppFig3E_", i, "dimensions_histo.png"), path = "../results", width = 15, height = 10, units = "cm")
  print(FeaturePlot(umaptest, features = "PC_1", order = T, cols = viridis(10)) + labs(title = paste0(i, " dimensions")) + coord_fixed())
  ggsave2(paste0("SuppFig3E_", i, "dimensions_PC1.png"), path = "../results", width = 15, height = 10, units = "cm")
  remove(umaptest)
}


### bin cells along differentation gradients

# tumor cell signatures "alveolar/club-like" = tumor_signature, "undifferentiated" = anti_tumor_signature

epi_pca <- AddModuleScore(epi_pca, features = list(c("SCGB3A1", "PIGR", "NAPSA", "C4BPA", "SCGB3A2", "HLA-DRA", "CD74", "ADGRF5", "C16orf89", "FOLR1", "SELENBP1", "HLA-DRB1", "ID4", "MGP", "AQP3", "CA2", "LHX9", "HLA-DPB1", "FMO5", "GKN2", "C5", "MUC1", "NPC2", "RNASE1", "PIK3C2G", "SFTA2", "SLC34A2", "HLA-DPA1", "FGFR3", "PGC")), name = "tumor_signature")
epi_pca <- AddModuleScore(epi_pca, features = list(c("DSG2", "CAMK2N1", "FAM3C", "KRT7", "IFI27", "SLC2A1", "MARCKS", "PLAU", "AHNAK2", "PERP", "S100A4", "KRT19", "COL6A1", "UACA", "COL17A1", "CDA", "TPM2", "S100A16", "KRT8", "PRSS23", "DST", "LAMC2", "S100P", "PRSS3", "LAMA3", "DSP", "ITGA3", "MDK", "FAM83A", "ITGB4")), name = "anti_tumor_signature")


epi_tumor_data <- FetchData(epi_pca, c(colnames(epi_pca@meta.data), paste0("PC_", 1:10)))

progeny_scores_data <- as.data.frame(t(GetAssayData(epi_tumor, assay = "progeny", slot = "scale.data")))
progeny_scores_data$cell_id <- rownames(progeny_scores_data)

epi_tumor_data <- full_join(epi_tumor_data, progeny_scores_data, by = "cell_id")

progeny_names <- rownames(epi_pca@assays$progeny)

epi_tumor_data$bin <- cut_number(epi_tumor_data$PC_1, n = 10, labels = c(1:10))

ggplot(epi_tumor_data, mapping = aes(x = bin)) +
  geom_bar()

# patients along PC1

ggplot(epi_tumor_data, mapping = aes(x = bin, fill = factor(patient_id, levels = c("p032", "p018", "p019", "p024", "p031", "p030", "p033", "p023", "p027", "p034")))) +
  geom_bar() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = use_colors)
ggsave2("SuppFig3G.pdf", path = "../results", width = 30, height = 30, units = "cm")

# histological subtypes along PC1

ggplot(epi_tumor_data, mapping = aes(x = bin, fill = factor(histo_subtype, levels = c("lepidic", "acinar", "mucinuous (papillary)", "(micro)papillary", "solid", "sarcomatoid")))) +
  geom_bar() +
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = use_colors)
ggsave2("Fig2F.pdf", path = "../results", width = 30, height = 30, units = "cm")

# cell cycle phase along PC1

ggplot(epi_tumor_data, mapping = aes(x = bin, fill = Phase)) +
  geom_bar() +
  scale_fill_manual(values = use_colors)
ggsave2("SuppFig3H.pdf", path = "../results", width = 30, height = 30, units = "cm")

# normal cell type signatures along PC1

epi_tumor_data %>%
  select(c(bin, VB_Basal.11, VB_Basal.22, VB_Ciliated.13, VB_Ciliated.24, VB_Club5, VB_Goblet.26, VB_Goblet.17, VB_Ionocytes8, VB_Type.2.alveolar9, VB_Type.1.alveolar10)) %>%
  pivot_longer(cols = c(VB_Basal.11, VB_Basal.22, VB_Ciliated.13, VB_Ciliated.24, VB_Club5, VB_Goblet.26, VB_Goblet.17, VB_Ionocytes8, VB_Type.2.alveolar9, VB_Type.1.alveolar10), names_to = "cell_type") %>%
  group_by(bin, cell_type) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ggplot(aes(x = as.numeric(bin), y = mean, color = cell_type)) +
  geom_line(size = 1) +
  scale_color_brewer(type = "qual", palette = "Paired") +
  scale_x_continuous(breaks = c(1:20))
ggsave2("Fig2G.pdf", path = "../results", width = 30, height = 30, units = "cm")

# progeny pathway scores along PC1

getcolors <- colorRampPalette(brewer.pal(14, "Dark2"))

epi_tumor_data %>%
  select(bin, progeny_names) %>%
  pivot_longer(cols = progeny_names, names_to = "progeny") %>%
  group_by(bin, progeny) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ggplot(aes(x = as.numeric(bin), y = mean, color = progeny)) +
  geom_line(size = 1) +
  scale_color_manual(values = getcolors(14)) +
  scale_x_continuous(breaks = c(1:20))
ggsave2("Fig2H.pdf", path = "../results", width = 30, height = 30, units = "cm")

# tumor cell signatures along PC1

epi_tumor_data %>%
  select(bin, tumor_signature1, anti_tumor_signature1) %>%
  pivot_longer(cols = c(tumor_signature1, anti_tumor_signature1), names_to = "genes") %>%
  group_by(bin, genes) %>%
  summarise(mean = mean(value), sd = sd(value)) %>%
  ggplot(aes(x = as.numeric(bin), y = mean, color = genes)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width = 0.3) +
  scale_x_continuous(breaks = c(1:20))
ggsave2("SuppFig3I.pdf", path = "../results", width = 30, height = 30, units = "cm")


#######################################stromal analyses


### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(readr)
library(stringr)
library(gplots)
library(grid)
library(rlang)
library(tibble)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467",
  Endothelial1 = "#FED976",
  Endothelial2 = "#FEB24C",
  Endothelial3 = "#fd8d3C",
  Endothelial4 = "#FC4E2A",
  Endothelial5 = "#E31A1C",
  Endothelial6 = "#BD0026",
  Endothelial7 = "#800026",
  Lymphaticendothelial = "salmon",
  Fibroblast1 = "#2166AC",
  Fibroblast2 = "#4393C3",
  Myofibroblast1 = "#5AAE61",
  Myofibroblast2 = "#1B7837",
  Smoothmuscle1 = "#9970AB",
  Smoothmuscle2 = "#762A83",
  Mesothelial = "#40004B")


#str_anno <- readRDS("seurat_objects/str_anno.RDS")

str_anno@meta.data$cell_type_str <- factor(str_anno@meta.data$cell_type_str, levels = c("Endothelial1",
                                                                                        "Endothelial2",
                                                                                        "Endothelial3",
                                                                                        "Endothelial4",
                                                                                        "Endothelial5",
                                                                                        "Endothelial6",
                                                                                        "Endothelial7",
                                                                                        "Lymphaticendothelial",
                                                                                        "Fibroblast1",
                                                                                        "Fibroblast2",
                                                                                        "Myofibroblast1",
                                                                                        "Myofibroblast2",
                                                                                        "Smoothmuscle1",
                                                                                        "Smoothmuscle2",
                                                                                        "Mesothelial"))



DimPlot(str_anno, group.by = "tissue_type", cols = use_colors)
#ggsave2("DimPlot_str_Normal_Tumor.pdf", path = "output/fig3", width = 15, height = 15, units = "cm")

DimPlot(str_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
ggsave2("SuppFig1C_str_patients.pdf", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(str_anno, group.by = "cell_type_str", label = F, split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
ggsave2("Fig3A_umap.pdf", path = "../results", width = 30, height = 15, units = "cm")

DotPlot(str_anno, features = c("WT1", "UPK3B", "MYH11", "PDGFRB", "ACTA2", "MYLK", "LUM", "PDGFRA", "CCL21", "PROX1", "PECAM1", "VWF"), group.by = "cell_type_str") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
ggsave2("Fig3B.pdf", path = "../results", width = 16, height = 12, units = "cm")


###subsetting

str_endo <- subset(str_anno, subset = cell_type_str %in% c("Endothelial1",
                                                           "Endothelial2",
                                                           "Endothelial3",
                                                           "Endothelial4",
                                                           "Endothelial5",
                                                           "Endothelial6",
                                                           "Endothelial7",
                                                           "Lymphaticendothelial"))
str_endo <- ScaleData(str_endo)

str_fibro <- subset(str_anno, subset = cell_type_str %in% c("Fibroblast1",
                                                            "Fibroblast2",
                                                            "Myofibroblast1",
                                                            "Myofibroblast2",
                                                            "Smoothmuscle1",
                                                            "Smoothmuscle2",
                                                            "Mesothelial"))

str_fibro <- ScaleData(str_fibro)

endo_counts <- FetchData(str_endo, vars = c("tissue_type", "cell_type_str", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

endo_counts_tbl <- endo_counts %>%
  dplyr::count(cell_type_str, patient_id, tissue_type)
write_csv(endo_counts_tbl, path = "../results/SuppTable1.csv")

fibro_counts <- FetchData(str_fibro, vars = c("tissue_type", "cell_type_str", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal"))) 

fibro_counts_tbl <- fibro_counts %>%
  dplyr::count(cell_type_str, patient_id, tissue_type)
write_csv(fibro_counts_tbl, path = "../results/SuppTable2.csv")

ggplot(data = endo_counts, aes(x = tissue_type, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig3A_barplot_endothelial.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = fibro_counts, aes(x = tissue_type, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig3A_barplot_fibroblastic.pdf", path = "../results", width = 20, height = 5, units = "cm")

endo_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig3A_barplot_endothelial_per_patient.pdf", path = "../results", width = 30, height = 30, units = "cm")

fibro_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig3A_barplot_fibroblastic_per_patient.pdf", path = "../results", width = 30, height = 30, units = "cm")

endo_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("SuppFig6A_endothelial.pdf", path = "../results", width = 30, height = 30, units = "cm")

fibro_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_str)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("SuppFig6A_fibroblastic.pdf", path = "../results", width = 30, height = 30, units = "cm")

#heatmap wrapper function
###code from https://github.com/satijalab/seurat/issues/2201

DoMultiBarHeatmap <- function (object, 
                               features = NULL, 
                               cells = NULL, 
                               group.by = "ident", 
                               additional.group.by = NULL, 
                               additional.group.sort.by = NULL, 
                               cols.use = NULL,
                               group.bar = TRUE, 
                               disp.min = -2.5, 
                               disp.max = NULL, 
                               slot = "scale.data", 
                               assay = NULL, 
                               label = TRUE, 
                               size = 5.5, 
                               hjust = 0, 
                               angle = 45, 
                               raster = TRUE, 
                               draw.lines = TRUE, 
                               lines.width = NULL, 
                               group.bar.height = 0.02, 
                               combine = TRUE) 
{
  cells <- cells %||% colnames(x = object)
  if (is.numeric(x = cells)) {
    cells <- colnames(x = object)[cells]
  }
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  features <- features %||% VariableFeatures(object = object)
  ## Why reverse???
  features <- rev(x = unique(x = features))
  disp.max <- disp.max %||% ifelse(test = slot == "scale.data", 
                                   yes = 2.5, no = 6)
  possible.features <- rownames(x = GetAssayData(object = object, 
                                                 slot = slot))
  if (any(!features %in% possible.features)) {
    bad.features <- features[!features %in% possible.features]
    features <- features[features %in% possible.features]
    if (length(x = features) == 0) {
      stop("No requested features found in the ", slot, 
           " slot for the ", assay, " assay.")
    }
    warning("The following features were omitted as they were not found in the ", 
            slot, " slot for the ", assay, " assay: ", paste(bad.features, 
                                                             collapse = ", "))
  }
  
  if (!is.null(additional.group.sort.by)) {
    if (any(!additional.group.sort.by %in% additional.group.by)) {
      bad.sorts <- additional.group.sort.by[!additional.group.sort.by %in% additional.group.by]
      additional.group.sort.by <- additional.group.sort.by[additional.group.sort.by %in% additional.group.by]
      if (length(x = bad.sorts) > 0) {
        warning("The following additional sorts were omitted as they were not a subset of additional.group.by : ", 
                paste(bad.sorts, collapse = ", "))
      }
    }
  }
  
  data <- as.data.frame(x = as.matrix(x = t(x = GetAssayData(object = object, 
                                                             slot = slot)[features, cells, drop = FALSE])))
  
  object <- suppressMessages(expr = StashIdent(object = object, 
                                               save.name = "ident"))
  group.by <- group.by %||% "ident"
  groups.use <- object[[c(group.by, additional.group.by[!additional.group.by %in% group.by])]][cells, , drop = FALSE]
  plots <- list()
  for (i in group.by) {
    data.group <- data
    if (!is_null(additional.group.by)) {
      additional.group.use <- additional.group.by[additional.group.by!=i]  
      if (!is_null(additional.group.sort.by)){
        additional.sort.use = additional.group.sort.by[additional.group.sort.by != i]  
      } else {
        additional.sort.use = NULL
      }
    } else {
      additional.group.use = NULL
      additional.sort.use = NULL
    }
    
    group.use <- groups.use[, c(i, additional.group.use), drop = FALSE]
    
    for(colname in colnames(group.use)){
      if (!is.factor(x = group.use[[colname]])) {
        group.use[[colname]] <- factor(x = group.use[[colname]])
      }  
    }
    
    if (draw.lines) {
      lines.width <- lines.width %||% ceiling(x = nrow(x = data.group) * 
                                                0.0025)
      placeholder.cells <- sapply(X = 1:(length(x = levels(x = group.use[[i]])) * 
                                           lines.width), FUN = function(x) {
                                             return(Seurat:::RandomName(length = 20))
                                           })
      placeholder.groups <- data.frame(rep(x = levels(x = group.use[[i]]), times = lines.width))
      group.levels <- list()
      group.levels[[i]] = levels(x = group.use[[i]])
      for (j in additional.group.use) {
        group.levels[[j]] <- levels(x = group.use[[j]])
        placeholder.groups[[j]] = NA
      }
      
      colnames(placeholder.groups) <- colnames(group.use)
      rownames(placeholder.groups) <- placeholder.cells
      
      group.use <- sapply(group.use, as.vector)
      rownames(x = group.use) <- cells
      
      group.use <- rbind(group.use, placeholder.groups)
      
      for (j in names(group.levels)) {
        group.use[[j]] <- factor(x = group.use[[j]], levels = group.levels[[j]])
      }
      
      na.data.group <- matrix(data = NA, nrow = length(x = placeholder.cells), 
                              ncol = ncol(x = data.group), dimnames = list(placeholder.cells, 
                                                                           colnames(x = data.group)))
      data.group <- rbind(data.group, na.data.group)
    }
    
    order_expr <- paste0('order(', paste(c(i, additional.sort.use), collapse=','), ')')
    group.use = with(group.use, group.use[eval(parse(text=order_expr)), , drop=F])
    
    plot <- Seurat:::SingleRasterMap(data = data.group, raster = raster, 
                                     disp.min = disp.min, disp.max = disp.max, feature.order = features, 
                                     cell.order = rownames(x = group.use), group.by = group.use[[i]])
    
    if (group.bar) {
      pbuild <- ggplot_build(plot = plot)
      group.use2 <- group.use
      cols <- list()
      na.group <- Seurat:::RandomName(length = 20)
      for (colname in rev(x = colnames(group.use2))) {
        if (colname == i) {
          colid = paste0('Identity (', colname, ')')
        } else {
          colid = colname
        }
        
        # Default
        cols[[colname]] <- c(scales::hue_pal()(length(x = levels(x = group.use[[colname]]))))  
        
        #Overwrite if better value is provided
        if (!is_null(cols.use[[colname]])) {
          req_length = length(x = levels(group.use))
          if (length(cols.use[[colname]]) < req_length){
            warning("Cannot use provided colors for ", colname, " since there aren't enough colors.")
          } else {
            if (!is_null(names(cols.use[[colname]]))) {
              if (all(levels(group.use[[colname]]) %in% names(cols.use[[colname]]))) {
                cols[[colname]] <- as.vector(cols.use[[colname]][levels(group.use[[colname]])])
              } else {
                warning("Cannot use provided colors for ", colname, " since all levels (", paste(levels(group.use[[colname]]), collapse=","), ") are not represented.")
              }
            } else {
              cols[[colname]] <- as.vector(cols.use[[colname]])[c(1:length(x = levels(x = group.use[[colname]])))]
            }
          }
        }
        
        # Add white if there's lines
        if (draw.lines) {
          levels(x = group.use2[[colname]]) <- c(levels(x = group.use2[[colname]]), na.group)  
          group.use2[placeholder.cells, colname] <- na.group
          cols[[colname]] <- c(cols[[colname]], "#FFFFFF")
        }
        names(x = cols[[colname]]) <- levels(x = group.use2[[colname]])
        
        y.range <- diff(x = pbuild$layout$panel_params[[1]]$y.range)
        y.pos <- max(pbuild$layout$panel_params[[1]]$y.range) + y.range * 0.015
        y.max <- y.pos + group.bar.height * y.range
        pbuild$layout$panel_params[[1]]$y.range <- c(pbuild$layout$panel_params[[1]]$y.range[1], y.max)
        
        plot <- suppressMessages(plot + 
                                   annotation_raster(raster = t(x = cols[[colname]][group.use2[[colname]]]),  xmin = -Inf, xmax = Inf, ymin = y.pos, ymax = y.max) + 
                                   annotation_custom(grob = grid::textGrob(label = colid, hjust = 0, gp = gpar(cex = 0.75)), ymin = mean(c(y.pos, y.max)), ymax = mean(c(y.pos, y.max)), xmin = Inf, xmax = Inf) +
                                   coord_cartesian(ylim = c(0, y.max), clip = "off")) 
        
        if ((colname == i) && label) {
          x.max <- max(pbuild$layout$panel_params[[1]]$x.range)
          x.divs <- pbuild$layout$panel_params[[1]]$x.major %||% pbuild$layout$panel_params[[1]]$x$break_positions()
          group.use$x <- x.divs
          label.x.pos <- tapply(X = group.use$x, INDEX = group.use[[colname]],
                                FUN = median) * x.max
          label.x.pos <- data.frame(group = names(x = label.x.pos), 
                                    label.x.pos)
          plot <- plot + geom_text(stat = "identity", 
                                   data = label.x.pos, aes_string(label = "group", 
                                                                  x = "label.x.pos"), y = y.max + y.max * 
                                     0.03 * 0.5, angle = angle, hjust = hjust, 
                                   size = size)
          plot <- suppressMessages(plot + coord_cartesian(ylim = c(0, 
                                                                   y.max + y.max * 0.002 * max(nchar(x = levels(x = group.use[[colname]]))) * 
                                                                     size), clip = "off"))
        }
      }
    }
    plot <- plot + theme(line = element_blank())
    plots[[i]] <- plot
  }
  if (combine) {
    plots <- CombinePlots(plots = plots)
  }
  return(plots)
}


###DEGs endothelial

Idents(str_endo) <- str_endo@meta.data$cell_type_str

endo_markers <- FindAllMarkers(str_endo, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_endo_markers <- endo_markers %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

DoMultiBarHeatmap(str_endo, features = top_endo_markers$gene, group.by = "cell_type_str", additional.group.by = "tissue_type",additional.group.sort.by = "tissue_type", cols.use = list(tissue_type = use_colors), draw.lines = F) +
  scale_fill_viridis()
ggsave2("SuppFig6B.png", path = "../results", width = 30, height = 40, units = "cm")
#ggsave2("HeatMap_Endo.pdf", path = "output/fig3", width = 30, height = 40, units = "cm")
#ggsave2("HeatMap_Endo.emf", path = "output/fig3", width = 30, height = 40, units = "cm")


###DEGs fibroblastic

Idents(str_fibro) <- str_fibro@meta.data$cell_type_str

fibro_markers <- FindAllMarkers(str_fibro, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_fibro_markers <- fibro_markers %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

DoMultiBarHeatmap(str_fibro, features = top_fibro_markers$gene, group.by = "cell_type_str", additional.group.by = "tissue_type",additional.group.sort.by = "tissue_type", cols.use = list(tissue_type = use_colors), draw.lines = F) +
  scale_fill_viridis()
ggsave2("Fig6C.png", path = "../results", width = 30, height = 40, units = "cm")
#ggsave2("HeatMap_Fibro.pdf", path = "output/fig3", width = 30, height = 40, units = "cm")
#ggsave2("HeatMap_Fibro.emf", path = "output/fig3", width = 30, height = 40, units = "cm")


###progeny scores

str_fibro2 <- subset(str_anno, subset = cell_type_str %in% c("Fibroblast1",
                                                             "Fibroblast2",
                                                             "Myofibroblast1",
                                                             "Myofibroblast2",
                                                             "Smoothmuscle1",
                                                             "Smoothmuscle2"))

progeny_scores <- as.data.frame(t(GetAssayData(str_fibro2, assay = "progeny", slot = "scale.data")))
progeny_scores$cell_id <- rownames(progeny_scores)
progeny_scores <- gather(progeny_scores, Pathway, Activity, -cell_id)

cells_clusters <- FetchData(str_anno, c("cell_type_str"))
cells_clusters$cell_id <- rownames(cells_clusters)

progeny_scores <- inner_join(progeny_scores, cells_clusters)

summarized_progeny_scores <- progeny_scores %>% 
  group_by(Pathway, cell_type_str) %>% 
  summarise(avg = mean(Activity), std = sd(Activity)) %>%
  pivot_wider(id_cols = Pathway, names_from = cell_type_str, values_from = avg) %>%
  column_to_rownames("Pathway") %>%
  as.matrix()

pdf("../results/Fig3D.pdf", width = 7, height = 10)
heatmap.2(summarized_progeny_scores, trace = "none", density.info = "none", col = bluered(100), margins = c(10,10))
dev.off()


#############################################immune analyses


### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(stringr)
library(cowplot)
library(scales)
library(readr)
library(progeny)
library(gplots)
library(tibble)
library(grid)
library(rlang)

theme_set(theme_cowplot())

use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467",
  Alveolar_Macrophages1 = "#6bAEd6",
  Alveolar_Macrophages2 = "#3182BD",
  Alveolar_Macrophages3 = "#08519C",
  CD14_Macrophages1= "#fff500",
  CD14_Macrophages2= "#FE9929",
  CD14_Macrophages3= "#EC7014",
  CD14_Macrophages4= "#CC4C02",
  CD14_Macrophages5= "#8C2D04",
  Macrophages_Proliferating= "#E31A1C",
  Monocytes= "#FA9FB5",
  Myeloid_Dendritic= "#DD3497",
  Plasmacytoid_Dendritic= "#7A0177",
  T_conv1= "#c2e699",
  T_conv2= "#78c679",
  T_reg= "#006837",
  T_CD8_1= "#bcbddc",
  T_CD8_2= "#9e9ac8",
  T_CD8_3= "#807dba",
  T_CD8_Proliferating= "#6a51a3",
  NK_cells= "#4a1486",
  B_cells= "#969696",
  Plasma= "#636363",
  Mast= "#252525")


#load data

#imm_anno <- readRDS("seurat_objects/imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages1",
                                                                                         "Alveolar_Macrophages2",
                                                                                         "Alveolar_Macrophages3",
                                                                                         "CD14_Macrophages1",
                                                                                         "CD14_Macrophages2",
                                                                                         "CD14_Macrophages3",
                                                                                         "CD14_Macrophages4",
                                                                                         "CD14_Macrophages5",
                                                                                         "Macrophages_Proliferating",
                                                                                         "Monocytes",
                                                                                         "Myeloid_Dendritic",
                                                                                         "Plasmacytoid_Dendritic",
                                                                                         "Mast",
                                                                                         "T_conv1",
                                                                                         "T_conv2",
                                                                                         "T_reg",
                                                                                         "T_CD8_1",
                                                                                         "T_CD8_2",
                                                                                         "T_CD8_3",
                                                                                         "T_CD8_Proliferating",
                                                                                         "NK_cells",
                                                                                         "B_cells",
                                                                                         "Plasma"))


DimPlot(imm_anno, group.by = "tissue_type", cols = use_colors)
#ggsave2("DimPlot_imm_Normal_Tumor.pdf", path = "output/fig4", width = 15, height = 15, units = "cm")
#ggsave2("DimPlot_imm_Normal_Tumor.png", path = "output/fig4", width = 35, height = 15, units = "cm")

DimPlot(imm_anno, group.by = "patient_id", cols = use_colors, pt.size = 0.5)
#ggsave2("DimPlot_imm_patients.pdf", path = "output/fig4", width = 30, height = 15, units = "cm")
ggsave2("SuppFig1C_imm_patients.png", path = "../results", width = 15, height = 15, units = "cm")

DimPlot(imm_anno, group.by = "cell_type_imm", split.by = "tissue_type", cols = use_colors, pt.size = 0.5)
#ggsave2("DimPlot_imm_celltype.pdf", path = "output/fig4", width = 35, height = 15, units = "cm")
ggsave2("Fig4A_umap.png", path = "../results", width = 35, height = 15, units = "cm")

DotPlot(imm_anno, features = c("CD68", "LYZ", "FABP4", "MARCO", "LGMN", "CSF1R", "CD14", "S100A12", "FCN1", "CD1C", "FCER1A", "LILRA4", "IL3RA", "KIT", "GATA2", "CD3E", "CD4", "FOXP3", "IL2RA", "CD8A", "NKG7", "KLRD1", "MS4A1", "CD79A", "JCHAIN", "IGKC", "MKI67"), group.by = "cell_type_imm") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  scale_color_viridis()
ggsave2("Fig4B.pdf", path = "../results", width = 20, height = 20, units = "cm")


###subsetting

imm_lympho <- subset(imm_anno, subset = cell_type_imm %in% c("T_conv1",
                                                             "T_conv2",
                                                             "T_reg",
                                                             "T_CD8_1",
                                                             "T_CD8_2",
                                                             "T_CD8_3",
                                                             "T_CD8_Proliferating",
                                                             "NK_cells",
                                                             "B_cells",
                                                             "Plasma"))

imm_lympho <- ScaleData(imm_lympho)

imm_myelo <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages1",
                                                            "Alveolar_Macrophages2",
                                                            "Alveolar_Macrophages3",
                                                            "CD14_Macrophages1",
                                                            "CD14_Macrophages2",
                                                            "CD14_Macrophages3",
                                                            "CD14_Macrophages4",
                                                            "CD14_Macrophages5",
                                                            "Macrophages_Proliferating",
                                                            "Monocytes",
                                                            "Myeloid_Dendritic",
                                                            "Plasmacytoid_Dendritic",
                                                            "Mast"))

imm_myelo <- ScaleData(imm_myelo)

lympho_counts <- FetchData(imm_lympho, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

lympho_counts_tbl <- lympho_counts %>%
  dplyr::count(cell_type_imm, patient_id, tissue_type)
write_csv(lympho_counts_tbl, path = "../results/SuppTable4.csv")

myelo_counts <- FetchData(imm_myelo, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal"))) 

myelo_counts_tbl <- myelo_counts %>%
  dplyr::count(cell_type_imm, patient_id, tissue_type)
write_csv(myelo_counts_tbl, path = "../results/SuppTable3.csv")


ggplot(data = lympho_counts, aes(x = tissue_type, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig4A_barplot_lymphoid.pdf", path = "../results", width = 20, height = 5, units = "cm")

ggplot(data = myelo_counts, aes(x = tissue_type, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig4A_barplot_myeloid.pdf", path = "../results", width = 20, height = 5, units = "cm")

lympho_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig4A_barplot_lymphoid_per_patient.pdf", path = "../results", width = 30, height = 30, units = "cm")

myelo_counts %>%
  filter(tissue_type == "Tumor") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("Fig4A_barplot_myeloid_per_patient.pdf", path = "../results", width = 30, height = 30, units = "cm")

lympho_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("SuppFig7A_lymphoid.pdf", path = "../results", width = 30, height = 30, units = "cm")

myelo_counts %>%
  filter(tissue_type == "Normal") %>%
  ggplot(aes(x = sample_id, fill = cell_type_imm)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = use_colors) +
  coord_flip() +
  scale_y_reverse()
ggsave2("SuppFig7A_myeloid.pdf", path = "../results", width = 30, height = 30, units = "cm")


###DEGs macrophages

imm_macro <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages1",
                                                            "Alveolar_Macrophages2",
                                                            "Alveolar_Macrophages3",
                                                            "CD14_Macrophages1",
                                                            "CD14_Macrophages2",
                                                            "CD14_Macrophages3",
                                                            "CD14_Macrophages4",
                                                            "CD14_Macrophages5",
                                                            "Macrophages_Proliferating"))

imm_macro <- ScaleData(imm_macro)

Idents(imm_macro) <- imm_macro@meta.data$cell_type_imm

macro_markers <- FindAllMarkers(imm_macro, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_macro_markers <- macro_markers %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

DoMultiBarHeatmap(imm_macro, features = top_macro_markers$gene, group.by = "cell_type_imm", additional.group.by = "tissue_type",additional.group.sort.by = "tissue_type", cols.use = list(tissue_type = use_colors), draw.lines = F) +
  scale_fill_viridis()
#ggsave2("HeatMap_Macro.pdf", path = "output/fig4", width = 30, height = 40, units = "cm")
ggsave2("SuppFig7B.png", path = "../results", width = 30, height = 40, units = "cm")
#ggsave2("HeatMap_Macro.emf", path = "output/fig4", width = 30, height = 40, units = "cm")


###DEGs T cells

imm_T <- subset(imm_anno, subset = cell_type_imm %in% c("T_conv1",
                                                        "T_conv2",
                                                        "T_reg",
                                                        "T_CD8_1",
                                                        "T_CD8_2",
                                                        "T_CD8_3",
                                                        "T_CD8_Proliferating"))

imm_T <- ScaleData(imm_T)

Idents(imm_T) <- imm_T@meta.data$cell_type_imm

markers_T <- FindAllMarkers(imm_T, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)

top_markers_T <- markers_T %>% group_by(cluster) %>% top_n(10, wt = avg_logFC)

DoMultiBarHeatmap(imm_T, features = top_markers_T$gene, group.by = "cell_type_imm", additional.group.by = "tissue_type",additional.group.sort.by = "tissue_type", cols.use = list(tissue_type = use_colors), draw.lines = F) +
  scale_fill_viridis()
#ggsave2("HeatMap_T.pdf", path = "output/fig4", width = 30, height = 35, units = "cm")
ggsave2("SuppFig7C.png", path = "../results", width = 30, height = 35, units = "cm")
#ggsave2("HeatMap_T.emf", path = "output/fig4", width = 30, height = 35, units = "cm")



###selected hallmark signatures and M1vsM2 signatures

m1m2_pws <- read_lines("../data/data/CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.gmt") %>%
  lapply(str_split, "//t") %>% 
  unlist(recursive = F) %>% 
  lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
  unlist(recursive = F)

m1m2_pws <- append(m1m2_pws, read_lines("../data/data/CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.gmt") %>%
                     lapply(str_split, "//t") %>% 
                     unlist(recursive = F) %>% 
                     lapply(function(x) setNames(list(x[-c(1:2)]), x[1])) %>% 
                     unlist(recursive = F))

imm_anno <- AddModuleScore(object = imm_anno, features = m1m2_pws, name = c("m1up", "m1dn"), nbin = 12)

VlnPlot(imm_anno, features = c("HALLMARK_INFLAMMATORY_RESPONSE31",
                               "HALLMARK_ALLOGRAFT_REJECTION46",
                               "HALLMARK_INTERFERON_GAMMA_RESPONSE19",
                               "HALLMARK_TNFA_SIGNALING_VIA_NFKB1",
                               "m1up1", 
                               "m1dn2"), 
        group.by = "cell_type_imm", pt.size = 0, ncol = 3, idents = c("Alveolar_Macrophages1",
                                                                      "Alveolar_Macrophages2",
                                                                      "Alveolar_Macrophages3",
                                                                      "CD14_Macrophages1",
                                                                      "CD14_Macrophages2",
                                                                      "CD14_Macrophages3",
                                                                      "CD14_Macrophages4",
                                                                      "CD14_Macrophages5",
                                                                      "Macrophages_Proliferating"), cols = use_colors)
ggsave2("Fig4C.pdf", path = "../results", width = 30, height = 20, units = "cm")


###T cell signatures

T_exhausted <- read_excel("../data/data/CD8_T_cells_exhausted.xlsx", skip = 1)
cytotoxicity <- c("PRF1", "IFNG", "GNLY", "NKG7", "GZMB", "GZMA", "GZMH", "KLRK1", "KLRB1", "KLRD1", "CTSW", "CST7")

T_cell_markers <- list(T_exhausted$GeneSymbol, cytotoxicity)

imm_T <- AddModuleScore(imm_T, features = T_cell_markers, name = c("exhaustion", "cytotoxicity"), nbin = 12)

VlnPlot(imm_T, features = c("cytotoxicity2", "exhaustion1"), pt.size = 0, group.by = "cell_type_imm", cols = use_colors, idents = c("T_CD8_1", "T_CD8_2", "T_CD8_3", "T_CD8_Proliferating"), ncol = 1)
ggsave2("Fig4F.pdf", path = "../results", width = 10, height = 20, units = "cm")


#############################################correlation cell types


### load libraries
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
library(cowplot)
library(gplots)
library(RColorBrewer)
library(plotrix)
library(corrplot)
library(pcaMethods)
library(readr)

theme_set(theme_cowplot())

#color scheme
use_colors <- c(
  Tumor = "brown2",
  Normal = "deepskyblue2",
  G1 = "#46ACC8",
  G2M = "#E58601",
  S = "#B40F20",
  Epithelial = "seagreen",
  Immune = "darkgoldenrod2",
  Stromal = "steelblue",
  p018 = "#E2D200",
  p019 = "#46ACC8",
  p023 = "#E58601",
  p024 = "#B40F20",
  p027 = "#0B775E",
  p028 = "#E1BD6D",
  p029 = "#35274A",
  p030 = "#F2300F",
  p031 = "#7294D4",
  p032 = "#5B1A18",
  p033 = "#9C964A",
  p034 = "#FD6467",
  lepidic = "#0B775E",
  acinar = "#74A089",
  `mucinuous (papillary)` = "#E2D200",
  `(micro)papillary` = "#CEAB07",
  solid = "#B40F20",
  sarcomatoid = "#5B1A18")


###load data and subsetting

#epi_anno <- readRDS("seurat_objects/epi_anno.RDS")

epi_tumor <- SubsetData(SubsetData(epi_anno, subset.name = "cluster_type", accept.value = "Tumor"), subset.name = "tissue_type", accept.value = "Tumor")
epi_tumor <- ScaleData(epi_tumor)

epi_pca <- epi_tumor
epi_pca <- RunPCA(epi_pca)



#imm_anno <- readRDS("seurat_objects/imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages1",
                                                                                         "Alveolar_Macrophages2",
                                                                                         "Alveolar_Macrophages3",
                                                                                         "CD14_Macrophages1",
                                                                                         "CD14_Macrophages2",
                                                                                         "CD14_Macrophages3",
                                                                                         "CD14_Macrophages4",
                                                                                         "CD14_Macrophages5",
                                                                                         "Macrophages_Proliferating",
                                                                                         "Monocytes",
                                                                                         "Myeloid_Dendritic",
                                                                                         "Plasmacytoid_Dendritic",
                                                                                         "Mast",
                                                                                         "T_conv1",
                                                                                         "T_conv2",
                                                                                         "T_reg",
                                                                                         "T_CD8_1",
                                                                                         "T_CD8_2",
                                                                                         "T_CD8_3",
                                                                                         "T_CD8_Proliferating",
                                                                                         "NK_cells",
                                                                                         "B_cells",
                                                                                         "Plasma"))

imm_lympho <- subset(imm_anno, subset = cell_type_imm %in% c("T_conv1",
                                                             "T_conv2",
                                                             "T_reg",
                                                             "T_CD8_1",
                                                             "T_CD8_2",
                                                             "T_CD8_3",
                                                             "T_CD8_Proliferating",
                                                             "NK_cells",
                                                             "B_cells",
                                                             "Plasma"))

imm_myelo <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages1",
                                                            "Alveolar_Macrophages2",
                                                            "Alveolar_Macrophages3",
                                                            "CD14_Macrophages1",
                                                            "CD14_Macrophages2",
                                                            "CD14_Macrophages3",
                                                            "CD14_Macrophages4",
                                                            "CD14_Macrophages5",
                                                            "Macrophages_Proliferating",
                                                            "Monocytes",
                                                            "Myeloid_Dendritic",
                                                            "Plasmacytoid_Dendritic",
                                                            "Mast"))

lympho_counts <- FetchData(imm_lympho, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

myelo_counts <- FetchData(imm_myelo, vars = c("tissue_type", "cell_type_imm", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal"))) 




#str_anno <- readRDS("seurat_objects/str_anno.RDS")

str_anno@meta.data$cell_type_str <- factor(str_anno@meta.data$cell_type_str, levels = c("Endothelial1",
                                                                                        "Endothelial2",
                                                                                        "Endothelial3",
                                                                                        "Endothelial4",
                                                                                        "Endothelial5",
                                                                                        "Endothelial6",
                                                                                        "Endothelial7",
                                                                                        "Lymphaticendothelial",
                                                                                        "Fibroblast1",
                                                                                        "Fibroblast2",
                                                                                        "Myofibroblast1",
                                                                                        "Myofibroblast2",
                                                                                        "Smoothmuscle1",
                                                                                        "Smoothmuscle2",
                                                                                        "Mesothelial"))

str_endo <- subset(str_anno, subset = cell_type_str %in% c("Endothelial1",
                                                           "Endothelial2",
                                                           "Endothelial3",
                                                           "Endothelial4",
                                                           "Endothelial5",
                                                           "Endothelial6",
                                                           "Endothelial7",
                                                           "Lymphaticendothelial"))

str_fibro <- subset(str_anno, subset = cell_type_str %in% c("Fibroblast1",
                                                            "Fibroblast2",
                                                            "Myofibroblast1",
                                                            "Myofibroblast2",
                                                            "Smoothmuscle1",
                                                            "Smoothmuscle2",
                                                            "Mesothelial"))


endo_counts <- FetchData(str_endo, vars = c("tissue_type", "cell_type_str", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal")))

fibro_counts <- FetchData(str_fibro, vars = c("tissue_type", "cell_type_str", "sample_id", "patient_id")) %>%  
  mutate(tissue_type = factor(tissue_type, levels = c("Tumor", "Normal"))) 


###count immune and stromal cells

myelo_counts_rel <- myelo_counts %>%
  filter(tissue_type == "Tumor") %>%
  dplyr::count(cell_type_imm, patient_id) %>%
  group_by(patient_id) %>% 
  mutate(n_rel = n/sum(n))


myelo_counts_rel <- myelo_counts_rel %>%
  pivot_wider(id_cols = patient_id, names_from = cell_type_imm, values_from = n_rel)


lympho_counts_rel <- lympho_counts %>%
  filter(tissue_type == "Tumor") %>%
  dplyr::count(cell_type_imm, patient_id) %>%
  group_by(patient_id) %>% 
  mutate(n_rel = n/sum(n)) 


lympho_counts_rel <- lympho_counts_rel %>%
  pivot_wider(id_cols = patient_id, names_from = cell_type_imm, values_from = n_rel)


fibro_counts_rel <- fibro_counts %>%
  filter(tissue_type == "Tumor") %>%
  dplyr::count(cell_type_str, patient_id) %>%
  group_by(patient_id) %>% 
  mutate(n_rel = n/sum(n))


fibro_counts_rel <- fibro_counts_rel %>%
  pivot_wider(id_cols = patient_id, names_from = cell_type_str, values_from = n_rel)


endo_counts_rel <- endo_counts %>%
  filter(tissue_type == "Tumor") %>%
  dplyr::count(cell_type_str, patient_id) %>%
  group_by(patient_id) %>% 
  mutate(n_rel = n/sum(n))


endo_counts_rel <- endo_counts_rel %>%
  pivot_wider(id_cols = patient_id, names_from = cell_type_str, values_from = n_rel)


cell_counts_rel <- full_join(myelo_counts_rel, lympho_counts_rel, by = "patient_id")
cell_counts_rel <- full_join(cell_counts_rel, endo_counts_rel, by = "patient_id")
cell_counts_rel <- full_join(cell_counts_rel, fibro_counts_rel, by = "patient_id")

cell_counts_rel[is.na(cell_counts_rel)] <- 0

###cell proportion statistics

test <- cell_counts_rel %>%
  mutate(pattern = ifelse(patient_id %in% c("p032", "p018", "p019", "p024", "p031", "p033"), "N3MC", "CP2E")) %>%
  mutate(T_CD8_exhausted = T_CD8_1 + T_CD8_2)

#for Kim et al. dataset
#test <- cell_counts_rel %>%
#  mutate(pattern = ifelse(patient_id %in% c("P0009", "P0018", "P0020", "P0028"), "CP2E", "N3MC")) %>%
#  mutate(T_CD8_exhausted = T_CD8_1 + T_CD8_2)

ggplot(test, aes(x = pattern, y = CD14_Macrophages1)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = CD14_Macrophages1~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_CD14_Macrophages1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = CD14_Macrophages2)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = CD14_Macrophages2~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_CD14_Macrophages2.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = Myeloid_Dendritic)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = Myeloid_Dendritic~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_Myeloid_Dendritic.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = Plasmacytoid_Dendritic)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = Plasmacytoid_Dendritic~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_Plasmacytoid_Dendritic.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = Myofibroblast1)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = Myofibroblast1~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_Myofibroblast1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = Myofibroblast2)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = Myofibroblast2~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_Myofibroblast2.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = T_conv1)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = T_conv1~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_T_conv1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = NK_cells)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = NK_cells~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_NK_cells.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(test, aes(x = pattern, y = T_CD8_exhausted)) +
  geom_boxplot() +
  ggtitle(paste0("p = ", wilcox.test(formula = T_CD8_exhausted~pattern, data = test, alternative = "two.sided", paired = F)$p.value))
ggsave2("SuppFig8A_Group1vs2_T_CD8_exhausted.pdf", path = "../results", width = 10, height = 10, units = "cm")



###select only cell types that occured in at least three patients

celltype_selected <- cell_counts_rel
celltype_selected[celltype_selected == 0] <- NA
celltype_selected <- cell_counts_rel[, which(colMeans(!is.na(celltype_selected)) >= 0.3)]

cell_counts_rel_selected <- cell_counts_rel[, colnames(cell_counts_rel) %in% colnames(celltype_selected)]

cell_counts_rel_mtrx <- as.matrix(cell_counts_rel_selected[,-1])
rownames(cell_counts_rel_mtrx) <- cell_counts_rel_selected$patient_id

###principal component analysis

cell_counts_pca <- as.data.frame(scores(pca(cell_counts_rel_mtrx, nPcs = 10)))
cell_counts_pca$patient_id <- rownames(cell_counts_pca)


histo_data <- unique(FetchData(epi_tumor, c("patient_id", "histo_subtype")))
cell_counts_pca <- full_join(cell_counts_pca, histo_data, by = "patient_id")

ggplot(cell_counts_pca, aes(x = PC1, y = PC2, color = histo_subtype, label = patient_id)) +
  geom_point(size = 2) +
  geom_text(aes(label = patient_id), color = "black", vjust = 1.5, size = 2) +
  scale_color_manual(values = use_colors) +
  scale_x_continuous(limits = c(-1, 1.2)) +
  scale_y_continuous(limits = c(-1.2, 0.8))
ggsave2("Fig5A.pdf", path = "../results", width = 15, height = 9, units = "cm")


###add mean tumor cell signature score

epi_pca <- AddModuleScore(epi_pca, features = list(c("SCGB3A1", "PIGR", "NAPSA", "C4BPA", "SCGB3A2", "HLA-DRA", "CD74", "ADGRF5", "C16orf89", "FOLR1", "SELENBP1", "HLA-DRB1", "ID4", "MGP", "AQP3", "CA2", "LHX9", "HLA-DPB1", "FMO5", "GKN2", "C5", "MUC1", "NPC2", "RNASE1", "PIK3C2G", "SFTA2", "SLC34A2", "HLA-DPA1", "FGFR3", "PGC")), name = "tumor_signature")
epi_pca <- AddModuleScore(epi_pca, features = list(c("DSG2", "CAMK2N1", "FAM3C", "KRT7", "IFI27", "SLC2A1", "MARCKS", "PLAU", "AHNAK2", "PERP", "S100A4", "KRT19", "COL6A1", "UACA", "COL17A1", "CDA", "TPM2", "S100A16", "KRT8", "PRSS23", "DST", "LAMC2", "S100P", "PRSS3", "LAMA3", "DSP", "ITGA3", "MDK", "FAM83A", "ITGB4")), name = "anti_tumor_signature")

epi_tumor_data <- FetchData(epi_pca, c(colnames(epi_pca@meta.data)))

epi_tumor_signature <- epi_tumor_data %>%
  group_by(patient_id) %>% 
  mutate(diff_mean = mean(tumor_signature1)) %>%
  mutate(undiff_mean = mean(anti_tumor_signature1)) %>%
  dplyr::select(patient_id, diff_mean, undiff_mean) %>%
  unique()

###heatmap of normalized immune and stromal cell counts, and mean tumor cell signatures

cell_counts_rel_heatmap <- full_join(cell_counts_rel_selected, cell_counts_pca, by = "patient_id")
cell_counts_rel_heatmap <- full_join(cell_counts_rel_heatmap, epi_tumor_signature, by = "patient_id")

cell_counts_rel_heatmap_ordered <- cell_counts_rel_heatmap[order(cell_counts_rel_heatmap$PC1),]

cell_counts_rel_heatmap1 <- cell_counts_rel_heatmap_ordered[,2:ncol(celltype_selected)]
cell_counts_rel_heatmap2 <- cell_counts_rel_heatmap_ordered[,c("diff_mean", "undiff_mean")]

rownames(cell_counts_rel_heatmap1) <- cell_counts_rel_heatmap_ordered$patient_id
rownames(cell_counts_rel_heatmap2) <- cell_counts_rel_heatmap_ordered$patient_id

cell_counts_rel_heatmap1 <- scale(cell_counts_rel_heatmap1)

pdf("../results/Fig5B_immune_and_stromal_cells.pdf", width = 7, height = 7)
heatmap.2(as.matrix(t(cell_counts_rel_heatmap1)), trace = "none", density.info = "none", scale = "none", symbreaks = F, col = hcl.colors(100, palette = "orrd", rev = T), Rowv = T, Colv = F, dendrogram = "row", margins = c(5,10))
dev.off()

pdf("../results/Fig5B_tumor_cell_signatures.pdf", width = 7, height = 5)
heatmap.2(as.matrix(t(cell_counts_rel_heatmap2)), trace = "none", density.info = "none", scale = "none", col = bluered(100), Rowv = F, Colv = F, dendrogram = "none")
dev.off()

###correlation plots

# matrix of the p-value of the correlation (from http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(cell_counts_rel_mtrx)

test <- cor(cell_counts_rel_mtrx, method = "spearman")

png("../results/SuppFig8A.png", width = 30, height = 30, units = "cm", res = 600)
corrplot(test, p.mat = p.mat, sig.level = 0.05, insig = "blank", type = "full", tl.col = "black")
dev.off()

#pdf("output/fig5/corrplot.pdf", width = 7, height = 7)
#corrplot(test, p.mat = p.mat, sig.level = 0.05, insig = "blank", type = "full", tl.col = "black")
#dev.off()

### correlation network plot

test_net <- test

for (i in 1:nrow(test_net)){
  for (j in 1:ncol(test_net)){
    if(p.mat[i,j] > 0.05){
      test_net[i,j] <- 0
    }
    if(test_net[i,j] == 1){
      test_net[i,j] <- 0
    }
  }
}

library(qgraph)

pdf("../results/Fig5C.pdf", width = 6, height = 6)
qgraph(test_net, layout = "spring", threshold = 0.7, labels = colnames(test_net), label.scale = F, label.cex = 0.7, theme = "colorblind", vsize = 3, color = "grey", borders = F)
dev.off()


#####scatter plots

cell_counts_rel <- mutate(cell_counts_rel, T_CD8_exhausted = T_CD8_1+T_CD8_2)

###plots stromal

ggplot(cell_counts_rel, aes(x = Myofibroblast1, y = Myofibroblast2)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.8)) +
  scale_y_continuous(limits = c(0,0.8)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$Myofibroblast1, cell_counts_rel$Myofibroblast2, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$Myofibroblast1, cell_counts_rel$Myofibroblast2, method = "spearman")$estimate))
ggsave2("Fig3E.pdf", path = "../results", width = 10, height = 9, units = "cm")


###plots immune

ggplot(cell_counts_rel, aes(x = CD14_Macrophages1, y = CD14_Macrophages2)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.6)) +
  scale_y_continuous(limits = c(0,0.8)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$CD14_Macrophages1, cell_counts_rel$CD14_Macrophages2, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$CD14_Macrophages1, cell_counts_rel$CD14_Macrophages2, method = "spearman")$estimate))
ggsave2("Fig4D_Macrophages1_2.pdf", path = "../results", width = 10, height = 9, units = "cm")

ggplot(cell_counts_rel, aes(x = CD14_Macrophages2, y = Plasmacytoid_Dendritic)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.8)) +
  scale_y_continuous(limits = c(0,0.06)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$CD14_Macrophages2, cell_counts_rel$Plasmacytoid_Dendritic, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$CD14_Macrophages2, cell_counts_rel$Plasmacytoid_Dendritic, method = "spearman")$estimate))
ggsave2("Fig4D_Macrophages2_pDC.pdf", path = "../results", width = 10, height = 9, units = "cm")

ggplot(cell_counts_rel, aes(x = T_reg, y = T_CD8_exhausted)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.5)) +
  scale_y_continuous(limits = c(0,0.85)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$T_reg, cell_counts_rel$T_CD8_exhausted, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$T_reg, cell_counts_rel$T_CD8_exhausted, method = "spearman")$estimate))
ggsave2("Fig4G_T_reg_exhausted.pdf", path = "../results", width = 10, height = 9, units = "cm")

ggplot(cell_counts_rel, aes(x = CD14_Macrophages2, y = T_CD8_exhausted)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.8)) +
  scale_y_continuous(limits = c(0,0.85)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$CD14_Macrophages2, cell_counts_rel$T_CD8_exhausted, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$CD14_Macrophages2, cell_counts_rel$T_CD8_exhausted, method = "spearman")$estimate))
ggsave2("Fig4G_Macrophages_exhausted.pdf", path = "../results", width = 10, height = 9, units = "cm")

ggplot(cell_counts_rel, aes(x = Plasmacytoid_Dendritic, y = T_CD8_exhausted)) +
  geom_point(aes(color = patient_id), size = 4) +
  scale_x_continuous(limits = c(0,0.06)) +
  scale_y_continuous(limits = c(0,0.85)) +
  scale_color_manual(values = use_colors) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(cell_counts_rel$Plasmacytoid_Dendritic, cell_counts_rel$T_CD8_exhausted, method = "spearman")$p.value, 
                 "/nr = ",cor.test(cell_counts_rel$Plasmacytoid_Dendritic, cell_counts_rel$T_CD8_exhausted, method = "spearman")$estimate))
ggsave2("Fig4G_pDC_exhausted.pdf", path = "../results", width = 10, height = 9, units = "cm")



###compare cell counts with IHC quantification

CD123 <- read_excel("../data/IHC_quantification/Quantification_CD123.xlsx")
CD123_scrnaseq <- cell_counts_rel %>% select(Plasmacytoid_Dendritic)
CD123 <- full_join(CD123, CD123_scrnaseq, by = "patient_id")
CD123 <- pivot_longer(CD123, cols = starts_with("image"), names_to = "image_nr", values_to = "count")
CD123_average <- CD123 %>% group_by(patient_id) %>% mutate(mean = mean(count)) %>% mutate(sd = sd(count)) %>% select(c(patient_id, Plasmacytoid_Dendritic, mean, sd)) %>% unique()
ggplot() +
  geom_jitter(data = CD123, mapping = aes(x = Plasmacytoid_Dendritic, y = count, color = patient_id), size = 0.5, alpha = 0.3) +
  geom_smooth(data = CD123_average, mapping = aes(x = Plasmacytoid_Dendritic, y = mean, color = patient_id), method = "lm", se = FALSE, size = 1, color = "black") +
  geom_errorbar(data = CD123_average, aes(x = Plasmacytoid_Dendritic, ymin=mean-sd, ymax=mean+sd, color = patient_id), size = 0.5, width = 0.002) +
  geom_point(data = CD123_average, aes(x = Plasmacytoid_Dendritic, y=mean, color = patient_id), shape=15, size = 3) +
  scale_color_manual(values = use_colors) +
  ggtitle(paste0("p = ",cor.test(CD123_average$mean, CD123_average$Plasmacytoid_Dendritic, method = "pearson")$p.value, 
                 "/nr = ",cor.test(CD123_average$mean, CD123_average$Plasmacytoid_Dendritic, method = "pearson")$estimate))
ggsave2("Fig4E_pDC_CD123.pdf", path = "../results", width = 10, height = 9, units = "cm")

CTHRC1 <- read_excel("../data/IHC_quantification/Quantification_CTHRC1.xlsx")
CTHRC1_scrnaseq <- cell_counts_rel %>% select(Myofibroblast2)
CTHRC1 <- full_join(CTHRC1, CTHRC1_scrnaseq, by = "patient_id")
CTHRC1 <- pivot_longer(CTHRC1, cols = starts_with("image"), names_to = "image_nr", values_to = "count")
CTHRC1_average <- CTHRC1 %>% group_by(patient_id) %>% mutate(mean = mean(count)) %>% mutate(sd = sd(count)) %>% select(c(patient_id, Myofibroblast2, mean, sd)) %>% unique()
ggplot() +
  geom_jitter(data = CTHRC1, mapping = aes(x = Myofibroblast2, y = count, color = patient_id), size = 0.5, alpha = 0.3) +
  geom_smooth(data = CTHRC1_average, mapping = aes(x = Myofibroblast2, y = mean, color = patient_id), method = "lm", se = FALSE, size = 1, color = "black") +
  geom_errorbar(data = CTHRC1_average, aes(x = Myofibroblast2, ymin=mean-sd, ymax=mean+sd, color = patient_id), size = 0.5) +
  geom_point(data = CTHRC1_average, aes(x = Myofibroblast2, y=mean, color = patient_id), shape=15, size = 3) +
  scale_color_manual(values = use_colors) +
  ggtitle(paste0("p = ",cor.test(CTHRC1_average$mean, CTHRC1_average$Myofibroblast2, method = "pearson")$p.value, 
                 "/nr = ",cor.test(CTHRC1_average$mean, CTHRC1_average$Myofibroblast2, method = "pearson")$estimate))
ggsave2("Fig3F_Myofibroblast2_CTHRC1.pdf", path = "../results", width = 10, height = 9, units = "cm")


CXCL9 <- read_excel("../data/IHC_quantification/Quantification_CXCL9.xlsx")
CXCL9_scrnaseq <- cell_counts_rel %>% select(CD14_Macrophages2)
CXCL9 <- full_join(CXCL9, CXCL9_scrnaseq, by = "patient_id")
CXCL9 <- pivot_longer(CXCL9, cols = starts_with("image"), names_to = "image_nr", values_to = "count")
CXCL9_average <- CXCL9 %>% group_by(patient_id) %>% mutate(mean = mean(count)) %>% mutate(sd = sd(count)) %>% select(c(patient_id, CD14_Macrophages2, mean, sd)) %>% unique()
ggplot() +
  geom_jitter(data = CXCL9, mapping = aes(x = CD14_Macrophages2, y = count, color = patient_id), size = 0.5, alpha = 0.3, width = 0.02) +
  geom_smooth(data = CXCL9_average, mapping = aes(x = CD14_Macrophages2, y = mean, color = patient_id), method = "lm", se = FALSE, size = 1, color = "black") +
  geom_errorbar(data = CXCL9_average, aes(x = CD14_Macrophages2, ymin=mean-sd, ymax=mean+sd, color = patient_id), size = 0.5, width = 0.03) +
  geom_point(data = CXCL9_average, aes(x = CD14_Macrophages2, y=mean, color = patient_id), shape=15, size = 3) +
  scale_color_manual(values = use_colors) +
  ggtitle(paste0("p = ",cor.test(CXCL9_average$mean, CXCL9_average$CD14_Macrophages2, method = "pearson")$p.value, 
                 "/nr = ",cor.test(CXCL9_average$mean, CXCL9_average$CD14_Macrophages2, method = "pearson")$estimate))
ggsave2("Fig4E_CD14_Macrophages_2_CXCL9.pdf", path = "../results", width = 10, height = 9, units = "cm")




########################################TCGA dataset


library(TCGAbiolinks)
library(GSVA)
library(SummarizedExperiment)
library(Seurat)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(corrplot)
library(stringr)
library(survival)
library(survminer)
library(plotrix)
library(readr)
library(HDF5Array)

theme_set(theme_cowplot())



###load gene expression data

luad_se <- loadHDF5SummarizedExperiment(dir="../data/TCGA LUAD")
luad_mut <- read_csv("../data/data/TCGA_LUAD_mutations.csv")



###replace gene symbols, calculate log2-transformed matrix

luad_data <- as.data.frame(assay(luad_se))

luad_data$hgnc_symbol <- luad_se@rowRanges$external_gene_name

luad_data <- luad_data[!(duplicated(luad_data$hgnc_symbol)),]
rownames(luad_data) <- luad_data$hgnc_symbol

luad_data <- luad_data[,1:(ncol(luad_data)-1)]
luad_data <- as.matrix(luad_data)

log_luad_data <- log2(luad_data+1)


###calculate cell cluster marker genes

#imm_anno <- readRDS("seurat_objects/imm_anno.RDS")

imm_anno@meta.data$cell_type_imm <- ordered(imm_anno@meta.data$cell_type_imm, levels = c("Alveolar_Macrophages1",
                                                                                         "Alveolar_Macrophages2",
                                                                                         "Alveolar_Macrophages3",
                                                                                         "CD14_Macrophages1",
                                                                                         "CD14_Macrophages2",
                                                                                         "CD14_Macrophages3",
                                                                                         "CD14_Macrophages4",
                                                                                         "CD14_Macrophages5",
                                                                                         "Macrophages_Proliferating",
                                                                                         "Monocytes",
                                                                                         "Myeloid_Dendritic",
                                                                                         "Plasmacytoid_Dendritic",
                                                                                         "Mast",
                                                                                         "T_conv1",
                                                                                         "T_conv2",
                                                                                         "T_reg",
                                                                                         "T_CD8_1",
                                                                                         "T_CD8_2",
                                                                                         "T_CD8_3",
                                                                                         "T_CD8_Proliferating",
                                                                                         "NK_cells",
                                                                                         "B_cells",
                                                                                         "Plasma"))

imm_lympho <- subset(imm_anno, subset = cell_type_imm %in% c("T_conv1",
                                                             "T_conv2",
                                                             "T_reg",
                                                             "T_CD8_1",
                                                             "T_CD8_2",
                                                             "T_CD8_3",
                                                             "T_CD8_Proliferating",
                                                             "NK_cells",
                                                             "B_cells",
                                                             "Plasma"))
imm_lympho <- ScaleData(imm_lympho)

imm_myelo <- subset(imm_anno, subset = cell_type_imm %in% c("Alveolar_Macrophages1",
                                                            "Alveolar_Macrophages2",
                                                            "Alveolar_Macrophages3",
                                                            "CD14_Macrophages1",
                                                            "CD14_Macrophages2",
                                                            "CD14_Macrophages3",
                                                            "CD14_Macrophages4",
                                                            "CD14_Macrophages5",
                                                            "Macrophages_Proliferating",
                                                            "Monocytes",
                                                            "Myeloid_Dendritic",
                                                            "Plasmacytoid_Dendritic",
                                                            "Mast"))
imm_myelo <- ScaleData(imm_myelo)

#str_anno <- readRDS("seurat_objects/str_anno.RDS")

str_anno@meta.data$cell_type_str <- factor(str_anno@meta.data$cell_type_str, levels = c("Endothelial1",
                                                                                        "Endothelial2",
                                                                                        "Endothelial3",
                                                                                        "Endothelial4",
                                                                                        "Endothelial5",
                                                                                        "Endothelial6",
                                                                                        "Endothelial7",
                                                                                        "Lymphaticendothelial",
                                                                                        "Fibroblast1",
                                                                                        "Fibroblast2",
                                                                                        "Myofibroblast1",
                                                                                        "Myofibroblast2",
                                                                                        "Smoothmuscle1",
                                                                                        "Smoothmuscle2",
                                                                                        "Mesothelial"))

str_endo <- subset(str_anno, subset = cell_type_str %in% c("Endothelial1",
                                                           "Endothelial2",
                                                           "Endothelial3",
                                                           "Endothelial4",
                                                           "Endothelial5",
                                                           "Endothelial6",
                                                           "Endothelial7",
                                                           "Lymphaticendothelial"))
str_endo <- ScaleData(str_endo)

str_fibro <- subset(str_anno, subset = cell_type_str %in% c("Fibroblast1",
                                                            "Fibroblast2",
                                                            "Myofibroblast1",
                                                            "Myofibroblast2",
                                                            "Smoothmuscle1",
                                                            "Smoothmuscle2",
                                                            "Mesothelial"))
str_fibro <- ScaleData(str_fibro)

Idents(imm_myelo) <- imm_myelo@meta.data$cell_type_imm
Idents(imm_lympho) <- imm_lympho@meta.data$cell_type_imm
Idents(str_endo) <- str_endo@meta.data$cell_type_str
Idents(str_fibro) <- str_fibro@meta.data$cell_type_str

myelo_markers <- FindAllMarkers(imm_myelo, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)
write_csv(myelo_markers, path = "../results/SuppTable9.csv")

lympho_markers <- FindAllMarkers(imm_lympho, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)
write_csv(lympho_markers, path = "../results/SuppTable10.csv")

endo_markers <- FindAllMarkers(str_endo, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)
write_csv(endo_markers, path = "../results/SuppTable7.csv")

fibro_markers <- FindAllMarkers(str_fibro, only.pos = T, min.pct = 0.25, min.diff.pct = 0.25)
write_csv(fibro_markers, path = "../results/SuppTable8.csv")

clusters_myelo <- as.character(unique(myelo_markers$cluster))
clusters_lympho <- as.character(unique(lympho_markers$cluster))
clusters_endo <- as.character(unique(endo_markers$cluster))
clusters_fibro <- as.character(unique(fibro_markers$cluster))

markers <- list()

for (i in seq_along(clusters_myelo)){
  markers[paste0(clusters_myelo[i])] <- myelo_markers %>% filter(cluster == clusters_myelo[i]) %>% dplyr::select(gene)
}

for (i in seq_along(clusters_lympho)){
  markers[paste0(clusters_lympho[i])] <- lympho_markers %>% filter(cluster == clusters_lympho[i]) %>% dplyr::select(gene)
}

for (i in seq_along(clusters_endo)){
  markers[paste0(clusters_endo[i])] <- endo_markers %>% filter(cluster == clusters_endo[i]) %>% dplyr::select(gene)
}

for (i in seq_along(clusters_fibro)){
  markers[paste0(clusters_fibro[i])] <- fibro_markers %>% filter(cluster == clusters_fibro[i]) %>% dplyr::select(gene)
}


### genes from PCA of tumor epithelial cells

markers[["PC1_pos"]] <- c("SCGB3A1", "PIGR", "NAPSA", "C4BPA", "SCGB3A2", "HLA-DRA", "CD74", "ADGRF5", "C16orf89", "FOLR1", "SELENBP1", "HLA-DRB1", "ID4", "MGP", "AQP3", "CA2", "LHX9", "HLA-DPB1", "FMO5", "GKN2", "C5", "MUC1", "NPC2", "RNASE1", "PIK3C2G", "SFTA2", "SLC34A2", "HLA-DPA1", "FGFR3", "PGC")
markers[["PC1_neg"]] <- c("DSG2", "CAMK2N1", "FAM3C", "KRT7", "IFI27", "SLC2A1", "MARCKS", "PLAU", "AHNAK2", "PERP", "S100A4", "KRT19", "COL6A1", "UACA", "COL17A1", "CDA", "TPM2", "S100A16", "KRT8", "PRSS23", "DST", "LAMC2", "S100P", "PRSS3", "LAMA3", "DSP", "ITGA3", "MDK", "FAM83A", "ITGB4")



###short signatures

for (i in seq_along(clusters_myelo)){
  markers[paste0(clusters_myelo[i], "_short")] <- myelo_markers %>% filter(cluster == clusters_myelo[i]) %>% top_n(5, wt = avg_logFC) %>% dplyr::select(gene)
}

for (i in seq_along(clusters_lympho)){
  markers[paste0(clusters_lympho[i], "_short")] <- lympho_markers %>% filter(cluster == clusters_lympho[i]) %>% top_n(5, wt = avg_logFC)%>% dplyr::select(gene)
}

for (i in seq_along(clusters_endo)){
  markers[paste0(clusters_endo[i], "_short")] <- endo_markers %>% filter(cluster == clusters_endo[i]) %>% top_n(5, wt = avg_logFC)%>% dplyr::select(gene)
}

for (i in seq_along(clusters_fibro)){
  markers[paste0(clusters_fibro[i], "_short")] <- fibro_markers %>% filter(cluster == clusters_fibro[i]) %>% top_n(5, wt = avg_logFC)%>% dplyr::select(gene)
}

markers[["PC1_pos_short"]] <- c("SCGB3A1", "PIGR", "NAPSA", "C4BPA", "SCGB3A2")
markers[["PC1_neg_short"]] <- c("DSG2", "CAMK2N1", "FAM3C", "KRT7", "IFI27")



###complete pattern signatures

#long
markers[["N3MC_complete"]] <- unique(c(markers[["CD14_Macrophages1"]], markers[["Myofibroblast1"]], markers[["NK_cells"]], markers[["Myeloid_Dendritic"]], markers[["T_conv1"]], markers[["PC1_pos"]]))
markers[["CP2E_complete"]] <- unique(c(markers[["CD14_Macrophages2"]], markers[["Myofibroblast2"]], markers[["Plasmacytoid_Dendritic"]], markers[["T_CD8_1"]], markers[["T_CD8_2"]], markers[["PC1_neg"]]))

#short
markers[["N3MC_complete_short"]] <- unique(c(markers[["CD14_Macrophages1_short"]], markers[["Myofibroblast1_short"]], markers[["NK_cells_short"]], markers[["Myeloid_Dendritic_short"]], markers[["T_conv1_short"]], markers[["PC1_pos_short"]]))
markers[["CP2E_complete_short"]] <- unique(c(markers[["CD14_Macrophages2_short"]], markers[["Myofibroblast2_short"]], markers[["Plasmacytoid_Dendritic_short"]], markers[["T_CD8_1_short"]], markers[["T_CD8_2_short"]], markers[["PC1_neg_short"]]))



###simple pattern signatures

#long
markers[["N3MC_simple"]] <- unique(c(markers[["CD14_Macrophages1"]], markers[["Myofibroblast1"]]))
markers[["CP2E_simple"]] <- unique(c(markers[["Myofibroblast2"]], markers[["PC1_neg"]]))

#short
markers[["N3MC_simple_short"]] <- unique(c(markers[["CD14_Macrophages1_short"]], markers[["Myofibroblast1_short"]]))
markers[["CP2E_simple_short"]] <- unique(c(markers[["Myofibroblast2_short"]], markers[["PC1_neg_short"]]))



###save gene sets
write_csv(as_data_frame(markers["N3MC_complete"]), path = "../results/SuppTable11_1.csv")
write_csv(as_data_frame(markers["CP2E_complete"]), path = "../results/SuppTable11_2.csv")
write_csv(as_data_frame(markers["N3MC_simple_short"]), path = "../results/SuppTable11_3.csv")
write_csv(as_data_frame(markers["CP2E_simple_short"]), path = "../results/SuppTable11_4.csv")



###single-sample Gene set Enrichment Analysis (ssGSEA)

luad_gsva <- gsva(expr = log_luad_data, gset.idx.list = markers, method = "ssgsea", verbose = TRUE, kcdf = "Gaussian")

luad_gsva <- as.data.frame(t(luad_gsva))


###calculate ratios

luad_gsva <- mutate(luad_gsva, PC1_ratio = PC1_pos / PC1_neg)
luad_gsva <- mutate(luad_gsva, Myofibro_ratio = Myofibroblast1 / Myofibroblast2)
luad_gsva <- mutate(luad_gsva, Macro_ratio = CD14_Macrophages1 / CD14_Macrophages2)
luad_gsva <- mutate(luad_gsva, Pattern_complete_ratio = N3MC_complete / CP2E_complete)
luad_gsva <- mutate(luad_gsva, Pattern_simple_ratio = N3MC_simple / CP2E_simple)

luad_gsva <- mutate(luad_gsva, PC1_ratio_short = PC1_pos_short / PC1_neg_short)
luad_gsva <- mutate(luad_gsva, Myofibro_ratio_short = Myofibroblast1_short / Myofibroblast2_short)
luad_gsva <- mutate(luad_gsva, Macro_ratio_short = CD14_Macrophages1_short / CD14_Macrophages2_short)
luad_gsva <- mutate(luad_gsva, Pattern_complete_short_ratio = N3MC_complete_short / CP2E_complete_short)
luad_gsva <- mutate(luad_gsva, Pattern_simple_short_ratio = N3MC_simple_short / CP2E_simple_short)


###add follow-up data

luad_gsva$days_to_death <- luad_se$days_to_death
luad_gsva$vital_status <- luad_se$vital_status
luad_gsva$days_to_last_follow_up <- luad_se$days_to_last_follow_up

luad_gsva$overall_survival <- ifelse(luad_gsva$vital_status == "Dead", luad_gsva$days_to_death, luad_gsva$days_to_last_follow_up)


###binarize data

luad_gsva <- mutate(luad_gsva, vital_status_bin = ifelse(vital_status == "Alive", 0, 1))

ES_names <- colnames(luad_gsva[1:98])

for (i in seq_along(ES_names)) {
  luad_gsva[,paste0(ES_names[i], "_bin")] <- ifelse(luad_gsva[,ES_names[i]] > median(luad_gsva[,ES_names[i]]), 1, 0)
}


###correlation plot all clusters

# matrix of the p-value of the correlation (from http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram)
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], method = "spearman", ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

p.mat <- cor.mtest(as.matrix(luad_gsva[,1:40]))

corr_tcga <- cor(as.matrix(luad_gsva[,1:40]), method = "spearman")

png(file = "../results/SuppFig8B.png", width = 30, height = 30, units = "cm", res = 600)
corrplot(corr_tcga, p.mat = p.mat, sig.level = 0.05, insig = "blank", type = "full", tl.col = "black")
dev.off()

#pdf(file = "output/fig5/corrplot_tcga.pdf", width = 7, height = 7)
#corrplot(corr_tcga, p.mat = p.mat, sig.level = 0.05, insig = "blank", type = "full", tl.col = "black")
#dev.off()


###correlation plot for selected cell clusters

### PC1 ~ Macrophages

ggplot(luad_gsva, aes(y = PC1_pos, x = CD14_Macrophages1)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_pos, luad_gsva$CD14_Macrophages1, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_pos, luad_gsva$CD14_Macrophages1, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_pos_CD14_Macrophages1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_pos, x = CD14_Macrophages2)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_pos, luad_gsva$CD14_Macrophages2, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_pos, luad_gsva$CD14_Macrophages2, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_pos_CD14_Macrophages2.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_neg, x = CD14_Macrophages1)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_neg, luad_gsva$CD14_Macrophages1, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_neg, luad_gsva$CD14_Macrophages1, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_neg_CD14_Macrophages1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_neg, x = CD14_Macrophages2)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_neg, luad_gsva$CD14_Macrophages2, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_neg, luad_gsva$CD14_Macrophages2, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_neg_CD14_Macrophages2.pdf", path = "../results", width = 10, height = 10, units = "cm")


### PC1 ~ Myofibroblasts

ggplot(luad_gsva, aes(y = PC1_pos, x = Myofibroblast1)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_pos, luad_gsva$Myofibroblast1, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_pos, luad_gsva$Myofibroblast1, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_pos_Myofibroblast1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_pos, x = Myofibroblast2)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_pos, luad_gsva$Myofibroblast2, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_pos, luad_gsva$Myofibroblast2, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_pos_Myofibroblast2.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_neg, x = Myofibroblast1)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_neg, luad_gsva$Myofibroblast1, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_neg, luad_gsva$Myofibroblast1, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_neg_Myofibroblast1.pdf", path = "../results", width = 10, height = 10, units = "cm")

ggplot(luad_gsva, aes(y = PC1_neg, x = Myofibroblast2)) +
  geom_point(color = "grey", size = 0.1) +
  geom_smooth(method = "lm", se = FALSE, size = 1, color = "black") +
  ggtitle(paste0("p = ",cor.test(luad_gsva$PC1_neg, luad_gsva$Myofibroblast2, method = "spearman")$p.value, 
                 "/nrho = ",cor.test(luad_gsva$PC1_neg, luad_gsva$Myofibroblast2, method = "spearman")$estimate))
ggsave2("Fig5E_PC1_neg_Myofibroblast2.pdf", path = "../results", width = 10, height = 10, units = "cm")


###Survival analysis

# Kaplan-Meyer and Log-rank

for (i in c(104:201)) {
  ggsurvplot(survfit(Surv(as.numeric(overall_survival/30.4375), vital_status_bin)~luad_gsva[[i]], data = luad_gsva), pval = T, pval.method = T, title = paste0(colnames(luad_gsva)[i]), xlim = c(0,120), break.x.by = 30, palette = c("Red", "Blue"))
  ggsave2(paste0("Fig5F_SuppFig11_", colnames(luad_gsva)[i],".pdf"), path = "../results", width = 8, height = 10, units = "cm")
}

# Cox regression

summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Macro_ratio_bin, data = luad_gsva)) 
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Myofibro_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ PC1_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Pattern_complete_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Pattern_simple_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Pattern_complete_short_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Pattern_simple_short_ratio_bin, data = luad_gsva))
summary(coxph(Surv(as.numeric(overall_survival), vital_status_bin) ~ Macro_ratio_bin + Myofibro_ratio_bin + PC1_ratio_bin, data = luad_gsva))



###add mutation data

#genes from nNGMv2 Panel
mutations <- c("ALK", "BRAF", "CTNNB1", "EGFR", "ERBB2", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "HRAS", "IDH1", "KEAP1", "KRAS", "MAP2K1", "MET", "NRAS", "NTRK1", "NTRK2", "NTRK3", "PIK3CA", "PTEN", "RET", "ROS1", "STK11", "TP53")

patients_mut <- list()

for (i in seq_along(mutations)){
  patients_mut[[i]] <- filter(luad_mut, Hugo_Symbol == mutations[i] & Variant_Classification != "Silent")$Tumor_Sample_Barcode %>% as.vector() %>% substr(1,12)
}

names(patients_mut) <- mutations

luad_gsva_mut <- luad_gsva

luad_gsva_mut[, mutations] <- 0
luad_gsva_mut$samples <- rownames(luad_gsva_mut) %>% substr(1,12)

luad_gsva_mut$Pattern_ratio_complete_bin <- ifelse(luad_gsva_mut$Pattern_complete_ratio_bin == 1, "N3MC_complete", "CP2E_complete")
luad_gsva_mut$Pattern_ratio_simple_short_bin <- ifelse(luad_gsva_mut$Pattern_simple_short_ratio_bin == 1, "N3MC_simple_short", "CP2E_simple_short")

luad_gsva_mut$TMB <- as.character(luad_se$paper_Nonsilent.Mutations.per.Mb) %>% gsub(pattern = ",", replacement = ".", fixed = T) %>% as.numeric()

for(i in seq_along(mutations)){
  luad_gsva_mut[, mutations[i]] <- ifelse(luad_gsva_mut$sample %in% patients_mut[[i]], "mut", "wt")
}

luad_gsva_mut <- filter(luad_gsva_mut, samples %in% substr(luad_mut$Tumor_Sample_Barcode, 1, 12))


#Pattern_ratio_bin

for(i in seq_along(mutations)){
  mut_counts <- luad_gsva_mut %>% group_by(Pattern_ratio_complete_bin) %>% dplyr::count(eval(parse(text = mutations[i]))) %>% pivot_wider(id_cols = Pattern_ratio_complete_bin, names_from = "eval(parse(text = mutations[i]))", values_from = n)
  
  mut_counts[is.na(mut_counts)] <- 0
  
  chi <- chisq.test(mut_counts[,2:3])
  
  ggplot(luad_gsva_mut) + 
    geom_bar(aes(x = Pattern_ratio_complete_bin, fill = eval(parse(text = mutations[i]))), position = "fill") +
    theme(legend.title = element_blank()) +
    ggtitle(paste0(mutations[i], " X² = ", chi$statistic, " p = ", chi$p.value))
  ggsave2(paste0("Fig5G_mut_", mutations[[i]], ".pdf"), path = "../results", width = 20, height = 10, units = "cm")
}

#for(i in seq_along(mutations)){
#  mut_counts <- luad_gsva_mut %>% group_by(Pattern_ratio_simple_short_bin) %>% dplyr::count(eval(parse(text = mutations[i]))) %>% pivot_wider(id_cols = Pattern_ratio_simple_short_bin, names_from = "eval(parse(text = mutations[i]))", values_from = n)
#  
#  mut_counts[is.na(mut_counts)] <- 0
#  
#  chi <- chisq.test(mut_counts[,2:3])
#  
#  ggplot(luad_gsva_mut) + 
#    geom_bar(aes(x = Pattern_ratio_simple_short_bin, fill = eval(parse(text = mutations[i]))), position = "fill") +
#    theme(legend.title = element_blank()) +
#    ggtitle(paste0(mutations[i], " X² = ", chi$statistic, " p = ", chi$p.value))
#  ggsave2(paste0("Fig5G_mut_simple_short_", mutations[[i]], ".pdf"), path = "../results", width = 20, height = 10, units = "cm")
#}

###Tumor mutational burden

ggplot(luad_gsva_mut) +
  geom_boxplot(aes(x = Pattern_ratio_complete_bin, y = TMB)) +
  ggtitle(paste0("p = ", t.test(formula = TMB~Pattern_ratio_complete_bin, data = luad_gsva_mut, alternative = "two.sided", paired = F)$p.value))
ggsave2("Fig5G_TMB.pdf", path = "../results", width = 20, height = 10, units = "cm")

#ggplot(luad_gsva_mut) +
#  geom_boxplot(aes(x = Pattern_ratio_simple_short_bin, y = TMB)) +
#  ggtitle(paste0("p = ", t.test(formula = TMB~Pattern_ratio_simple_short_bin, data = luad_gsva_mut, alternative = "two.sided", paired = F)$p.value))
#ggsave2("Fig5G_TMB_simple_short.pdf", path = "../results", width = 20, height = 10, units = "cm")


###add fusion data

fusions.data <- data.frame(luad_se$paper_Fusions, luad_se$barcode)
fusions.data$barcode <- fusions.data$luad_se.barcode

luad_gsva_mut$barcode <- rownames(luad_gsva_mut)

luad_gsva_fus <- inner_join(fusions.data, luad_gsva_mut, by = "barcode")

luad_gsva_fus$alk_fusions <- ifelse(grepl("ALK", luad_gsva_fus$luad_se.paper_Fusions) == T, "fus", "wt")
luad_gsva_fus$ret_fusions <- ifelse(grepl("RET", luad_gsva_fus$luad_se.paper_Fusions) == T, "fus", "wt")
luad_gsva_fus$ros_fusions <- ifelse(grepl("ROS", luad_gsva_fus$luad_se.paper_Fusions) == T, "fus", "wt")

fusions <- c("alk_fusions", "ret_fusions", "ros_fusions")

#Pattern_ratio_bin

for(i in seq_along(fusions)){
  fus_counts <- luad_gsva_fus %>% group_by(Pattern_ratio_complete_bin) %>% dplyr::count(eval(parse(text = fusions[i]))) %>% pivot_wider(id_cols = Pattern_ratio_complete_bin, names_from = "eval(parse(text = fusions[i]))", values_from = n)
  
  fus_counts[is.na(fus_counts)] <- 0
  
  chi <- chisq.test(fus_counts[,2:3])
  
  ggplot(luad_gsva_fus) + 
    geom_bar(aes(x = Pattern_ratio_complete_bin, fill = eval(parse(text = fusions[i]))), position = "fill") +
    theme(legend.title = element_blank()) +
    ggtitle(paste0(fusions[i], " X? = ", chi$statistic, " p = ", chi$p.value))
  ggsave2(paste0("Fig5G_fus_", fusions[[i]], ".pdf"), path = "../results", width = 20, height = 10, units = "cm")      
}

#for(i in seq_along(fusions)){
#  fus_counts <- luad_gsva_fus %>% group_by(Pattern_ratio_simple_short_bin) %>% dplyr::count(eval(parse(text = fusions[i]))) %>% pivot_wider(id_cols = Pattern_ratio_simple_short_bin, names_from = "eval(parse(text = fusions[i]))", values_from = n)
#  
#  fus_counts[is.na(fus_counts)] <- 0
#  
#  chi <- chisq.test(fus_counts[,2:3])
#  
#  ggplot(luad_gsva_fus) + 
#    geom_bar(aes(x = Pattern_ratio_simple_short_bin, fill = eval(parse(text = fusions[i]))), position = "fill") +
#    theme(legend.title = element_blank()) +
#    ggtitle(paste0(fusions[i], " X? = ", chi$statistic, " p = ", chi$p.value))
#  ggsave2(paste0("Fig5_fus_simple_short_", fusions[[i]], ".pdf"), path = "../results", width = 20, height = 10, units = "cm")      
#}


print(nrow(luad_gsva))
print(nrow(luad_gsva_mut))
luad_gsva_mut %>% filter(TMB != "NA") %>% nrow() %>% print()
