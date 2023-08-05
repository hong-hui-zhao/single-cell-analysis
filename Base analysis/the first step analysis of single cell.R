###################################

### the first step analysis of single cell 

### author honghui Zhao et al

###################################

rm(list = ls())

library(Seurat)
library(tidyverse)
library(scater)
library(patchwork)
library(readxl)

setwd("~")

### Load data-----------------
### load metadata and add into samples 
samples <- read_excel("patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample_id

### import cellranger files from different data sets
for (i in seq_along(samples)){
  assign(paste0("scs_data", i), Read10X(data.dir = paste0(samples[i])))
}

### create seurat objects from cellranger files
for (i in seq_along(samples)){
  assign(paste0("seu_obj", i), CreateSeuratObject(counts = eval(parse(text = paste0("scs_data", i))), project = samples[i],
                                                  min.cells = 3,min.features = 200))
}
### merge data
sce <- merge(seu_obj1, y = c(seu_obj2, seu_obj3,seu_obj4),
             add.cell.ids = samples, project = "lung")

### calculate mitochondrial, hemoglobin and ribosomal gene counts
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
sce[["percent.HB"]] <- PercentageFeatureSet(sce, pattern = "^HBA|^HBB")
sce[["percent.rb"]] <- PercentageFeatureSet(sce, pattern = "^RPS|^RPL") 

### Data visualization
print(VlnPlot(sce, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb",
                                "percent.HB"), ncol = 3))

# ggsave("vlnplot_before_QC.pdf", path = "../results", width = 30, height = 20, units = "cm")
# VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, 
#         ncol = 1, log = T)

### saveRDS(sce,file = "sce_before_QC.rds")
# sce <- readRDS('sce_before_QC.rds')

### Amount of cell data prior to quality control
PreQCNumCells <- length(colnames(sce)) #22408

### Add metadata ------------------
metatable <- read_excel("patients_metadata.xlsx")
metadata <- FetchData(sce, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample_id <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample_id")
rownames(metadata) <- metadata$cell_id
sce <- AddMetaData(sce, metadata = metadata)

### QC quality control  -------------------------
### QC metrics: nCount_sce, nFeature_sce, and percent mitochrondial counts, hemoglobin and ribosomal gene counts
### Outliers are > 3 MADs
### log1p() = log(1 + number)
sce@meta.data$nCount_RNA_outlier_mad   <- isOutlier(log(sce@meta.data$nCount_RNA),
                                                    log = F,type = "both",nmads = 3)
sce@meta.data$nFeature_RNA_outlier_mad <- isOutlier(log(sce@meta.data$nFeature_RNA),
                                                    log = F,type = "both",nmads = 3)
sce@meta.data$percent_mt_outlier_mad   <- isOutlier(log1p(sce@meta.data$percent.mt),
                                                    log = F,type = "both",nmads = 3) 
sce@meta.data$percent_HB_outlier_mad   <- isOutlier(log1p(sce@meta.data$percent.HB),
                                                    log = F,type = "both",nmads = 3)
sce@meta.data$percent_rb_outlier_mad   <- isOutlier(log1p(sce@meta.data$percent.rb),
                                                    log = F,type = "both",nmads = 3)

sce <- subset(sce, subset = nCount_RNA_outlier_mad == 'FALSE'   &
                            nFeature_RNA_outlier_mad == 'FALSE' & 
                            percent_mt_outlier_mad == 'FALSE'   &
                            percent_HB_outlier_mad == 'FALSE'   &
                            percent_rb_outlier_mad == 'FALSE')

# sce <- subset(sce,subset= nCount_RNA < 40000 & nFeature_RNA < 7000 & percent.mt < 30)

# print(VlnPlot(sce, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),
#              ncol = 2),group.by = "patient.id" ) 
# ggsave("vln_after_QC.pdf", path = "C:/Users/ZHH/Desktop/results/QC", width = 30, height = 20, units = "cm")
### saveRDS(sce,file = "sce_after_QC.rds")

sce <- readRDS('sce_after_QC.rds')

### Amount of cell data after quality control
PostQCNumCells <- length(colnames(sce)) #25607

### Data visualization after quality control -------------------
### Data normalization 
sce <- NormalizeData(sce) %>%
       FindVariableFeatures(sce, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(sce)
sce <- sce %>%
       ScaleData( features = all.genes) %>%
       RunPCA(sce, npcs = 50) %>%
       FindNeighbors(sce, dims = 1:30)%>%
       FindClusters(sce)%>%
       RunUMAP(sce, dims = 1:30)%>%
       RunTSNE(sce, dims = 1:30)%>%

DimPlot(object = sce, reduction = "umap", pt.size = .1)
ggsave("umap_before_battch.pdf", path = "~", width = 30, height = 20, units = "cm")
### batch effect -----------------------
library(harmony)
sce <- sce %>% 
           RunHarmony(group.by.vars = "orig.ident", plot_convergence = T,assay.use = T) %>% 
           RunUMAP(reduction = "harmony", dims = 1:30) %>% 
           FindNeighbors(reduction = "harmony", dims = 1:30)

DimPlot(sce, reduction = "umap", pt.size = .1)
ggsave("umap_after_battch.pdf", path = "~", width = 30, height = 20, units = "cm")

library(clustree)
sce <- FindClusters(sce, resolution = c(0.1,0.3,0.5,0.7,0.9,1.0))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
colnames(sce@meta.data)
ggsave("clustree_all.pdf", path = "~", width = 30, height = 20, units = "cm")
### saveRDS(sce,file = "sce_after_processing.rds")

### the first step analysis of single cell was end --------------------------

