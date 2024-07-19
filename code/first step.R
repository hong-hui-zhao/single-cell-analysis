setwd('/home/zhh/KRAS/rawdata')
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
library(DoubletFinder)
library(decontX)
library(harmony)
library(cowplot)
library(clustree)
library(readxl)

nFeature_lower <- 500
nFeature_upper <- 6000
nCount_lower <- 1000
nCount_upper <- 40000
percent.mt_lower <- 0
percent.mt_upper <- 15
percent.HB_lower <- 0
percent.HB_upper <- 5

### Load data and data processings -----------------
samples <- read_excel("patients_metadata.xlsx", range = cell_cols("A:A")) %>% .$sample

for (i in seq_along(samples)){
  seu_obj <- CreateSeuratObject(counts = Read10X(data.dir = paste0(samples[i])), 
                                project = samples[i],
                                min.cells = 3,
                                min.features = 200)
  seu_obj[["percent.mt"]] <- PercentageFeatureSet(seu_obj, pattern = "^MT-")
  seu_obj[["percent.rb"]] <- PercentageFeatureSet(seu_obj, pattern = "^RPS|^RPL") 
  seu_obj <- custom_filter(seu_obj)
  seu_obj <- process_sample(seu_obj)
  seu_obj <- Find_doublet(seu_obj)
  seu_obj <- subset(seu_obj, subset = (DF.classify == "Singlet"))
  seu_obj <- CreateSeuratObject(
             counts = GetAssayData(seu_obj, assay = "RNA"),
             project = samples[i],
             meta.data = seu_obj@meta.data)
  seu_obj <- runDecontX(seu_obj)
  seu_obj <- seu_obj[,seu_obj$estConp < 0.3]
  assign(paste0("seurat_object_", i), seu_obj)
  
}


sce<- merge(seurat_object_1, y = c(seurat_object_2,seurat_object_3,seurat_object_4,seurat_object_5,seurat_object_6,seurat_object_7,
                                   seurat_object_8,seurat_object_9,seurat_object_10,seurat_object_11,seurat_object_12,seurat_object_13,
                                   seurat_object_14,seurat_object_15,seurat_object_16,seurat_object_17,seurat_object_18,seurat_object_19,seurat_object_20),
             add.cell.ids = samples, project = "lung")

rm(seurat_object_1,seurat_object_2,seurat_object_3,seurat_object_4,seurat_object_5,seurat_object_6,seurat_object_7,
   seurat_object_8,seurat_object_9,seurat_object_10,seurat_object_11,seurat_object_12,seurat_object_13,
   seurat_object_14,seurat_object_15,seurat_object_16,seurat_object_17,seurat_object_18,seurat_object_19,seurat_object_20,seu_obj)

sce <- NormalizeData(sce, normalization.method = "LogNormalize") 
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(sce, vars.to.regress = c("percent.mt", "percent.rb"),verbose = T)
sce <- RunPCA(sce, npcs = 50, verbose = FALSE)

Seurat::ElbowPlot(sce, ndims = 50)

sce <- sce %>% 
  RunHarmony(group.by.vars = "orig.ident", plot_convergence = T,assay.use = 'T')%>% 
  RunUMAP(reduction = "harmony", dims = 1:10) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:10)


sce <- FindClusters(sce, resolution = c(0.1,0.3,0.5,0.6,0.7,0.8,0.9,1.0))
clustree(sce@meta.data, prefix = "RNA_snn_res.")
ggsave("clustree_all.pdf", width = 8, height = 8)


p1 <- qc_std_plot(sce)
ggsave("feature_before_QC.pdf", p1, path = path, width = 20, height = 10)


mainmarkers <-c("PECAM1", "VWF", "ACTA2",#stromal
                "JCHAIN", "MS4A1", "PTPRC", "CD68", "KIT", #imm
                "EPCAM", "CDH1", "KRT7", "KRT19"#epi
)
cellmarkers <- c('CD19','CD79A','MS4A1',#B
                 'CD68','CD163','APOE','MARCO','MRC1',#MARCO
                 'TPSAB1','KIT','CPA3',#Mast
                 'CD14','FCN1',#Monocyte
                 'LYZ','CLEC9A', 'LAMP3',#Myeloid
                 'NKG7','KLRD1',#NK
                 'MZB1','IGHG1','JCHAIN',#Plasma
                 'CD3D','CD3','CD3E','CD8','CD4',#T
                 'PECAM1','VWF','CD31',#Endothlial
                 'EPCAM','KRT18','KRT19','CDH1','SFTPC',#Epithelial
                 'PDGFRA','COL1A1','COL1A2','ACTA2',#Fibroblast
                 'CD15','CSF3R','FCGR3B','S100A9'#Neutrophil
)
DotPlot(sce, features =cellmarkers, group.by = "RNA_snn_res.0.6") + 
  coord_flip() + 
  scale_color_viridis()

metatable <- read_excel("patients_metadata.xlsx")
metadata <- FetchData(sce, "orig.ident")
metadata$cell_id <- rownames(metadata)
metadata$sample <- metadata$orig.ident
metadata <- left_join(x = metadata, y = metatable, by = "sample")
rownames(metadata) <- metadata$cell_id
sce <- AddMetaData(sce, metadata = metadata)


sce <- SetIdent(sce,value = "RNA_snn_res.0.6")
annotation_curated_main <- read_excel("curated_annotation_main.xlsx")
new_ids_main <- annotation_curated_main$celltype
names(new_ids_main) <- levels(sce)
sce <- RenameIdents(sce, new_ids_main)
sce@meta.data$celltype <- Idents(sce)
Idents(sce) <- sce@meta.data$celltype

DimPlot(sce,label = T,raster=FALSE) + NoLegend()
saveRDS(sce,file = "sce_after_processing.rds")

### the first step analysis of single cell was end --------------------------
custom_filter <- function(seu_obj){
  seu_obj@meta.data$nCount_RNA_outlier_mad   <- isOutlier(log(seu_obj@meta.data$nCount_RNA),
                                                          log = F,type = "both",nmads = 3)
  seu_obj@meta.data$nFeature_RNA_outlier_mad <- isOutlier(log(seu_obj@meta.data$nFeature_RNA),
                                                          log = F,type = "both",nmads = 3)
  seu_obj@meta.data$percent_mt_outlier_mad   <- isOutlier(log1p(seu_obj@meta.data$percent.mt),
                                                          log = F,type = "both",nmads = 3) 
  seu_obj@meta.data$percent_rb_outlier_mad   <- isOutlier(log1p(seu_obj@meta.data$percent.rb),
                                                          log = F,type = "both",nmads = 3)
  
  seu_obj <- subset(seu_obj, subset = nCount_RNA_outlier_mad == 'FALSE'   &
                      nFeature_RNA_outlier_mad == 'FALSE' & 
                      percent_mt_outlier_mad == 'FALSE'   &
                      percent_rb_outlier_mad == 'FALSE'   & 
                      nCount_RNA > 1000 & nFeature_RNA > 500 & percent.mt < 20 & percent.rb < 40)
  return(seu_obj)}

process_sample <- function(sample) {
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
  sample <- ScaleData(sample, vars.to.regress = c("percent.mt", "percent.rb"), verbose = FALSE)
  sample <- RunPCA(sample, npcs = 50, verbose = FALSE)
  sample <- RunUMAP(sample, reduction = 'pca',dims = 1:20)
  sample <- FindNeighbors(sample,reduction = 'pca', dims = 1:20)
  sample <- FindClusters(sample,resolution = c(0.1,0.3,0.5,0.6,0.7,0.8,0.9,1.0))
  sample <- RunTSNE(sample, dims = 1:20)

}

Find_doublet <- function(data){
  sweep.res.list <- paramSweep(data, PCs=1:20, sct=F, num.cores=20)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  homotypic.prop <- modelHomotypic(data@meta.data$seurat_cluster)
  nExp_poi <- round(ncol(data)*8*1e-6 * length(Cells(data)))
  nExp_poi  <- round(nExp_poi*(1-homotypic.prop))
  data <- doubletFinder(data, PCs = 1:20, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = F)
  pann <- grep(pattern="^pANN", x=names(data@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(data@meta.data), value=TRUE)
  data$pANN <- data[[pann]]
  data$DF.classify <- data[[classify]]
  data[[pann]] <- NULL
  data[[classify]] <- NULL
  return(data)
}

runDecontX <- function(seurat_obj){
  counts <- seurat_obj@assays$RNA@counts
  x <- counts[rowSums(counts)>0,]
  decon <- decontX(x,verbose=TRUE, seed=1)
  seurat_obj[["origCounts"]] <- CreateAssayObject(counts = counts)
  newCounts <- decon$decontXcounts#这就是校正后的矩阵
  newCounts <- rbind(newCounts, counts[rowSums(counts)==0,])[rownames(counts),]
  seurat_obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  seurat_obj$estConp <- decon$contamination 
  return(seurat_obj)
}

qc_std_plot_helper <- function(x) x + 
  scale_color_viridis() +
  geom_point(size = 0.01, alpha = 0.3)

qc_std_plot <- function(seu_obj) {
  qc_data <- as_tibble(FetchData(seu_obj, c("nCount_RNA", "nFeature_RNA", "percent.mt",  "percent.rb")))
  plot_grid(
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = percent.mt))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),

    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), log2(nFeature_RNA), color = percent.rb))) + 
      geom_hline(yintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_hline(yintercept = log2(nFeature_upper), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), percent.mt, color = nFeature_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),

    qc_std_plot_helper(ggplot(qc_data, aes(log2(nCount_RNA), percent.rb, color = nFeature_RNA))) + 
      geom_vline(xintercept = log2(nCount_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nCount_upper), color = "red", linetype = 2),
    
    
    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), percent.mt, color = nCount_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),

    qc_std_plot_helper(ggplot(qc_data, aes(log2(nFeature_RNA), percent.rb, color = nCount_RNA))) + 
      geom_vline(xintercept = log2(nFeature_lower), color = "red", linetype = 2) +
      geom_vline(xintercept = log2(nFeature_upper), color = "red", linetype = 2),
    
    qc_std_plot_helper(ggplot(qc_data, aes(percent.rb, percent.mt, color = nCount_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2),
    qc_std_plot_helper(ggplot(qc_data, aes(percent.rb, percent.mt, color = nFeature_RNA))) + 
      geom_hline(yintercept = percent.mt_lower, color = "red", linetype = 2) +
      geom_hline(yintercept = percent.mt_upper, color = "red", linetype = 2),
    
    
    ggplot(gather(qc_data, key, value), aes(key, value)) +
      geom_violin() +
      facet_wrap(~key, scales = "free", ncol = 5),
    
    ncol = 3, align = "hv"
  )
}
