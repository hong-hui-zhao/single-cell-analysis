
setwd("/home/tq_ziv/data_analysis/WGCNA/单细胞WGCNA/")

# install.packages(c("Seurat", "WGCNA", "igraph", "devtools", "GeneOverlap))
devtools::install_github('smorabit/hdWGCNA', ref='dev')
library(hdWGCNA)
library(WGCNA)
library(Seurat) 
library(tidyverse) 
library(igraph)
library(cowplot)
library(patchwork)
library(dplyr)
library(ggplot2)
library(stringr)
theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 8)

# DimPlot(scRNA, label = T)
# Bra_rna <- subset(scRNA, celltype=="Vas", invert=T)
# Bra_rna@meta.data$celltype = droplevels(Bra_rna@meta.data$celltype, 
#                                        exclude = setdiff(levels(Bra_rna@meta.data$celltype),
#                                                          unique(Bra_rna@meta.data$celltype)))

DimPlot(Bra_rna, label = T)
DefaultAssay(Bra_rna) <- 'RNA'

#===========================================================================================
#                          1、Pre-precessing
#===========================================================================================

#1、Create a slot in a Seurat object to store hdWGCNA data
# This function also select genes for WGCNA analysis

#(选择哪些基因用于分析，取决于你自己的研究预期，没有优劣之分，切记！！！)

# three approches for select genes using gene_select parameter:
### variable: 使用存储在Seurat对象的VariableFeatures中的基因,就是高变基因
### fraction: 使用在整个数据集或每组细胞中特定部分细胞中表达的基因
### custom：使用自定义列表中指定的基因，就是自己选定基因集

Bra_rna <- SetupForWGCNA(Bra_rna,
                       wgcna_name = "Bra",
                       #Select genes that will be used for co-expression network analysis
                       gene_select = "fraction",
                       fraction = 0.05#基因表达百分比，0.05表示这个基因必在5%的细胞中表达
                       )


#2、Construct metacells
#简单来说，因为单细胞矩阵实稀疏矩阵，原始矩阵很大，而且不太适用于WGCNA分析
#这里作者k-Nearest Neighbors (KNN)算法将相似的细胞聚合，然后计算他们的平均基因表达量
#这样得到的矩阵稀疏性会大大降低。用这个矩阵代替原始矩阵进行分析！

# dim(Bra_rna)
# [1] 26830 19468

Bra_rna <- MetacellsByGroups(
  seurat_obj =  Bra_rna,
  group.by = c("celltype","orig.ident"),
  k = 25,#需要聚合的细胞数，一般20-75
  max_shared=20,
  reduction = 'umap',
  ident.group = 'celltype')


#Normalize metacell expression matrix:
#因为我们MetacellsByGroups默认的slot用的count，所以需要标准化，
#如果一开始使用data，那么这一步不用跑
Bra_rna <- NormalizeMetacells(Bra_rna)


# DefaultAssay(Bra_rna) <- "integrated"
# seurat_obj <- NormalizeMetacells(Bra_rna)
# seurat_obj <- ScaleMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# seurat_obj <- RunPCAMetacells(seurat_obj, features=VariableFeatures(seurat_obj))
# seurat_obj <- RunHarmonyMetacells(seurat_obj, group.by.vars='orig.ident')
# seurat_obj <- RunUMAPMetacells(seurat_obj, reduction='harmony', dims=1:20)
# 
# 
# p1 <- DimPlotMetacells(seurat_obj, group.by='celltype') + umap_theme() + ggtitle("Cell Type")
# p2 <- DimPlotMetacells(seurat_obj, group.by='orig.ident') + umap_theme() + ggtitle("Sample")
# 
# p1 | p2

#===========================================================================================
#                          2、Co-expression network analysis
#===========================================================================================
#这里我们纳入所有的celltype进行分析，所以在SetDatExpr中没有设置group_name。
#如果你需要关注某一个或者多个celltype，那么group_name设置你需要的celltype即可
#Set up the expression matrix
Bra_rna <- SetDatExpr(Bra_rna, assay = 'RNA', slot = 'data')


#Select soft-power threshold选择软阈值
#其他参数参照WGCNA包的pickSoftThreshold函数。
Bra_rna <- TestSoftPowers(Bra_rna, networkType = 'signed') # you can also use "unsigned" or "signed hybrid"


# plot the results:
plot <- PlotSoftPowers(Bra_rna)

# assemble with patchwork
#可视化软阈值
wrap_plots(plot, ncol=2)



#Construct co-expression network
#构建网络,其他参数参照WGCNA包的blockwiseConsensusModules
Bra_rna <- ConstructNetwork(
  Bra_rna,
  soft_power = 4,
  tom_name = "Bra_Test",
  setDatExpr=F)


#Module identification
#可视化tree
PlotDendrogram(Bra_rna, main='scRNA hdWGCNA Dendrogram')


#inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(Bra_rna)


#===========================================================================================
#                          3、compute module eigengenes and connectivity
#===========================================================================================

# Compute harmonized module eigengenes
Bra_rna <- ScaleData(Bra_rna)
Bra_rna <- ModuleEigengenes(Bra_rna, group.by.vars = 'orig.ident')
## Harmonize module eigengenes:
hMEs <- GetMEs(Bra_rna)
# module eigengenes:
MEs <- GetMEs(Bra_rna, harmonized=FALSE)


#Compute module connectivity
Bra_rna <- ModuleConnectivity(Bra_rna)

# rename the modules
Bra_rna <- ResetModuleNames(Bra_rna,new_name = "M")


# get the module table
modules <- GetModules(Bra_rna)


#修改module颜色

# get a table of just the module and it's unique color
mod_color_df <- GetModules(Bra_rna) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

#不要grey模块
n_mods <- nrow(mod_color_df) - 1

# reset the module colors(我们这里有9个module)
newcolor <- c("#f4c40f","#fe9b00","#d8443c","#de597c","#e87b89","#633372","#1f6e9c","#2b9b81","#92c051")
Bra_rna <- ResetModuleColors(Bra_rna, newcolor)


#重新getmodule
modules <- GetModules(Bra_rna)
write.csv(modules, file = 'modules.csv')


#重新做一下基因聚类图
PlotDendrogram(Bra_rna, main='scRNA hdWGCNA Dendrogram')


# get hub genes
hub_df <- GetHubGenes(Bra_rna, n_hubs = 25)
write.csv(hub_df, file = 'hub_df.csv')


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method

library(UCell)
Bra_rna <- ModuleExprScore(
  Bra_rna,
  n_genes = 25,
  method='UCell'
)

#保存文件
save(Bra_rna, file = 'Bra_rna_wgcna.RData')
#===========================================================================================
#                          4、Visualize the network  
#===========================================================================================

###plot1

# plot genes ranked by kME for each module
#按照kME排序，可视化每个module的基因
PlotKMEs(Bra_rna, ncol=3)


###plot2
#Plot module eigengenes--using hMEs
plot_hMEs <- ModuleFeaturePlot(
  Bra_rna,
  reduction = "umap",
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  raster = T
)


pdf("5_ModuleFeaturePlot_hMEs.pdf",width=10,height=8)

wrap_plots(plot_hMEs, ncol=3)

dev.off()



#Plot module eigengenes--using ucell score
plot_score <- ModuleFeaturePlot(
  Bra_rna,
  reduction = "umap",
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  raster = T,
  ucell = TRUE
)


pdf("5_ModuleFeaturePlot_score.pdf",width=10,height=8)

wrap_plots(plot_score, ncol=3)

dev.off()


###plot3

## add MEvalues to metadata
Bra_rna@meta.data <- cbind(
  Bra_rna@meta.data,
  GetMEs(Bra_rna, harmonized=TRUE)
)

MEs <- GetMEs(Bra_rna, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# plot with Seurat's DotPlot function

DotPlot(Bra_rna, features=mods, group.by = 'metacell_grouping')+coord_flip()+
  theme_bw()+theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))



###plot4

p1 <- VlnPlot(Bra_rna,
              features = c('M2'),
              group.by = 'celltype',
              pt.size = 0)
p1= p1+geom_boxplot(width=.25, fill='white')

## Change axis labels and remove legend:
p1 <- p1 + xlab('') + ylab('hME') + NoLegend()
## Plot output
p1



###plot5

ModuleNetworkPlot(Bra_rna,outdir = "./Bra_ModuleNetworks_final")


###plot6

#all module
# hubgene network

pdf("8-HubGeneNetworkPlot.pdf",width=8,height=5)

HubGeneNetworkPlot(
  Bra_rna,
  n_hubs = 10, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

dev.off()




###plot7

Bra_rna <- RunModuleUMAP(
  Bra_rna,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding 
  n_neighbors=10, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
  #supervised=TRUE,target_weight=0.5 
  
)

# # get the hub gene UMAP table from the seurat object
# umap_df <- GetModuleUMAP(Bra_rna)
# 
# # plot with ggplot
# ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
#   geom_point(
#     color=umap_df$color, # color each point by WGCNA module
#     size=umap_df$kME*2 # size of each point based on intramodular connectivity
#   ) +
#   umap_theme()

pdf("9_ModuleUMAPPlot_fianl.pdf",width=8,height=8)

ModuleUMAPPlot(
  Bra_rna,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=3 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE
  
)

dev.off()



#===========================================================================================
#                          5、enrichment analysis 
#===========================================================================================
#官网给的示例使用enrichR包进行富集分析。
#我们都已经拿到module基因的文件了，那么自己提取基因用其他包进行分析也是可以的

library(enrichR)
# enrichr databases to test，都2024了，那就用新的版本
dbs <- c('GO_Biological_Process_2023')

# perform enrichment tests
Bra_rna <- RunEnrichr(
  Bra_rna,
  dbs=dbs, # character vector of enrichr databases to test
  max_genes = 100 # number of genes per module to test
)

# retrieve the output table
enrich_df <- GetEnrichrTable(Bra_rna)
write.csv(enrich_df,"enrich_GO.csv")



# enrichr dotplot
pdf("10_module_GO.pdf",width = 8,height = 8)
EnrichrDotPlot(
  Bra_rna,
  mods = "all",
  database = "GO_Biological_Process_2023", # this has to be one of the lists we used above!!!
  n_terms=3 # number of terms for each module
)

dev.off()
save.image('Monocyte_WGCNA.RData')











