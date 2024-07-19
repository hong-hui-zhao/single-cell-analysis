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
library(harmony)
library(RColorBrewer)

theme_set(theme_cowplot())
set.seed(12345)
enableWGCNAThreads(nThreads = 20)

DefaultAssay(sce) <- 'RNA'
sce <- SetupForWGCNA(sce,
                         wgcna_name = "sce",
                         #Select genes that will be used for co-expression network analysis
                         gene_select = "fraction",
                         fraction = 0.05#基因表达百分比，0.05表示这个基因必在5%的细胞中表达
)

sce <- MetacellsByGroups(
  seurat_obj =  sce,
  group.by = c("celltype","mutation"),
  k = 25,#需要聚合的细胞数，一般20-75
  max_shared=20,
  reduction = 'umap',
  ident.group = 'celltype')




#这里我们纳入所有的celltype进行分析，所以在SetDatExpr中没有设置group_name。
#如果你需要关注某一个或者多个celltype，那么group_name设置你需要的celltype即可
#Set up the expression matrix
sce <- SetDatExpr(sce, assay = 'RNA', slot = 'data')



sce <- TestSoftPowers(sce, networkType = 'signed') # you can also use "unsigned" or "signed hybrid"

plot <- PlotSoftPowers(sce)
pdf("2_wrap_plots_hMEs.pdf",width=10,height=8)
# assemble with patchwork
#可视化软阈值
wrap_plots(plot, ncol=2)
dev.off()

#Construct co-expression network
#构建网络,其他参数参照WGCNA包的blockwiseConsensusModules
sce <- ConstructNetwork(
  sce,
  soft_power = 4,
  tom_name = "sce_Test",
  overwrite_tom = TRUE,
  setDatExpr=F)


#Module identification
#可视化tree
PlotDendrogram(sce, main='scRNA hdWGCNA Dendrogram')


#inspect the topoligcal overlap matrix (TOM)
TOM <- GetTOM(sce)


#===========================================================================================
#                          3、compute module eigengenes and connectivity
#===========================================================================================

# Compute harmonized module eigengenes
sce <- ScaleData(sce)
sce <- ModuleEigengenes(sce, group.by.vars = 'mutation')

## Harmonize module eigengenes:
hMEs <- GetMEs(sce)
# module eigengenes:
MEs <- GetMEs(sce, harmonized=FALSE)


#Compute module connectivity
sce <- ModuleConnectivity(sce)

# rename the modules
sce <- ResetModuleNames(sce,new_name = "M")


# get the module table
modules <- GetModules(sce)


#修改module颜色

# get a table of just the module and it's unique color
mod_color_df <- GetModules(sce) %>%
  dplyr::select(c(module, color)) %>%
  distinct %>% arrange(module)

#不要grey模块
n_mods <- nrow(mod_color_df) - 1

# reset the module colors(我们这里有9个module)
newcolor <- c("#00FFFF", "#7FFFD4", "#8EE5EE", "#7FFF00",  "#BF3EFF", "#FF3030")#,# "#FF69B4", "#00FF00", "#CDC673", "#FF00FF", "#9B30FF", "#C6E2FF", "#FFFF00", "#8A2BE2", "#CAFF70", "#CDAD00", "#66CDAA", "#FF83FA", "#FFAEB9")
sce <- ResetModuleColors(sce, newcolor)


#重新getmodule
modules <- GetModules(sce)
write.csv(modules, file = 'modules.csv')

pdf("2_scRNA hdWGCNA Dendrogram.pdf",width=10,height=8)
#重新做一下基因聚类图
PlotDendrogram(sce, main='scRNA hdWGCNA Dendrogram')
dev.off()

# get hub genes
hub_df <- GetHubGenes(sce, n_hubs = 25)
write.csv(hub_df, file = 'hub_df.csv')


# compute gene scoring for the top 25 hub genes by kME for each module
# with UCell method

library(UCell)
sce <- ModuleExprScore(
  sce,
  n_genes = 25,
  method='UCell'
)


#===========================================================================================
#                          4、Visualize the network  
#===========================================================================================

###plot1

# plot genes ranked by kME for each module
#按照kME排序，可视化每个module的基因
pdf("2_PlotKMEs.pdf",width=9,height=9)
PlotKMEs(sce, ncol=3)
dev.off()

###plot2
#Plot module eigengenes--using hMEs
plot_hMEs <- ModuleFeaturePlot(
  sce,
  reduction = "umap",
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  raster = T
)


pdf("5_ModuleFeaturePlot_hMEs.pdf",width=9,height=9)

wrap_plots(plot_hMEs, ncol=3)

dev.off()



#Plot module eigengenes--using ucell score
plot_score <- ModuleFeaturePlot(
  sce,
  reduction = "umap",
  features='hMEs', # plot the hMEs
  order=TRUE,# order so the points with highest hMEs are on top
  raster = T,
  ucell = TRUE
)


pdf("5_ModuleFeaturePlot_score.pdf",width=9,height=9)

wrap_plots(plot_score, ncol=3)

dev.off()


###plot3

## add MEvalues to metadata
sce@meta.data <- cbind(
  sce@meta.data,
  GetMEs(sce, harmonized=TRUE)
)

MEs <- GetMEs(sce, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']
mods <- c('M1','M2','M3','M4','M5','M6')#,'M7','M8','M9','M10','M12','M13','M14','M15','M16','M17','M18','M19')
# plot with Seurat's DotPlot function
pdf("5_dot_module.pdf",width=9,height=6)
DotPlot(sce, features=mods, group.by = 'metacell_grouping')+coord_flip()+
  theme_bw()+theme(axis.title = element_blank(),
                   axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

dev.off()

###plot4
# Loop through each mod and generate VlnPlot
for (mod in mods) {
  p <- VlnPlot(
    sce,
    features = mod,
    group.by = 'celltype',
    pt.size = 0
  )
  p <- p + geom_boxplot(width = 0.25, fill = 'white')
  p <- p + xlab('') + ylab('hME') + NoLegend()
  
  # Save plot as PDF
  ggsave(paste0("module_", mod, ".pdf"), plot = p, width = 5, height = 5)
}



###plot5

ModuleNetworkPlot(sce,outdir = "./Bra_ModuleNetworks_final")


###plot6

#all module
# hubgene network

pdf("8-HubGeneNetworkPlot.pdf",width=5,height=5)

HubGeneNetworkPlot(
  sce,
  n_hubs = 5, n_other=3,
  edge_prop = 0.75,
  mods = 'all'
)

dev.off()




###plot7

sce <- RunModuleUMAP(
  sce,
  n_hubs = 10, # number of hub genes to include for the UMAP embedding 
  n_neighbors=10, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
  #supervised=TRUE,target_weight=0.5 
  
)

# # get the hub gene UMAP table from the seurat object
# umap_df <- GetModuleUMAP(sce)
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
  sce,
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
# Define databases to test
dbs <- c('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023')

# Loop through each database and perform enrichment tests
for (db in dbs) {
  # Run enrichment analysis
  sce <- RunEnrichr(
    sce,
    dbs = db, # Use the current database in the loop
    max_genes = 100 # number of genes per module to test
  )
  
  # Retrieve the output table
  enrich_df <- GetEnrichrTable(sce)
  write.csv(enrich_df, paste0("enrich_", gsub("_", "", db), ".csv")) # Modify file name
  
  # Enrichr dotplot
  dotplot <- EnrichrDotPlot(
    sce,
    mods = "all",
    database = db, # Use the current database in the loop
    n_terms = 3 # number of terms for each module
  )
  
  # Save dotplot as PDF
  ggsave(filename = paste0("10_module_", gsub("_", "", db), ".pdf"), plot = dotplot, width = 15, height = 12)

}

saveRDS(sce,'T and NK_anno.RDS')
save.image('T and NK_wgcna.RData')









