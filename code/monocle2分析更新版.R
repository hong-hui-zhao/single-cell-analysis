setwd("C:/Users/ZHH/Desktop/Mar")
library(Seurat)
library(monocle)
library(dplyr)
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
  } 
  return(monocle_cds)
}
mouse_data <- subset(mouse_data, celltype %in% c('Malignant1','Malignant2','Malignant3','AT2_NO_G12D'))
mouse_monocle_fun <- seurat_to_monocle(mouse_data, assay = "RNA", slot = "counts")

#===========================================================================

#Estimate size factors and dispersions
mouse_monocle <- estimateSizeFactors(mouse_monocle_fun)
mouse_monocle <- estimateDispersions(mouse_monocle)


#质量控制
mouse_monocle <- detectGenes(mouse_monocle, min_expr = 0.1)#至少10%的表达，这一步完成之后，会出现num_cells_expressed这一列
print(head(fData(mouse_monocle)))
print(head(pData(mouse_monocle)))
expressed_genes <- row.names(subset(fData(mouse_monocle), num_cells_expressed >= 10))
pData(mouse_monocle)$Total_mRNAs <- Matrix::colSums(exprs(mouse_monocle))#将UMI加入cds
mouse_monocle <- mouse_monocle[,pData(mouse_monocle)$Total_mRNAs < 1e6]#阈值按照自己实际情况

#differentialGeneTest
cds_DGT <- mouse_monocle
diff_test_res <- differentialGeneTest(cds_DGT, fullModelFormulaStr = "~celltype")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)

cds_DGT <- reduceDimension(cds_DGT, max_components = 2, reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT)
cds_DGT <- orderCells(cds_DGT, root_state = NULL, num_paths = NULL, reverse = T)

GM_state <- function(cds){
   if (length(unique(pData(cds)$State)) > 1){
     T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"AT2"]
     return(as.numeric(names(T0_counts)[which
                                        (T0_counts == max(T0_counts))]))
   } else {
     return (1)
   }
}
 
cds_DGT <- orderCells(cds_DGT, root_state = 3)

#################################################################################
#初步可视化
plot_cell_trajectory(cds_DGT, cell_size = 2.2, color_by = "celltype")
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT, color_by = "orig.ident",cell_link_size = 1.5)
plot_cell_trajectory(cds_DGT, color_by = "State")

#基因表达拟时图
plot_cell_trajectory(cds_DGT, markers = "KRAS")
#各个分组拟时图
plot_cell_trajectory(cds_DGT, color_by = "orig.ident") + facet_wrap(~orig.ident, nrow = 2)



#################################################################################
#查看随着pseudotime变化的基因
cds_DGT_pseudotimegenes <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~sm.ns(Pseudotime)")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
write.csv(cds_DGT_pseudotimegenes,"cds_DGT_pseudotimegenes.csv")

#筛选一些基因去做图，这里完全是为了作图，实际展示基因按照自己按要求
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 1e-320)
Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds_DGT[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T)


saveRDS(cds_DGT,file = 'cds_DGT.rds')
save.image('epiother_pesudotime.RData')
#################################################################################
#查看在分叉处的差异基因
BEAM_res <- BEAM(cds_DGT, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
write.csv(BEAM_res,file = "BEAM_res.csv")
plot_genes_branched_heatmap(mouse_monocle[row.names(subset(BEAM_res, qval < 1e-120)),],
                            branch_point = 1,
                            num_clusters = 6,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

# colorRampPalette(c("steelblue", "white","darkorange2", "orangered4"))(100)