library(Seurat)
library(monocle)
library(dplyr)
#===========================================================================
setwd('~/KRAS/Epi')
mouse_data

mouse_data <- subset(mouse_data,idents = c('Malignant3','Malignant4','AT2_G12D'))
mouse_data[["origCounts"]] <- NULL
mouse_expr<-as.sparse(mouse_data[["RNA"]]@counts)
# #分组信息
mouse_idents<-mouse_data@meta.data
# #基因信息
mouse_features <-as.matrix(mouse_expr[,1])
mouse_features[,1] <-row.names(mouse_features)
colnames(mouse_features)[1] <-c("gene_short_name")
mouse_features<-data.frame(mouse_features)
# 
# 
# #创建monocle文件 
mouse_idents <- new("AnnotatedDataFrame", data = mouse_idents)
mouse_features <- new("AnnotatedDataFrame", data = mouse_features)
mouse_monocle <- newCellDataSet(as.matrix(mouse_expr),
                                 phenoData = mouse_idents, 
                                 featureData = mouse_features, 
                                 expressionFamily=negbinomial.size())


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


#选差异基因(可以使用多个方法),主要是用于后期order cell，可自定义
#====================method1====================================================
#monocle计算高变基因
cds_monocle <- mouse_monocle
disp_table <- dispersionTable(cds_monocle)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)#meanexpression筛选标准按自己情况
head(unsup_clustering_genes)
cds_monocle <- setOrderingFilter(cds_monocle, unsup_clustering_genes$gene_id)
save(cds_monocle, file = 'cds_monocle.RData')
#====================method2====================================================
#使用seurat计算
cds_seurat <- mouse_monocle
seurat_variable_genes <- VariableFeatures(FindVariableFeatures(mouse_data, assay = "RNA"), assay = "RNA")
cds_seurat <- setOrderingFilter(cds_seurat, seurat_variable_genes)
plot_ordering_genes(cds_DGT)
save(cds_seurat, file = 'cds_seurat.RData')

#====================method3====================================================
#differentialGeneTest
cds_DGT <- mouse_monocle
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~celltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]#或者选择前多少的基因
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)
save.image("epi_pesudotime.RData")

#################################################################################
#接下来我们以第三种方法进行
#降维,排序
cds_DGT <- reduceDimension(cds_DGT, max_components = 2,reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT)
cds_DGT <- orderCells(cds_DGT, root_state = NULL, num_paths = NULL, reverse = T) 


# #假设我们已经知道根了，就可以指定，这里我指定细胞群
# GM_state <- function(cds){
#   if (length(unique(pData(cds)$State)) > 1){
#     T0_counts <- table(pData(cds)$State, pData(cds)$celltype)[,"PMN(7)"]
#     return(as.numeric(names(T0_counts)[which
#                                        (T0_counts == max(T0_counts))]))
#   } else {
#     return (1)
#   }
# }
# 
# cds_DGT <- orderCells(cds_DGT, root_state = GM_state(cds_DGT))

#################################################################################
#初步可视化
plot_cell_trajectory(mouse_monocle, cell_size = 2.2, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2)
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT, color_by = "orig.ident",cell_link_size = 1.5)
plot_cell_trajectory(cds_DGT, color_by = "State")

#基因表达拟时图
plot_cell_trajectory(cds_DGT, markers = "Il1b",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
#各个分组拟时图
plot_cell_trajectory(cds_DGT, color_by = "orig.ident") + facet_wrap(~orig.ident, nrow = 2)



#################################################################################
#查看随着pseudotime变化的基因
cds_DGT_pseudotimegenes <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~sm.ns(Pseudotime)")
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 0.01)
write.csv(cds_DGT_pseudotimegenes_sig,"cds_DGT_pseudotimegenes_sig.csv")

#筛选一些基因去做图，这里完全是为了作图，实际展示基因按照自己按要求
cds_DGT_pseudotimegenes_sig <- subset(cds_DGT_pseudotimegenes, qval < 1e-320)
Time_genes <- cds_DGT_pseudotimegenes_sig %>% pull(gene_short_name) %>% as.character()
p <- plot_pseudotime_heatmap(cds_DGT[Time_genes,], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T)



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