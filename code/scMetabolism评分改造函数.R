# devtools::install_github("YosefLab/VISION")
# devtools::install_github("YosefLab/VISION@v2.1.0")
setwd('~/KRAS/Monocyte/metabolism')

library(scMetabolism)
library(Seurat)
library(ggplot2)
library(rsvd)
library(pheatmap)


human_countexp_Seurat<-sc.metabolism.Seurat(obj = human_data,
                                            method = "AUCell", 
                                            imputation =F, 
                                            ncores = 2, 
                                            metabolism.type = "KEGG")

T_metabolism <- subset(human_countexp_Seurat, celltype=='T cell')
df = data.frame(t(T_metabolism@assays[["METABOLISM"]][["score"]]))#提取细胞代谢通路score，即使用了subset提取的score还是所有细胞的，所以接下来提取T细胞的即可
names(T_metabolism$orig.ident)#可以看到，两者细胞命名不一样，所以先修改一致（自己视实际情况）
rownames(df) <- gsub(".", "-", rownames(df), fixed = TRUE)
df = df[names(T_metabolism$orig.ident),]#提取T细胞，因为我们这个数据有重复，到这里的这个矩阵其实还可以进行limma做差异分析，这里就不再演示了
df$orig.ident <- T_metabolism$orig.ident#添加上分组
#计算每个样的平均值做热图，所有细胞展现也可以，但是可能差异效果看不出来
avg_df =aggregate(df[,1:ncol(df)-1],list(df$orig.ident),mean)
rownames(avg_df) = avg_df$Group.1
avg_df=avg_df[,-1]

avg_df_sec <- as.data.frame(t(avg_df))
 
pdf('metabolism_Monocyte.pdf',width = 12,height = 12)
pheatmap(avg_df_sec, show_colnames = T,scale='row', cluster_rows = T,
         color=colorRampPalette(c('#1A5592','white',"#B83D3D"))(100),
         cluster_cols = T)
dev.off()


DotPlot.metabolism(obj = T_metabolism, pathway = rownames(T_metabolism@assays[["METABOLISM"]][["score"]])[1:20], 
                   phenotype = "orig.ident", norm = "y")




