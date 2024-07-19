library(infercnv)
library(Seurat)
library(SeuratObject)
library(Rtsne)
library(RColorBrewer)
library(ComplexHeatmap)

setwd("~/KRAS/Plasma")
data <- readRDS('PandB.RDS')
process_sample <- function(sample) {
  sample <- NormalizeData(sample)
  sample <- FindVariableFeatures(sample, selection.method = "vst", nfeatures = 3000)
  sample <- ScaleData(sample, vars.to.regress = c("percent.mt", "percent.rb"), verbose = FALSE)
  sample <- RunPCA(sample, npcs = 50, verbose = FALSE)
  sample <- RunUMAP(sample, reduction = 'pca',dims = 1:20)
  sample <- FindNeighbors(sample,reduction = 'pca', dims = 1:20)
  sample <- FindClusters(sample,resolution = c(0.1,0.3,0.5,0.6,0.7,0.8,0.9,1.0))
  sample <- RunTSNE(sample, dims = 1:20,check_duplicates = FALSE)
  
}
data <- process_sample(data)

matrix_counts <- as.matrix(data[["RNA"]]@counts)
data$newcelltype <- paste0(data$sample,"_",data$celltype)
unique(data$newcelltype)
write.table(data$newcelltype, "celltype.label.txt", sep = "\t", quote = F, col.names = F)
url <- "https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt"
local_file <- "hg38_gencode_v27.txt"
download.file(url, local_file)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts,#count矩阵
                                    annotations_file="./celltype.label.txt",#celltype信息
                                    delim="\t",
                                    gene_order_file="./hg38_gencode_v27.txt",
                                    ref_group_names=c("Normal1_B Cell","Normal2_B Cell"),
                                    chr_exclude=c("chrY", "chrM"))#选择不需要的染色体，查看帮助函数去除
# cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="infercnv",  #分析结果out文件夹名
                             no_prelim_plot = T,
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             min_cells_per_gene = 20,
                             num_threads=20,#线程数
                             write_expr_matrix = T,
                             write_phylo = FALSE)

expr <- infercnv_obj@expr.data

#top注释-染色体
gene_pos <- read.delim("hg38_gencode_v27.txt", header = F)
gene_pos <- gene_pos[gene_pos$V1 %in% rownames(expr),]

new_cluster <- unique(gene_pos$V2)

top_color <- HeatmapAnnotation(cluster = anno_block(labels = gsub("chr", "", new_cluster),
                                                    gp = gpar(col = "white"),
                                                    labels_gp = gpar(cex = 1, col = "black"),
                                                    height = unit(5,"mm"))) 



#提取reference cell位置
ref_cell <- c(infercnv_obj@reference_grouped_cell_indices$`Normal1_B Cell`,
              infercnv_obj@reference_grouped_cell_indices$`Normal2_B Cell`)

obs_cell <- c(infercnv_obj@observation_grouped_cell_indices$`G12D1_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`G12D2_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`G12D3_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`G12D4_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D1_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D2_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D3_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D4_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D5_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N2_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N3_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N4_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1N1_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1N2_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1N3_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1N4_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`N1N5_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Normal1_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`Normal2_Plasma`,
              infercnv_obj@observation_grouped_cell_indices$`G12D1_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`G12D2_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`G12D3_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`G12D4_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D1_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D2_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D3_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D4_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`Non_G12D5_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N2_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N3_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N4_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1N1_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1N2_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1N3_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1N4_B Cell`,
              infercnv_obj@observation_grouped_cell_indices$`N1N5_B Cell`)

cell_anno =data.frame(cell_id =c(colnames(expr)[ref_cell],colnames(expr)[obs_cell]),
                      group =c(rep("ref",length(ref_cell)),rep("obs",length(obs_cell))))



#CNV聚类，这里使用kmeans聚类，聚类多少个cluster是自己选择的，按照自己的情况调整
#是一个灵活的参数,这里的目的主要是将具有明显CNV的细胞聚在一起，这样方便我们识别哪些是恶性细胞
set.seed(123)
library(tidyverse)
kmeans.result <- kmeans(t(expr), 15)
kmeans_df <- data.frame(kmeans.result$cluster)
colnames(kmeans_df) <- "k_cluster"
kmeans_df <- as_tibble(cbind(cell_id = rownames(kmeans_df), kmeans_df))
kmeans_df=kmeans_df%>%inner_join(cell_anno,by="cell_id")%>%arrange(k_cluster)
kmeans_df$k_cluster=as.factor(kmeans_df$k_cluster) 

 #注释
annotation_row = data.frame(
  k_cluster =kmeans_df$k_cluster,
  group = kmeans_df$group
  
)
row.names(annotation_row) <- kmeans_df$cell_id

library(circlize)
color_cluster=c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CDAA7D", "#00CDCD", "#00FFFF", "#1E90FF", "#FF3030")
names(color_cluster)=as.character(1:10)
left_anno <- rowAnnotation(df = annotation_row,
                           col=list(group=c("N"="#00A0877F","T" = "#E64B357F"),
                                    k_cluster=color_cluster),
                           show_annotation_name = F)


#绘制热图
pdf("CNV_heatmap.pdf",width = 15,height = 10)
Heatmap(t(log2(expr))[rownames(annotation_row),],
             col = colorRamp2(c(-0.5,0,0.5), c("#2166ac","white","#b2182b")),
             cluster_rows = F,
             cluster_columns = F,
             show_column_names = F,
             show_row_names = F,
             column_split = factor(gene_pos$V2, new_cluster),
             heatmap_legend_param = list(title = "inferCNV",
                                         direction = "vertical",
                                         title_position = "leftcenter-rot",
                                         legend_height = unit(3, "cm")),
             
             left_annotation = left_anno, 
             row_title = NULL,
             column_title = NULL,
             top_annotation = top_color,
             border = T)
draw(ht, heatmap_legend_side = "right")
dev.off()



#method2------------------------------------------------------------------------
#计算CNV score

#这里CNV score的计算基本都是参考这篇文章方法的描述：Single-cell RNA-seq highlights intra-tumoral heterogeneity and malignant progression in pancreatic ductal adenocarcinoma

CNV_score=as.data.frame(colMeans((expr-1)^2))
colnames(CNV_score)="CNV_score"
CNV_score$cell_id=rownames(CNV_score)
CNV_score=CNV_score%>%inner_join(kmeans_df,by="cell_id")#结合分组和kemans聚类结果
library(ggplot2)
ggplot(CNV_score, aes(k_cluster,CNV_score))+
  geom_violin(aes(fill=k_cluster),color="black")+
  scale_fill_manual(values = c("#98FB98", "#FF83FA", "#FFFF00", "#C0FF3E", "#CDC673","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#CDAA7D", "#00CDCD", "#00FFFF", "#1E90FF", "#FF3030"))+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
               geom = "pointrange", color = "black", size=0.5)

#可以看到，cluster3，5有较高的CNV score
ggsave('CNV_score.pdf',width = 8,height = 5)
#我们将分组和cluster结合。这个示例数据确实有点问题，HC组的细胞CNV还挺高
CNV_score$newcluster <- paste0(CNV_score$group,"_",CNV_score$k_cluster)
ggplot(CNV_score, aes(newcluster,CNV_score))+
  geom_violin(aes(fill=newcluster),color="black")+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  stat_summary(fun.data = "mean_sdl",  fun.args = list(mult = 1), 
               geom = "pointrange", color = "black", size=0.5)

ggsave('CNV_score_group.pdf',width =15,height = 5 )
#method3------------------------------------------------------------------------
#CNV score和CNV correlation
#这也是文献中的展示方式：Microenvironment drives cell state, plasticity, and drug response in pancreatic cancer
obs <- read.table("./infercnv/infercnv.observations.txt", header=T, check.names=F)
ref <- read.table("./infercnv/infercnv.references.txt", header=T, check.names=F)
estimateCNV <- function(obs,#读入的obs矩阵
                        ref, #读入的ref矩阵
                        score_threshold, #CNV_score阈值，大于多少定义为Malignant cell
                        #个人人为这个参数阈值有很大的主观性
                        cor_threshold#相关性阈值，也是一个主观性参数
){
  cnv_obs <- colMeans((obs-1)^2)
  cnv_ref <- colMeans((ref-1)^2) 
  cell_top <- names(sort(cnv_obs, decreasing=T))[1:round(length(cnv_obs))]
  cnv_top <- rowMeans(obs[, cell_top])
  cor_obs <- apply(obs, 2, function(x)cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x)cor(x, cnv_top))
  
  cnv <- data.frame(score=c(cnv_obs, cnv_ref), 
                    cor=c(cor_obs, cor_ref), 
                    barcode=c(colnames(obs), colnames(ref)))
  cnv$type <- 'Middle CNV'
  cnv$type[cnv$score>score_threshold & cnv$cor>cor_threshold] <- 'High CNV'
  cnv$type[cnv$score<score_threshold & cnv$cor<cor_threshold] <- 'Low CNV'
  return(cnv)
}

#我们计算一下
CNV_score_cor <- estimateCNV(obs, ref, 
                             score_threshold= 0.0018, 
                             cor_threshold= 0.20)
colnames(CNV_score_cor)[3] <- "cell_id"
#与之前的分组数据结合
library(tidyverse)
table(CNV_score_cor$type)
CNV_score_cor=CNV_score_cor%>%inner_join(kmeans_df,by="cell_id")
#因为我们需要判断EEC肿瘤中的Malignant cell，所以我们可以只提取EEC组的数据作图

CNV_score_cor_EEC <- subset(CNV_score_cor)

#作图
library(ggplot2)
ggplot(CNV_score_cor_EEC, aes(x=score,y=cor,color=type))+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title = element_text(color = 'black', size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(x="CNV Score",
       y="CNV Correlation")+
  scale_color_manual(values =c("#009E73","grey","black") )

ggsave('CNV_score_cor.pdf',width =10,height = 5 )

#当然了，网上还有其他的方法，也有文献有别的方式，只不过比较小众
#至于其他的，可自行查阅文献，主要的是自己对于结果的理解

# 将这些细胞的 CNV_type 更新为 "Normal"
data@meta.data$CNV <- " Middle CNV"

#Malignant
Tumor_cells <- CNV_score_cor_EEC[CNV_score_cor_EEC$type %in% c("High CNV"),]
Tumor_cells <- Tumor_cells[Tumor_cells$group %in% c("obs"),]

Tumor_cells <- which((data@meta.data$cell_id) %in% (Tumor_cells$cell_id))
data@meta.data$CNV[Tumor_cells] <- "High CNV"

#other
Other_cells <- CNV_score_cor_EEC[CNV_score_cor_EEC$type %in% c("Low CNV"),]
Other_cells <- which((data@meta.data$cell_id) %in% (Other_cells$cell_id))
data@meta.data$CNV[Other_cells] <- "Low CNV"

table(data@meta.data$CNV)

unique(data@meta.data$CNV)
save.image('PB_infercnv.RData')
saveRDS(data,file='PB_CNV.RDS')
