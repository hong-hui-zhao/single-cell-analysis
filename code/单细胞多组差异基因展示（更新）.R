
setwd('~/KRAS/Basis')
library(msigdbr)
library(org.Hs.eg.db)
# 提取C6库
library(msigdbr)
library(clusterProfiler)
library(tidyverse)
library(Seurat)
library(scRNAseq)
library(GseaVis)
library(ReactomePA)
# 获取MSigDB hallmark基因集信息
df <- msigdbr(species = "Homo sapiens")
df <- df %>% dplyr::select(gs_name, gene_symbol)
# 获取所有唯一的细胞簇标识
cluster <- as.character(unique(Idents(data)))
celltype <- as.character(unique(Idents(data)))
# 初始化存储差异表达基因和富集结果的列表
DEGs <- list()
enrichment_results <- list()
for (i in 1:length(cluster)) {
  
  # 提取特定细胞亚群的数据
  sce <- subset(data, celltype == cluster[i])
  
  # 执行差异表达分析
  dges <- FindMarkers(sce, min.pct = 0.25, 
                      logfc.threshold = 0.25,
                      group.by = "mutation",
                      ident.1 ="G12D",
                      ident.2="Non_G12D")
  
  # 转换为tibble格式并添加细胞类型信息
  dges_tbl <- as_tibble(cbind(gene = rownames(dges), dges))
  dges_tbl$celltype <- cluster[i]
  
  # 将差异表达基因列表存储到DEGs中
  DEGs[[i]] <- dges_tbl
  names(DEGs)[i] <- paste0(cluster[i],"_DEGs")
  
  # 根据差异表达基因计算得分（这里假设使用log fold change作为得分依据）
  geneList <- dges_tbl$avg_log2FC
  names(geneList) = dges_tbl$gene
  geneList = sort(geneList, decreasing = TRUE)
  GSEA.res <- GSEA(geneList,minGSSize = 10,maxGSSize = 500,
                   pvalueCutoff = 0.15,pAdjustMethod = "BH",
                   verbose = FALSE,eps = 0,TERM2GENE = df)
  
  # 存储富集结果
  GSEA.result <- GSEA.res@result
  enrichment_results[[i]] <- GSEA.result
  
  # 绘制GSEA富集条形图并保存为PDF文件
  pdf(paste0(cluster[i], "_enrichment_gsea_barplot.pdf"),width = 16 ,height = 8)
  dotplotGsea(data = GSEA.result,topn = 10,order.by = 'NES',add.seg = T)
  dev.off()
}

for (i in 1:length(cluster)){
  pdf(paste0(cluster[9], "_enrichment_gsea_barplot.pdf"),width = 16 ,height = 8)
  dotplotGsea(data = enrichment_results[[9]],topn = 10,order.by = 'NES',add.seg = T)
  dev.off()
}
# 现在DEGs和enrichment_results列表分别包含了每个细胞亚群的差异表达基因和GSEA分析结果


diff <- do.call(rbind, DEGs)
#显著性
diff <- diff[which(abs(diff$avg_log2FC)>=0.4),]
diff$label <- ifelse(diff$p_val_adj<=0.05,"sig","unsig")

##label gene
top_gene = diff %>% group_by(celltype) %>% top_n(n = 10, wt = abs(avg_log2FC))

###bakground
max_fc <- c()
min_fc <- c()
for (i in 1:length(celltype)) {
  
  df <- subset(diff, celltype == cluster[i])
  df_max <- max(df$avg_log2FC)
  df_min <- min(df$avg_log2FC)
  max_fc <- append(max_fc, df_max)
  min_fc <- append(min_fc, df_min)
}

bakground_FC <- data.frame(x=cluster,
                           max_fc = max_fc+0.2,
                           min_fc = min_fc-0.2)

###cell label

cluster_label <-data.frame(x=cluster,
                           y=0,
                           label=c("SMC","Ly","UEC","SF","CEP","EC","MAC"))

#作图

diff$celltype <- factor(diff$celltype, levels = cluster)
bakground_FC$x <- factor(bakground_FC$x, levels = cluster)

ggplot()+
  geom_col(data = bakground_FC,
           mapping = aes(x = x,y = max_fc),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_col(data = bakground_FC,
           mapping = aes(x = x,y = min_fc),
           fill = "#dcdcdc",alpha = 0.6)+
  geom_jitter(data = diff,
              aes(x = celltype, y = avg_log2FC, color = label),
              size = 1,
              width =0.4)+
  scale_color_manual(values = c("#E64B357F","#00A0877F"))+
  geom_tile(data = cluster_label,
            aes(x=x,y=y),
            height=0.8,
            color = "black",
            alpha = 0.6,
            show.legend = F,
            fill=dittoColors()[1:7])+
  geom_text(data=cluster_label,
            aes(x=x,y=y,label=label),
            size =3,
            color ="white")+
  geom_text_repel(data=top_gene,
                  aes(x=celltype,y=avg_log2FC,label=gene),
                  color="black", size=2.5, fontface="bold.italic", 
                  arrow = arrow(ends="first", length = unit(0.01, "npc")),
                  box.padding = 0.2,
                  point.padding = 0.3, 
                  segment.color = 'black', 
                  segment.size = 0.3, force = 1, 
                  max.iter = 3e3)+
  theme_minimal()+
  theme(axis.title = element_text(size = 12,
                                  color = "black",
                                  face = "bold"),
        axis.line.y = element_line(color = "black",
                                   size = 0.8),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top",
        legend.direction = "vertical",
        legend.justification = c(1,0),
        legend.text = element_text(size = 12))