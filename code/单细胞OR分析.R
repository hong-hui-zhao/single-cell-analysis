
library("sscVis")
library("data.table")
library("grid")
library("cowplot")
library("ggrepel")
library("readr")
library("plyr")
library("ggpubr")
library("ggplot2")
library("dplyr")
library("tidyr")

setwd("D:/KS项目/公众号文章/张泽明code学习")
#加载函数及数据
mouse_data <- readRDS("./mouse_data.rds")
source("./test_function.R")
source("./draw_analysis.R")

#数据分析
A <- do.tissueDist(cellInfo.tb = mouse_data@meta.data,#这里的输入文件需要的是整理好的有分组和细胞类型的metadata文件
              out.prefix = "./ORplot",#设定导出的文件名，自己设置
              pdf.width = 5,#热图的宽设置
              pdf.height = 8,#热图的高度设置
              verbose=1,#设置为1文件以list形式存储
              meta.cluster = 'celltype',#这里是细胞类型，也可以是seurat_clusters，名称没有要求，就是你细胞类型的列名
              loc = 'SingleCellNet_Xie',#这里就是分组，metadata中分组的列名，至于命名没有要求
              z.hi=15)#热图legend最大值，根据实际情况自己设置

#查看并保存文件
A$OR.dist.mtx #做热图数据，OR值
A$p.dist.tb #p值
A$OR.dist.tb #OR值
A$count.dist.melt.ext.tb#组合表，adjust-pvalue等


###################################################################
#自己做图
data <- A$count.dist.melt.ext.tb
write.csv(data, file = 'data.csv')
library(ggplot2)
library(RColorBrewer)
data <- read.csv("data.csv", header = T)

ggplot(data, aes(cid, rid)) + 
  geom_tile(aes(fill = OR), colour = "black", size = 0.6)+
  scale_fill_gradientn(name='OR',
                       colours=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+
  theme_minimal() + 
  theme(axis.title.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(), 
        axis.text.y = element_text(size = 10,color = 'black')) + 
  scale_y_discrete(position = "right")
  
  
  
  