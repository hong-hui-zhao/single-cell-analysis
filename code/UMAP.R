library(tidydr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cols4all)
setwd('~/KRAS/Epi/basis')
UMAP <- as.data.frame(sce@reductions$umap@cell.embeddings)
mutation <- sce@meta.data$mutation
UMAP <- cbind(UMAP,mutation)
head(UMAP)

#建立自定义主题：
mytheme <- theme_void() + #空白主题，便于我们后期添加tSNE箭头
  theme(plot.margin = margin(5.5,15,5.5,5.5)) #画布空白页缘调整

#建立映射，添加散点：
p <- ggplot(data = UMAP, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(aes(color = mutation),
             size = 0.4,
             alpha = 0.8)
p

#添加实线椭圆置信区间：
p1 <- p +
  stat_ellipse(aes(color = mutation),
               level = 0.95, linetype = 1, show.legend = F) +
  mytheme
p1

#添加虚线椭圆置信区间：
p2 <- p +
  stat_ellipse(aes(color = mutation),
               level = 0.95, linetype = 2, show.legend = F) +
  mytheme
p2

#添加填充型置信区间：
p3 <- p +
  stat_ellipse(aes(color = mutation, fill = mutation),
               level = 0.95, linetype = 1, show.legend = F,
               geom = 'polygon', alpha = 0.1) +
  mytheme
p3

#添加tsne坐标轴箭头：
p4 <- p3 +
  theme_dr(xlength = 0.2, #x轴长度
           ylength = 0.2, #y轴长度
           arrow = grid::arrow(length = unit(0.1, "inches"), #箭头大小/长度
                               ends = 'last', type = "closed")) + #箭头描述信息
  theme(panel.grid = element_blank())
p4

#在图中增加亚群标签:
##计算每个亚群散点的中位数，作为图中标签的坐标
label <- UMAP %>%
  group_by(mutation)%>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))
head(label)

p5 <- p4 +
  geom_text(data = label,
            aes(x = UMAP_1, y = UMAP_2, label = mutation),
            fontface = "bold", #粗体强调
            color = 'black', size = 4)
p5
p6 <- p4 +
  guides(color = guide_legend(override.aes = list(size = 5))) #放大图例中的散点

pdf('UMAP_mutation_Epi.pdf',width = 8,height = 8)
p6
dev.off()
