
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)

#===============================================================================
#                           基因随拟时表达变化
#===============================================================================
colnames(pData(mouse_monocle))#mouse_monocle是已经做好拟时的monocle2对象
#将需要展示的基因表达量添加到cds， 这里我使用的是log2处理，也可以使用log10进行标准化
genes <- c("Anxa1","Xist", "Ncf1", "Ltf")
genes_exp <- list()

for(i in 1:length(genes)){
  A <- log2(exprs(mouse_monocle)[genes[i],]+1)
  A <- as.data.frame(A)
  genes_exp[[i]] <- A
}

gene_exp <- do.call(cbind, genes_exp)
colnames(gene_exp) <- genes
#将上述几个基因的拟时表达添加到monocle
pData(mouse_monocle) = cbind(pData(mouse_monocle), gene_exp)


#提取作图数据，只需要游基因表达和拟时即可
data <- pData(mouse_monocle)
colnames(data)
#选择需要的列即可，我这里的origin.ident就是分组
data<-data %>% select("orig.ident","Pseudotime","Anxa1","Xist", "Ncf1", "Ltf")

#ggplot作图

#使用分屏，应该就是文献种的办法
#首先将data宽数据转化为长数据
data_long_m<-melt(data, id.vars = c("orig.ident", "Pseudotime"), #需保留的不参与聚合的变量列名
                  measure.vars = 3:6,#选择需要转化的列
                  variable.name = c('gene'),#聚合变量的新列名
                  value.name = 'value')#聚合值的新列名
colnames(data_long_m)

ggplot(data_long_m, aes(x=Pseudotime, y=value, color=orig.ident))+
  geom_smooth(aes(fill=orig.ident))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~variable, scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 14),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色



#===============================================================================
#                           通路随拟时表达变化
#===============================================================================
library(dplyr)
library(monocle)
library(ggplot2)
library(Seurat)
library(fgsea)
library(nichenetr)
mouse_data <- load(./mouse_data.rds)
#将拟时添加到seurat对象
mouse_sc <- mouse_data[, colnames(mouse_monocle@reducedDimS)]
mouse_sc$pseudotime <- mouse_monocle@phenoData@data$Pseudotime

#文章里面使用的是代谢通路，当让实际情况中并不一定是代谢通路
#所有通路基因的gmt文件可以去GSEA官网下载，这里我们以代谢通路为例子

pathway_file <- "./pathway.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
#选择需要的通路查看基因（自己关注的或者感兴趣的通路）
#这里我随便找两个举例子
Glutathione <- pathways[["Glutathione metabolism"]]
Caffeine <- pathways[["Caffeine metabolism"]]


#我用的是鼠的示例数据，所以需要将人的代谢基因转化为鼠的。如果是人的
#就不必操作了
Glutathione <- as.data.frame(Glutathione)
Glutathione$gene = Glutathione$Glutathione %>% convert_human_to_mouse_symbols()
#去除NA
Glutathione <- na.omit(Glutathione)
Glutathione_gene <- Glutathione$gene


Caffeine <- as.data.frame(Caffeine)
Caffeine$gene = Caffeine$Caffeine %>% convert_human_to_mouse_symbols()
#去除NA
Caffeine <- na.omit(Caffeine)
Caffeine_gene <- Caffeine$gene


#计算代谢通路评分，添加到seurat对象
DefaultAssay(mouse_sc) <- "RNA"
mouse_sc <- AddModuleScore(mouse_sc,
                           features=list('Glutathione' = Glutathione_gene,
                                         'Caffeine'=Caffeine_gene),
                           pool = rownames(mouse_sc), k=F, nbin=24,
                           name = c('Glutathione', 'Caffeine'))


#提取作图数据
data1 <- mouse_sc@meta.data
colnames(data1)
data1 <- data1[, c("orig.ident", "pseudotime", "Glutathione1","Caffeine2")]

data1_long   <-melt(data1, id.vars = c("orig.ident", "pseudotime"), #需保留的不参与聚合的变量列名
                    measure.vars = 3:4,#选择需要转化的列
                    variable.name = c('pathway'),#聚合变量的新列名
                    value.name = 'value')#聚合值的新列名


data1_new <- data1_long
levels(data1_new$pathway) <- c('Glutathione metabolism',
                               'Caffeine metabolism')


#作图和上述一样
ggplot(data1_new, aes(x=pseudotime, y=value, color=orig.ident))+
  geom_smooth(aes(fill=orig.ident))+ #平滑的填充
  xlab('pseudotime') + 
  ylab('Relative expression') +
  facet_wrap(~levels(pathway), scales = "free_y")+ #分面，y轴用各自数据
  theme(axis.text = element_text(color = 'black',size = 12),
        axis.title = element_text(color = 'black',size = 12),
        strip.text = element_text(color = 'black',size = 14))+ #分面标题
  scale_color_manual(name=NULL, values = c("#089E86","#3D5387"))+#修改颜色
  scale_fill_manual(name=NULL, values = c("#089E86","#3D5387"))#修改颜色





