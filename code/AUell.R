
library(Seurat)
BiocManager::install("UCell")
library(UCell)
DimPlot(uterus,label = T)+NoLegend()

#UCell评分函数AddModuleScore_UCell可以提供方seurat对象
#评分基因set是以list的形式提供
markers <- list()
markers$SM <- c("ACTA2", "RGS5")
markers$Marc <- c("MS4A6A", "CD68","LYZ")
markers$Ly <- c("CCL5", "STK17B","PTPRC")
markers$SF <- c("DCN", "COL6A3", "LUM")
markers$Endo <- c("PECAM1","PCDH17", "VWF")
markers$unEP <- c("EPCAM", "CDH1")
markers$cEP <- c("FOXJ1","CDHR3","DYDC2")

#评分计算
marker_score <- AddModuleScore_UCell(uterus,
                                     features=markers)
#可视化是Ucell score
library(stringr)
library(ggplot2)
library(viridis)
a <- colnames(marker_score@meta.data) %>% str_subset("_UCell")#提取包含_UCell的列名
FeaturePlot(marker_score,features = a,order = T, ncol = 4, cols = viridis(256))


#===================================================================================
#                                    应用1 通路评分
#===================================================================================
#假设我们要看代谢通路的评分，我们可以下辖代谢通路gmt文件
#如果是其他的KEGG通路，或者GO，Hallmarker等等，都可以在GSEA官网去下载对应的gmt

#这里的通路评分是按照通路基因集进行的，当然也可以自己选定基因进行评分
library(clusterProfiler)
metabolism <- read.gmt("KEGG_metabolism_nc.gmt") 
unique(metabolism$term)
#我们选择其中一条通路进行评分
Oxidative <- subset(metabolism, term=="Oxidative phosphorylation")
Oxidative <- list(Oxidative$gene)#将基因整成list
names(Oxidative)[1] <- 'Oxidative'
DefaultAssay(uterus) <- 'RNA'
metabolism_score <- AddModuleScore_UCell(sce,
                                         features=Oxidative,
                                         name="_metabolism_score")

#可视化所有细胞
FeaturePlot(metabolism_score,features = "Oxidative_metabolism_score",
            order = T,cols = viridis(256))

ggsave('Oxidative_metabolism_score.pdf',width = 8,height = 8)
FeaturePlot(metabolism_score,features = "Oxidative_metabolism_score",
            order = T,cols = viridis(256), split.by = 'main_cell_type')
ggsave('Oxidative_metabolism_score_cell.pdf',width = 30,height = 3)
#箱线图图展示每个细胞中的代谢评分（或者细胞亚群）
#或者展示不用分组的评分，只需提取相应的数据即可
library(ggrastr)
library(dplyr)
data<- FetchData(metabolism_score,vars = c("main_cell_type","Oxidative_metabolism_score"))
data$cellid <- case_when(data$celltype ==unique(data$celltype)[1] ~ "SMC",
                         data$celltype ==unique(data$celltype)[2] ~ 'Ly',
                         data$celltype ==unique(data$celltype)[3] ~ 'unEP',
                         data$celltype ==unique(data$celltype)[4] ~ 'SF',
                         data$celltype ==unique(data$celltype)[5] ~ 'cEP',
                         data$celltype ==unique(data$celltype)[6] ~ 'Endo',
                         data$celltype ==unique(data$celltype)[7] ~ 'Macro')
colors <- c('#507BA8','#F38D37','#5D9F53','#B5972D','#48998E','#E05758','#F1CE60')



ggplot(data, aes(x=cellid,y=Oxidative_metabolism_score,fill=cellid,color=cellid)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=12),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Oxidative_metabolism_score")+ 
  geom_jitter_rast(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colors)+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)


#添加细胞数--------------------------------------------------------------------------
cell_number <- as.data.frame(table(data$cellid))
cell_number$number <- paste0('n=',cell_number$Freq)
cell_number$number 
#获取每一组最大数值
max_score <- c()
for (i in cell_number$Var1){
  A = subset(data, cellid==i)
  a = max(A$Oxidative_metabolism_score)
  max_score <- append(max_score,a)
}

max_score
max_score <- max_score+0.05
cell_number$y <- max_score
colnames(cell_number)[1] <- 'cellid'


library(ggrepel)
ggplot(data, aes(x=cellid,y=Oxidative_metabolism_score,fill=cellid,color=cellid)) +
  theme_bw()+RotatedAxis()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(size=12),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none')+
  labs(x=NULL,y=NULL,title = "Oxidative_metabolism_score")+ 
  geom_jitter_rast(col="#00000033", pch=19,cex=2, position = position_jitter(0.2))+
  geom_boxplot(position=position_dodge(0))+
  scale_fill_manual(values = colors)+
  geom_boxplot(position=position_dodge(0),color='black',
               outlier.colour  = NA,outlier.fill=NA,outlier.shape=NA)+
  geom_text(data=cell_number,aes(cellid,y,label=number), color='black')


#===================================================================================
#                                    应用2 基因评分
#===================================================================================
#也可以做某一个基因集的评分，这个基因集可视化致病基因，可以是关注基因集合
#看看这些评分在我们样本之间的差异

#做基因集评分
genes <- list(c("PTEN","PIK3CA","KRAS","ARID1A","RCA1","WNT5A"))
names(genes) <- 'gene'
gene_score <- AddModuleScore_UCell(uterus,features=genes,name="_score")


#提取数据
library(ggpubr)
df<- FetchData(gene_score,vars = c("orig.ident","gene_score"))
df$orig.ident <- factor(df$orig.ident,levels = c("HC","EEC","AEH"))#设置顺序


#设置比较组
my_comparisons1 <- list(c("HC", "EEC")) 
my_comparisons2 <- list(c("AEH", "EEC"))
my_comparisons3 <- list(c("HC", "AEH"))

#做小提琴图
ggplot(df,aes(x=orig.ident,y=gene_score,fill=orig.ident))+
  geom_violin(color='black',size=1)+#小提琴
  theme_classic() + 
  theme(text = element_text(size=10, colour = "black")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        axis.text.x = element_text(colour = "black", size = 12),
        axis.text.y = element_text(colour = "black", size = 10),
        axis.title.y = element_text(color = 'black', size = 12),
        axis.line = element_line(size = 1))+ 
  theme(legend.position="none") +  
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +#箱线图
  stat_compare_means(method="t.test",hide.ns = F,
                     comparisons =c(my_comparisons1,my_comparisons2,my_comparisons3),
                     label="p.signif",
                     bracket.size=0.8,
                     size=6)#添加显著性比较

