setwd('~/KRAS/TCell/enrich')


library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(Seurat)

#------------------------------------------------------------------------------
#1、差异基因
# unique(sce$celltype)
celltypes = unique(sce$celltype)
#结合celltype与分组，做每种celltype两组之间差异基因
deg = c()
Idents(sce) = paste(sce$celltype, sce$mutation, sep = '_')

for (cell in celltypes)
  try({
    tmp = FindMarkers(sce, ident.1 = paste0(cell, '_G12D'), ident.2 = paste0(cell, '_Non_G12D'))
    tmp$gene = rownames(tmp)
    tmp$celltype = cell
    deg = rbind(deg, tmp)
  })

#显著性基因
deg = subset(deg, p_val_adj < 0.05)
deg$DE =''
deg$DE = ifelse(deg$avg_log2FC>0, 'up', 'down')#上下调分组


#------------------------------------------------------------------------------
#2、富集分析
group <- data.frame(gene=deg$gene,group=deg$celltype)#分组情况
#gene转化为ID
Gene_ID <- bitr(deg$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#构建文件并分析
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
diff_Pathway <- compareCluster(ENTREZID~group,
                          data=data,
                          fun = "enrichPathway",#函数选择什么定义什么分析
                          #ont = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff = 0.01,
                          qvalueCutoff = 0.01,
                          #OrgDb='org.Hs.eg.db'，
                         )#物种

#将gene ID转化为gene symbol
diff_Pathway = setReadable(diff_Pathway,OrgDb = "org.Hs.eg.db", keyType = "ENTREZID")
Pathway <- diff_Pathway@compareClusterResult


#------------------------------------------------------------------------------
#3、DE score

#df就是每一个富集terms里面，上调基因数-下调基因数， 除以总的count得到的score

celltypes = unique(sce$celltype)
Pathway.all = c()
for (cell in celltypes)
  try({
    
    #提取每种celltype富集分析terms基因
    Pathway.raw = subset(Pathway, Cluster==cell)
    Pathway.genes = as.character(Pathway.raw$geneID)
    Pathway.genes = str_split(Pathway.genes, '/')
    
    all.count = c()#表示富集terms的基因数Count
    up.count = c()#富集terms里面上调基因数
    down.count = c()#富集terms里面下调基因数
    de.gene = c()
    de.count = c()
    
    for (gene.list in Pathway.genes){
      tmp = length(gene.list)
      all.count = c(all.count, tmp)
      
      tmp = subset(deg, celltype == cell & DE == 'up' & gene %in% gene.list)
      up.gene = as.character(tmp$gene)
      up.count = c(up.count, length(up.gene))
      up.gene = paste(up.gene, collapse = '/')
      
      tmp = subset(deg, celltype == cell & DE == 'down' & gene %in% gene.list)
      down.gene = as.character(tmp$gene)
      down.count = c(down.count, length(down.gene))
      down.gene = paste(down.gene, collapse = '/')
      
      de.gene = rbind(de.gene, 
                      data.frame(Upregulated=up.gene,
                                 Downreguated=down.gene))
    }
    
    de.count = rbind(de.count, data.frame(upCount=up.count,
                                          downCount=down.count,
                                          allCount=all.count))
    
    Pathway.out = cbind(Pathway.raw, de.gene, de.count)
    colnames(Pathway.out)[1] = 'GroupID'
    Pathway.out$DE.score = (Pathway.out$upCount - Pathway.out$downCount) / (Pathway.out$downCount + Pathway.out$upCount)
    Pathway.out$celltype = cell
    
    Pathway.all = rbind(Pathway.all, Pathway.out)
    
    write.csv(Pathway.out, paste0(cell, '_out.csv'), row.names = F)
  })


write.csv(Pathway.all, file = 'Pathway.all.csv')

df <- Pathway.all %>%
  group_by(celltype) %>%
  top_n(3, wt = DE.score)
#------------------------------------------------------------------------------
#4、可视化
#可视化可以挑选自己需要的terms进行展示
#读入整理的数据
df <- read.csv("df.csv", header = T, row.names = 1)
df$Description <- factor(df$Description, levels = unique(df$Description))


#plot
p1 = ggplot(df,aes(x=-log10(qvalue),y=Description))+
  geom_col(aes(fill=group),width = 0.1)+ #连线
  geom_point(data=df,aes(size=allCount,color=DE.score))+#点
  theme_classic()+
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_line(colour = 'black', linewidth =0.5),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.ticks.x = element_line(colour = 'black'),
        axis.title.x = element_text(colour = 'black', size = 12))+
  scale_x_continuous(expand = c(0,0),limits = c(0, 12.5))+ #设置坐标轴，去除y轴与图形之间的间隙，修改X轴范围
  scale_fill_manual(values = c("#D51F26", "#00A08A", "#F98400", "#5BBCD6","#7FFFD4", "#98F5FF", "#9AFF9A", "#00C5CD", "#87CEFF"))+#修改颜色
  guides(fill=FALSE)+#不要fill的legend
  scale_color_gradientn(colors=c('#810f7c','#88419d','#8c6bb1','#f7fcfd',
                                 '#bfd3e6','#9ebcda','#023858'))#修改连续变量颜色




#注释

p2 = ggplot(df, aes(x=1,y=Description,fill=group))+
  geom_tile() + 
  theme_classic()+
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 10),
        legend.position = "none")+
  scale_fill_manual(values = c("#D51F26", "#00A08A", "#F98400", "#5BBCD6","#7FFFD4", "#98F5FF", "#9AFF9A", "#00C5CD", "#87CEFF"))+
  scale_x_continuous(expand = c(0,0))+
  geom_text(aes(1, 3, label = "Mast"), #添加标注文字，坐标轴的位置需要根据实际情况自行调整
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 7, label = "T and NK"),
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 10, label = "MARCO"),
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 13, label = "Plasma"),
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 3, label = "Monocyte"), 
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 3, label = "B Cell"), 
            size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 3, label = "Epithelial"), 
           size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 3, label = "Fibroblast"), 
           size = 4, fontface="italic",angle = 90)+
  geom_text(aes(1, 3, label = "Endothlial"), 
           size = 4, fontface="italic",angle = 90)
#拼图  
library(deeptime)
pdf("2_wrap_plots_hMEs.pdf",width=100,height=80)
ggarrange2(p2, p1, nrow = 1)
dev.off()



