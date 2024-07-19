
# setwd('D:/KS项目/公众号文章/KEGG通路注释')
library(dplyr)
library(ggplot2)
library(forcats)
library(dittoSeq)
kegg <- read.csv('KEGG.csv', header = T)
#其实做一个柱状图没有什么难度
ggplot(kegg,aes(Description, Gene.Count))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=Gene.Count, y=Gene.Count+2),size=3)+
  coord_flip()+
  labs(x='',y='Gene count')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black', size = 10))

#注释文件
kegg_ann <- read.csv('kegg_annotation.csv')
rownames(kegg_ann) <- kegg_ann$directory3
df <- kegg_ann[kegg$Description,]
kegg1 <- cbind(kegg, df)
table(kegg1$directory1)


# Cellular Processes Environmental Information Processing 
# 1                                    3 
# Human Diseases                           Metabolism 
# 3                                    3 
# Organismal Systems 
# 2 
kegg1 <- kegg1[order(kegg1$directory1),]
kegg1$Description <- as.factor(kegg1$Description)
kegg1$Description <- fct_reorder(kegg1$Description)

p <- ggplot(kegg1,aes(Gene.Count,Description))+
  geom_bar(stat = "identity", aes(fill=directory1))+
  geom_text(aes(label=Gene.Count, x=Gene.Count+2),size=3)+
  labs(y='',x='Gene count')+
  theme_classic()+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = 'black', size = 10),
        plot.margin = margin(0,0,0,-0.05, "cm"))+
  scale_fill_manual(values = dittoColors())



#构建注释
df1<-data.frame(x="A", y=kegg1$Description, group=kegg1$directory1)

df1$group <- factor(df1$group, levels = c("Organismal Systems",
                                          "Metabolism",
                                          "Human Diseases",
                                          "Environmental Information Processing",
                                          "Cellular Processes"))

p1 <- ggplot(df1,aes(x,y, fill=group))+
  geom_tile(show.legend = F)+
  facet_grid(group~.,scales = 'free',space = 'free')+
  labs(x=NULL,y=NULL)+
  scale_x_discrete(expand = c(0,0))+
  scale_y_discrete(expand = c(0,0))+
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        plot.margin = margin(0,-0.05,0,0, "cm"),
        strip.text.y = element_text(angle=0,size=10,color = "black",
                                    hjust = 0,margin = margin(b = 3,t=3)),
        strip.background = element_rect(colour=NULL,fill = 'white'),
        panel.spacing=unit(0, "mm"))+
  scale_fill_manual(values = c("#0072B2","#F0E442","#009E73","#56B4E9","#E69F00"))
  
 
#拼图
library(deeptime)
ggarrange2(p, p1,nrow = 1,widths =c(2,0.05))

