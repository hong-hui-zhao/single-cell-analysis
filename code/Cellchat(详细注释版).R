library(devtools)
library(CellChat)
library(patchwork)

# 假设您有以下包含所有mutation标签的向量
mutation_tags <- c("G12D", "Non_G12D")  # 替换省略号为实际需要处理的所有mutation标签

# 初始化一个空列表来存储不同mutation标签对应的cellchat对象
cellchat_list <- list()

for (mutation_tag in mutation_tags) {
  
  # 根据当前mutation标签筛选数据
  current_subset <- subset(sce, mutation == mutation_tag)
  
  # 处理RNA数据和元数据
  current_datainput <- current_subset[["RNA"]]@data
  current_meta <- current_subset@meta.data
  current_meta <- current_meta[, c(27, 28)]
  colnames(current_meta) <- c("group", "labels")
  
  # 创建并配置CellChat对象
  current_cellchat <- createCellChat(object = current_datainput, 
                                     meta = current_meta, 
                                     group.by = "labels")
  
  current_cellchat <- addMeta(current_cellchat, meta = current_meta)
  current_cellchat <- setIdent(current_cellchat, ident.use = "labels")
  CellChatDB.use <- CellChatDB.human
  # 继续后续的CellChat分析步骤
  current_cellchat@DB <- CellChatDB.use
  current_cellchat <- subsetData(current_cellchat, features = NULL)
  current_cellchat <- identifyOverExpressedGenes(current_cellchat)
  current_cellchat <- identifyOverExpressedInteractions(current_cellchat)
  current_cellchat <- projectData(current_cellchat, PPI.human)
  current_cellchat <- computeCommunProb(current_cellchat, raw.use = TRUE)
  current_cellchat <- filterCommunication(current_cellchat, min.cells = 10)
  current_net <- subsetCommunication(current_cellchat)
  write.csv(current_net, file = paste0(mutation_tag, "_net_inter.csv"))
  
  current_cellchat <- computeCommunProbPathway(current_cellchat)
  current_cellchat <- aggregateNet(current_cellchat)
  current_cellchat <- netAnalysis_computeCentrality(current_cellchat)
  # 将当前mutation标签对应的cellchat对象添加到列表中
  cellchat_list[[mutation_tag]] <- current_cellchat
}
other_mutation <- cellchat_list[[2]]
G12D <- cellchat_list[[1]]
# 合并单独组的cellchat对象
object.list <- list(other_mutation = other_mutation,G12D= G12D)
merged_cellchat <- mergeCellChat(object.list, add.names = names(object.list))
#比较两组互作数目
gg1 <- compareInteractions(merged_cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(merged_cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("compareInteractions.pdf",height = 5,width = 8)
gg1 + gg2
dev.off()
devtools::install_github("jinworks/merged_cellchat")
library(merged_cellchat)
merged_cellchat <- updatemerged_cellchat(merged_cellchat)

pdf("netVisual_diffInteraction.pdf",height = 10,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(merged_cellchat, weight.scale = T)
netVisual_diffInteraction(merged_cellchat, weight.scale = T, measure = "weight")
dev.off()

gg1 <- netVisual_heatmap(merged_cellchat)
gg2 <- netVisual_heatmap(merged_cellchat, measure = "weight")
pdf("netVisual_heatmap.pdf",height = 5,width = 8)
gg1 + gg2
dev.off()


pdf("netVisual_diffInteraction.pdf",height = 10,width = 10)
weight.max <- getMaxWeight(cellchat_list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(cellchat_list)) {
  netVisual_circle(cellchat_list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(cellchat_list)[i]))
}
dev.off()


pdf("netVisual_bubble.pdf",height = 36,width = 36)
netVisual_bubble(merged_cellchat, comparison = c(1, 2), angle.x = 45)
dev.off()


gg1 <- rankNet(merged_cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(merged_cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 10)


pdf("netVisual_singnaling.pdf",height = 4,width = 8)
num.link <- sapply(cellchat_list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(cellchat_list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(cellchat_list[[i]], title = names(cellchat_list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)
dev.off()

pdf("out_singnaling.pdf",height = 10,width = 10)
library(ComplexHeatmap)
pathway.union <- union(cellchat_list[[1]]@netP$pathways, cellchat_list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cellchat_list[[1]], pattern = "outgoing", signaling = pathway.union, title = names(cellchat_list)[1], width = 5, height = 20)
ht2 = netAnalysis_signalingRole_heatmap(cellchat_list[[2]], pattern = "outgoing", signaling = pathway.union, title = names(cellchat_list)[2], width = 5, height = 20)
ht1+ht2
dev.off()

pdf("incoming_singnaling.pdf",height = 10,width = 10)
library(ComplexHeatmap)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 20, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 20, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()


pdf("overall_singnaling.pdf",height = 12,width = 10)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[1]], pattern = "all", signaling = pathway.union, title = names(object.list)[1], width = 5, height = 20, color.heatmap = "OrRd")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[2]], pattern = "all", signaling = pathway.union, title = names(object.list)[2], width = 5, height = 20, color.heatmap = "OrRd")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf("MHC_singnaling.pdf",height = 6,width = 8)
pathways.show <- c("MHC-II") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}
draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
dev.off()

save.image('cellchat.RData')
 