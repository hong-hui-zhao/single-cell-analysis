# 读取数据
data <- read.csv("cellper.csv")

# 将数据乘以100
data[, 2:17] <- data[, 2:17] * 100

# 设置颜色
colors <- c("#1f78b4", "#33a02c")  # Normal: blue, Tumor: green

# 创建一个文件夹保存生成的箱线图
dir.create("boxplots", showWarnings = FALSE)

# 创建一个空的列表来存储箱线图
plots_list <- list()

# 遍历每个细胞类型
for (cell_type in colnames(data)[2:17]) {
  # 绘制箱线图
  p <- ggplot(data, aes(x = group, y = !!as.name(cell_type), fill = group)) +
    geom_boxplot() +
    scale_fill_manual(values = colors) +
    labs(title = paste("Boxplot of", cell_type, "by SampleType"),
         fill = "SampleType") +
    theme_minimal() +
    ylab(paste(cell_type, "(%)"))  # 添加y轴标签，表示百分比
  

  
  # 将绘制好的图加入列表
  plots_list[[cell_type]] <- p
}

# 合并箱线图
merged_plot <- plot_grid(plotlist = plots_list, nrow = 4, ncol = 4)

# 保存合并的箱线图为PDF
ggsave("boxplots/merged_boxplots_with_pvalues.pdf", merged_plot, width = 40, height = 25, units = "cm")

# 打印消息
print("Merged boxplots with p-values saved as 'merged_boxplots_with_pvalues.pdf' in the 'boxplots' folder.")
