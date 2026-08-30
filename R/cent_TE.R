library(tidyverse)
library(DescTools)
# 0 preparation -----------------------------------------------------------

library(showtext)
#font_path <- "C:/Users/hanlab/AppData/Local/Microsoft/Windows/Fonts/" %>% font_paths()
#font_add('helvetica','helvetica.ttf')
showtext_auto()
# 如没有helvetica字体，可以考虑使用Califri，二者很相似


# windowsFonts(helvetica=windowsFont("Helvetica CE 55 Roman"),
#              Times=windowsFont("Times New Roman"),
#              Arial=windowsFont("Arial"))
# 将windows字体映射到ggplot中
# 

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(size = 12), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = margin(unit(0, "cm")),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



# 也可以使用下面两个函数
scale_fill_Publication <- function(color_name = 'Set2'){
  library(scales)
  library(dplyr)
  library(RColorBrewer)
  test_color <- as.data.frame(brewer.pal.info)
  test_color$name <- rownames(test_color)
  special_color <- filter(test_color,name == color_name)
  discrete_scale("fill","Publication",manual_pal(values = brewer.pal(special_color[,1],special_color[,4])))
  
}

scale_colour_Publication <- function(color_name = 'Set2'){
  library(scales)
  library(RColorBrewer)
  library(RColorBrewer)
  test_color <- as.data.frame(brewer.pal.info)
  test_color$name <- rownames(test_color)
  special_color <- filter(test_color,name == color_name)
  discrete_scale("colour","Publication",manual_pal(values = brewer.pal(special_color[,1],special_color[,4])))
  
}




# 1 cent TE  ----------------------------------------------------------------

TE <- read_delim('../data/CENH3/TE/HiTE.cent.family.bed',
                 delim = '\t',col_names = F) %>% 
  mutate(X4=str_extract(X4,"^[^/]+"),
         Sub = str_extract(X1,str_sub(X1,5,5)),
         TE = case_when(
           X4 %in% c("LTR", "SINE", "LINE") ~ "Retrotransposon",
           X4 == "RC" ~ "DNA",
           TRUE ~ X4  # 其他情况直接等于 X4
         ),
         len = X3-X2)
ggplot(TE) + 
  geom_bar(aes(y = X1,fill= TE), position = 'fill') +
  scale_fill_Publication('Set1') + 
  theme_Publication() +
  theme(
    legend.position = 'right',
    legend.direction = "vertical",
    # 调整分面标签大小（避免字体太大显得拥挤）
    strip.text = element_text(size = 10)
  ) +
  labs(x = 'Proportion', y = 'Chromosome', fill = 'TE Type')
ggsave('../figure/figure2/centromere.all_chr.bar.pdf',width = 6,height =3.6 )



cent.TE <- read_delim('../data/CENH3/WEW.cent.TE_type.bed',
                      delim = '\t',col_names = F) %>% 
  mutate(X4=str_extract(X4,"^[^/]+"),
         Sub = str_extract(X1,str_sub(X1,5,5)),
         TE = case_when(
           X4 %in% c("LTR", "SINE", "LINE") ~ "Retrotransposon",
           X4 == "RC" ~ "DNA",
           TRUE ~ X4  # 其他情况直接等于 X4
         ),
         len = X3-X2)


ggplot(cent.TE) + geom_bar(aes(x = Sub, fill = TE),position = 'fill') +
  scale_fill_Publication('Set1') + theme_Publication() + 
  theme(
    legend.position = 'right',
    legend.direction = "vertical"
  ) + labs(x = 'Subgenome',y = 'Proportion',fill = 'TE Type')

ggsave('../figure/Sf1/WEW.cent.TE.bar.pdf',width = 4.5,height = 4.5)
GTest(table(cent.TE$Sub, cent.TE$TE))
# G = 447.23, X-squared df = 4, p-value < 2.2e-16

data <- table(cent.TE$Sub, cent.TE$TE)
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v # 0.08843722


ggplot(cent.TE) +
  geom_bar(aes(x = Sub, y = len,fill = TE),
           position = 'fill',stat = 'identity') +
  scale_fill_Publication('Set1') + theme_Publication() + 
  theme(
    legend.position = 'right',
    legend.direction = "vertical"
  ) + labs(x = 'Subgenome',y = 'Proportion',fill = 'TE Type')
ggsave('../figure/Sf1/WEW.cent.TE.bar_len.pdf',width = 4.5,height = 4.5)

data_summary <- cent.TE %>%
  group_by(Sub, TE) %>%
  summarise(total_len = sum(len, na.rm = TRUE)) %>%
  ungroup()
data_table <- data_summary %>%
  pivot_wider(names_from = TE, values_from = total_len, values_fill = list(total_len = 0))

# 转换为矩阵格式（用于卡方检验）
data <- as.matrix(data_table[,-1])  # 去掉第一列（Sub）
rownames(data) <- data_table$Sub
GTest(data)
# G = 1623806, X-squared df = 4, p-value < 2.2e-16
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v # 0.1192097  

# 2. cent Retrotransposons ---------------------------------------------------

cent.TE <- read_delim('../data/CENH3/WEW.cent.TE_type.bed',
                      delim = '\t',col_names = F) %>% 
  filter(str_detect(X4,'LTR')) %>% 
  mutate(X4 = str_extract(X4, "(?<=/)[^/]+$"),
         Sub = str_extract(X1,str_sub(X1,5,5)),
         len = X3-X2)

cent.TE_prop <- cent.TE %>%
  group_by(Sub) %>%
  count(X4)  %>%              # 统计每个类别的频数
  mutate(proportion = n / sum(n),  # 计算比例
         label = paste0(round(proportion * 100, 1), "%")) 


ggplot(filter(cent.TE_prop,Sub == 'A'), aes(x = factor(1), y = proportion, fill = X4)) +
  geom_bar(stat = "identity", width = 1, color = "white") +           # 绘制条形图，添加白色边框
  coord_polar(theta = "y") +                                          # 转换为极坐标系
  geom_text(aes(label = label, x = 1.2, y = proportion), 
            size = 4, color = "black", position = position_stack(vjust = 0.5)) +  # 手动指定标签位置
  scale_fill_brewer(palette = "Set3") +                               # 使用更美观的调色板
  theme_void() +                                                      # 移除背景和坐标轴
  theme(
    legend.position = "right",                                        # 图例位置
    legend.direction = "vertical",                                    # 图例方向
    legend.text = element_text(size = 10),                            # 图例文字大小
    legend.title = element_text(size = 12, face = "bold"),            # 图例标题大小
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 图表标题
  ) +
  labs(
    fill = "TE Type",                                                 # 图例标题
    title = "Proportion of TE Types by count\n(A subgenome)"                      # 图表标题
  )
ggsave('../figure/Sf1/WEW.cent.Retro.bar.A.pdf',width = 5,height = 3.5)

ggplot(filter(cent.TE_prop,Sub == 'B'), aes(x = factor(1), y = proportion, fill = X4)) +
  geom_bar(stat = "identity", width = 1, color = "white") +           # 绘制条形图，添加白色边框
  coord_polar(theta = "y") +                                          # 转换为极坐标系
  geom_text(aes(label = label, x = 1.2, y = proportion), 
            size = 4, color = "black", position = position_stack(vjust = 0.5)) +  # 手动指定标签位置
  scale_fill_brewer(palette = "Set3") +                               # 使用更美观的调色板
  theme_void() +                                                      # 移除背景和坐标轴
  theme(
    legend.position = "right",                                        # 图例位置
    legend.direction = "vertical",                                    # 图例方向
    legend.text = element_text(size = 10),                            # 图例文字大小
    legend.title = element_text(size = 12, face = "bold"),            # 图例标题大小
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 图表标题
  ) +
  labs(
    fill = "TE Type",                                                 # 图例标题
    title = "Proportion of TE Types by count\n(B subgenome)"                      # 图表标题
  )
ggsave('../figure/Sf1/WEW.cent.Retro.bar.B.pdf',width = 5,height = 3.5)




data_summary <- cent.TE %>%
  group_by(Sub, X4) %>%
  count(X4) %>%
  ungroup()
data_table <- data_summary %>%
  pivot_wider(names_from = X4, values_from = n, values_fill = list(total_len = 0))
# 转换为矩阵格式（用于卡方检验）
data <- as.matrix(data_table[,-1])  # 去掉第一列（Sub）
rownames(data) <- data_table$Sub
GTest(data)
# G = 365.06, X-squared df = 3, p-value < 2.2e-16
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v # 0.08959512







cent.TE_prop <- cent.TE %>%
  group_by(Sub, X4) %>%
  summarise(total_len = sum(len, na.rm = TRUE)) %>% # 统计每个类别的频数
  mutate(proportion = total_len / sum(total_len),  # 计算比例
         label = paste0(round(proportion * 100, 1), "%")) 
ggplot(filter(cent.TE_prop,Sub =='A'), aes(x = factor(1), y = proportion, fill = X4)) +
  geom_bar(stat = "identity", width = 1, color = "white") +           # 绘制条形图，添加白色边框
  coord_polar(theta = "y") +                                          # 转换为极坐标系
  geom_text(aes(label = label, x = 1.2, y = proportion), 
            size = 4, color = "black", position = position_stack(vjust = 0.5)) +  # 手动指定标签位置
  scale_fill_brewer(palette = "Set3") +                               # 使用更美观的调色板
  theme_void() +                                                      # 移除背景和坐标轴
  theme(
    legend.position = "right",                                        # 图例位置
    legend.direction = "vertical",                                    # 图例方向
    legend.text = element_text(size = 10),                            # 图例文字大小
    legend.title = element_text(size = 12, face = "bold"),            # 图例标题大小
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 图表标题
  ) +
  labs(
    fill = "TE Type",                                                 # 图例标题
    title = "Proportion of TE Types by length\n(A subgenome)"                      # 图表标题
  )
ggsave('../figure/Sf1/WEW.cent.Retro.bar_len.A.pdf',width = 5,height = 3.5)

ggplot(filter(cent.TE_prop,Sub =='B'), aes(x = factor(1), y = proportion, fill = X4)) +
  geom_bar(stat = "identity", width = 1, color = "white") +           # 绘制条形图，添加白色边框
  coord_polar(theta = "y") +                                          # 转换为极坐标系
  geom_text(aes(label = label, x = 1.2, y = proportion), 
            size = 4, color = "black", position = position_stack(vjust = 0.5)) +  # 手动指定标签位置
  scale_fill_brewer(palette = "Set3") +                               # 使用更美观的调色板
  theme_void() +                                                      # 移除背景和坐标轴
  theme(
    legend.position = "right",                                        # 图例位置
    legend.direction = "vertical",                                    # 图例方向
    legend.text = element_text(size = 10),                            # 图例文字大小
    legend.title = element_text(size = 12, face = "bold"),            # 图例标题大小
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5)  # 图表标题
  ) +
  labs(
    fill = "TE Type",                                                 # 图例标题
    title = "Proportion of TE Types by length\n(B subgenome)"                      # 图表标题
  )
ggsave('../figure/Sf1/WEW.cent.Retro.bar_len.B.pdf',width = 5,height = 3.5)


data_summary <- cent.TE %>%
  group_by(Sub, X4) %>%
  summarise(total_len = sum(len, na.rm = TRUE)) %>%
  ungroup()
data_table <- data_summary %>%
  pivot_wider(names_from = X4, values_from = total_len, values_fill = list(total_len = 0))

# 转换为矩阵格式（用于卡方检验）
data <- as.matrix(data_table[,-1])  # 去掉第一列（Sub）
rownames(data) <- data_table$Sub
GTest(data)
# G = 652689, X-squared df = 3, p-value < 2.2e-16
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v # 0.0783784



# insertion time ----------------------------------------------------------
library(ggpubr)
cent.TE <- read_delim('../data/fl-LTR/WEW.cent.flTE.insert.bed',
                      delim = '\t',col_names = F) %>%
  mutate(Sub = str_extract(X1,str_sub(X1,5,5)),
         len = X3-X2)

wilcox.test(filter(cent.TE,Sub == 'A')$X6,filter(cent.TE,Sub == 'B')$X6)
# W = 4270878, p-value = 0.00272
ggplot(cent.TE,aes(x = Sub,y = X6,fill = Sub)) + 
  geom_violin(trim = FALSE, alpha = 0.5, width = 0.7) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.1) + 
  theme_Publication() + 
  scale_fill_Publication('Dark2') + 
  theme(
    legend.position = 'None'
  )  +
  stat_compare_means(
    method = "wilcox.test",          # 选择检验方法: "t.test", "wilcox.test"
    comparisons = list(c("A", "B")), # 比较组别
    label = "p.signif",              # 显示符号（***/<0.001, **/<0.01, */<0.05, ns/不显著）
    tip.length = 0.01,               # 调整标记线条长度
    size = 5,                        # 标注文字大小
    vjust = 0.5                      # 文字垂直位置
  ) + 
  labs(x = 'Subgenome', y = 'Insertion time(Mya)')
ggsave('../figure/figure1/cent.insertTime.subGenome.box.pdf',width = 3.3,height = 3.3)

library(cowplot)
main_plot <- ggplot(cent.TE, aes(x = X6, fill = Sub)) +
  geom_density(alpha = 0.4) +
  theme_Publication() +                      # 使用自定义主题（需提前定义或从ggthemes加载）
  scale_fill_Publication("Dark2") +          # 使用自定义颜色方案
  labs(x = "Insertion time (Mya)", y = "Density") +
  # 可选：在主图上标记放大区域（红色方框）
  annotate("rect", 
           xmin = 0, xmax = 0.1, 
           ymin = 9, ymax = 12.5,
           alpha = 0.3, color = "red", fill = "transparent")
inset_plot <- ggplot(cent.TE, aes(x = X6, fill = Sub)) +
  geom_density(alpha = 0.4) +
  coord_cartesian(                           # 放大坐标范围（不删除原始数据）
    xlim = c(0, 0.1), 
    ylim = c(9, 12.5)
  ) +
  theme_Publication() +
  scale_fill_Publication("Dark2") +
  # 简化子图样式：去除图例、调整字体、透明背景
  theme(
    legend.position = "none",                # 移除图例
    axis.title = element_blank(),            # 移除坐标轴标题
    axis.text = element_text(size = 6),      # 缩小坐标轴文字
    plot.background = element_rect(fill = "transparent")  # 透明背景
  ) + labs(x = '', y = '')

# ----------------------------------------
# 3. 组合主图和子图
# ----------------------------------------
final_plot <- ggdraw() +
  draw_plot(main_plot) +                     # 先画主图
  draw_plot(
    inset_plot, 
    x = 0.35, y = 0.48,                      # 调整子图位置（右上角）
    width = 0.5, height = 0.5             # 调整子图大小
  )

print(final_plot)
ggsave('../figure/figure1/cent.insertTime.subGenome.density.pdf',final_plot,
       width = 3.3,height = 3.3)

# A: 0.028Myr
# B: 0.035Myr˜

# 2. 提取A/B亚组的密度估计
Sub.A <- density(filter(cent.TE, Sub == "A")$X6, 
                 bw = "SJ",   # 推荐使用Sheather-Jones带宽
                 na.rm = TRUE) # 处理可能的NA

# 3. 检测峰值
tp <- turnpoints(Sub.A$y)  # 分析密度值序列

# 4. 提取重要峰值
peak_indices <- which(tp$peaks)            # 所有候选峰索引
peak_heights <- Sub.A$y[peak_indices]      # 对应峰高度
significant_peaks <- peak_indices[peak_heights > quantile(peak_heights, 0.75)]  # 取前25%高度的峰

# 获取峰值对应时间（Mya）
peak_times <- Sub.A$x[significant_peaks]

# 5. 可视化确认
ggplot(filter(cent.TE, Sub == "A"), aes(X6)) +
  geom_density(fill = "#1b9e77", alpha = 0.4) +  # Dark2调色板中A组的颜色
  geom_vline(xintercept = peak_times, color = "red", linetype = 2) +
  ggtitle("A Subgenome Insertion Time Peaks") +
  theme_Publication() +
  labs(x = "Insertion time (Mya)", y = "Density")

# 输出结果
cat("检测到", length(peak_times), "个显著峰:", "\n")
print(data.frame(Peak_Time_Mya = peak_times, 
                 Relative_Height = peak_heights[peak_heights > quantile(peak_heights, 0.75)]))

