library(tidyverse)
library(ggsignif)
library(reshape2)
library(RColorBrewer)
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




# cenh3 -------------------------------------------------------------------

cenh3.1 <- read_delim('../data/CENH3/SRR23029304.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type = 'Rep1') 

cenh3.2 <- read_delim('../data/CENH3/SRR23029306.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type = 'Rep2')
cenh3.3 <- read_delim('../data/CENH3/zavitan-CENH3.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type = 'Rep3')
cenh3 <- rbind(cenh3.1,cenh3.2,cenh3.3)%>%
  mutate(
    X4 = case_when(
      X1 == "chr3A" & X2 == 1.79e+08 & Type == 'Rep1' ~ 0.06140510,
      X1 == "chr4B" & X2 == 1.10e+07 & Type == 'Rep1' ~ 0.06140510,
      X1 == "chr5A" & X2 == 2.92e+08 & Type == 'Rep1' ~ 0.1051305,
      X1 == "chr7B" & X2 == 3.53e+08 & Type == 'Rep2' ~ 0.09321214,
      TRUE ~ X4  # 其他情况保留原值
    )
  )

ggplot(cenh3, aes(x = (X2 + X3) / 2/1e6, y = X4,color = Type)) + ylim(0,2.2) +
  geom_line() +  # 绘制折线图
  facet_wrap(~ X1,ncol = 2) +  # 根据 X1 分面
  theme_Publication() + theme(
    axis.text.x = element_text(angle = 60,hjust = 1,vjust = 1),
    strip.text.x = element_text(size = 10),
    strip.background = element_blank(),
    axis.text.y = element_text(),
    legend.position = 'right',
    legend.direction = 'vertical'
  ) + 
  labs(x = "Position (Mb)",y = 'Relative Signal') +
  scale_colour_Publication()
ggsave('../figure/figure2/CENH3.densiyt.pdf',height = 7,width = 5.5)


# CENH3 TE type -----------------------------------------------------------
cent_TE <- read_delim('../data/CENH3/WEW_new.cent.TE.proportion',delim = '\t',
                      col_names = c('Chr','RT','DNA','simple/low-complexity repeats','Others')) %>% 
  melt(id = 'Chr')
  
cent_len <- read_delim('../data/CENH3/WEW_new.cent.bed',delim = '\t',
                       col_names = c('Chr','Start','End')) %>% mutate(
                         Len = (End - Start)/1e6,
                         Sub = str_sub(Chr,5,5)
                       )
cent_len$ratio <- signif(c(
  227500000/(611773303-227500000),
  262750000/(735254030-262750000),
  380250000/(790741867-380250000),
  378185000/(842231976-378185000),
  328200000/(766661948-328200000),
  361200000/(872156323-361200000),
  313350000/(751530861-313350000),
  289000000/(696107520-289000000),
  247500000/(722676616-247500000),
  221850000/(738030801-221850000),
  295000000/(634596046-295000000),
  328350000/(759895744-328350000),
  364350000/(752998143-364350000),
  321250000/(773842621-321250000)
),2)
ggplot(cent_len) + geom_bar(aes(x = Len, y =Chr, fill = Sub),stat = 'identity') + 
  geom_text(
    aes(x = Len, y = Chr, label = ratio),  # 使用 Len 列作为标签
    position = position_stack(vjust = 1.05),  # 堆叠柱状图的中心位置
    size = 3,  # 字体大小
    color = "black"  # 字体颜色
  ) + 
  theme_Publication() + scale_fill_Publication('Paired') 
ggsave('../figure/figure2/cent_len.bar.pdf',width = 5.5,height = 4)

cent_TE$variable <- factor(cent_TE$variable,
                           levels = c('DNA','RT','simple/low-complexity repeats','Others'))
ggplot(cent_TE) + geom_bar(aes(y = Chr, x = value, fill = variable),
                           stat = 'identity',position = 'fill') + 
  theme_Publication() + scale_fill_Publication('Dark2')
ggsave('../figure/figure2/cent_TE_types.bar.pdf',width = 7.2,height = 4.5)




# metaplot chromosome -----------------------------------------------------
cenh_slop <- read_delim('../data/CENH3/all_CENH3.100k.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type ='CENH3')
Cereba_slop <- read_delim('../data/CENH3/Cereba_WEW_new.100k.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type = 'Cereba')
Quinta_slop <- read_delim('../data/CENH3/Quinta_WEW_new.100k.bg',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% mutate(Type = 'Quinta')

slop_signal <- rbind(cenh_slop,Cereba_slop,Quinta_slop)
# 针对chr1A计算缩放因子
chr1A_data <- filter(slop_signal, X1 == 'chr1A')
scale_factor <- max(chr1A_data$X4[chr1A_data$Type == "CENH3"], na.rm = TRUE) / 
  max(chr1A_data$X4[chr1A_data$Type != "CENH3"], na.rm = TRUE)

# 创建新列用于缩放后的值
chr1A_data <- chr1A_data %>%
  mutate(scaled_value = ifelse(Type == "CENH3", X4, X4 * scale_factor))

# 设置颜色映射
line_colors <- c("CENH3" = "#D95F02", "Cereba" = "#1B9E77", "Quinta" = "#7570B3")

# 绘制双坐标轴图
ggplot(chr1A_data) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type == "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type != "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  scale_y_continuous(
    name = "CENH3 Signal",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cereba/Quinta Signal")
  ) +
  scale_color_manual(values = line_colors) +
  xlim(160, 280) + 
  labs(x = "Position (Mb)") +
  theme_Publication() + 
  scale_colour_Publication() +
  theme(
    axis.title.y.left = element_text(color = line_colors["CENH3"]),
    axis.title.y.right = element_text(color = "grey35"),
    axis.text.y.left = element_text(color = line_colors["CENH3"]),
    axis.text.y.right = element_text(color = "grey35")
  )
ggsave('../figure/figure2/Cereba_CENH3.chr1A.metaplot.pdf',width =5,height = 3)

chr2A_data <- filter(slop_signal, X1 == 'chr2A')
scale_factor <- max(chr2A_data$X4[chr2A_data$Type == "CENH3"], na.rm = TRUE) / 
  max(chr2A_data$X4[chr2A_data$Type != "CENH3"], na.rm = TRUE)
chr2A_data <- chr2A_data %>%
  mutate(scaled_value = ifelse(Type == "CENH3", X4, X4 * scale_factor))
ggplot(chr2A_data) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type == "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type != "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  scale_y_continuous(
    name = "CENH3 Signal",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cereba/Quinta Signal")
  ) +
  scale_color_manual(values = line_colors) +
  xlim(350, 420) + 
  labs(x = "Position (Mb)") +
  theme_Publication() + 
  scale_colour_Publication() +
  theme(
    axis.title.y.left = element_text(color = line_colors["CENH3"]),
    axis.title.y.right = element_text(color = "grey35"),
    axis.text.y.left = element_text(color = line_colors["CENH3"]),
    axis.text.y.right = element_text(color = "grey35")
  )
ggsave('../figure/figure2/Cereba_CENH3.chr2A.metaplot.pdf',width =5,height = 3)


chr3A_data <- filter(slop_signal, X1 == 'chr3A')
scale_factor <- max(chr3A_data$X4[chr3A_data$Type == "CENH3"], na.rm = TRUE) / 
  max(chr3A_data$X4[chr3A_data$Type != "CENH3"], na.rm = TRUE)
chr3A_data <- chr3A_data %>%
  mutate(scaled_value = ifelse(Type == "CENH3", X4, X4 * scale_factor))
ggplot(chr3A_data) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type == "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  geom_line(
    aes(x = (X2 + X3) / 2 / 1e6, y = ifelse(Type != "CENH3", scaled_value, NA), color = Type),
    stat = 'identity'
  ) +
  scale_y_continuous(
    name = "CENH3 Signal",
    sec.axis = sec_axis(~ . / scale_factor, name = "Cereba/Quinta Signal")
  ) +
  scale_color_manual(values = line_colors) +
  xlim(275, 375) + 
  labs(x = "Position (Mb)") +
  theme_Publication() + 
  scale_colour_Publication() +
  theme(
    axis.title.y.left = element_text(color = line_colors["CENH3"]),
    axis.title.y.right = element_text(color = "grey35"),
    axis.text.y.left = element_text(color = line_colors["CENH3"]),
    axis.text.y.right = element_text(color = "grey35")
  )
ggsave('../figure/figure2/Cereba_CENH3.chr3A.metaplot.pdf',width =5,height = 3)

# cent Insertime ----------------------------------------------------------
cent_LTR <- read_delim('../data/fl-LTR/LTR.cent.insertTime.list',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Obs',Time = (1-X1/100)/(2*1.3e-8)/1e6)
cent_LTR_null <- read_delim('../data/fl-LTR/LTR.cent.null.insertTime.list',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Exp',Time = (1-X1/100)/(2*1.3e-8)/1e6)

cent_sig <- rbind(cent_LTR,cent_LTR_null)

ggplot(cent_sig,aes(y = Time,  x = Type)) +
  geom_violin(aes(color = Type,),width = 1.2) + 
  geom_boxplot(aes(color = Type),outlier.size = 0.1,width = 0.05) + 
  theme_Publication() + scale_fill_Publication('Dark2') + scale_colour_Publication('Dark2') + 
  geom_signif(comparisons = list(c("Exp", "Obs")), y_position = 6,
              map_signif_level=TRUE) + ylim(c(0,6.5)) + 
  labs(y = 'Time(Mya)')
ggsave('../figure/figure2/cent_Insertime.violin.pdf',width = 2.4 ,height = 3.8)




cent_LTR <- read_delim('../data/fl-LTR/fl_LTR.cent.bed',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% 
  mutate(Time = (1-X5/100)/(2*1.3e-8)/1e6,Type = 'Cent',Sub = str_sub(X1,5,5))
arm_LTR <- read_delim('../data/fl-LTR/fl_LTR.non_cent.bed',delim = '\t',col_names = F) %>% 
  filter(str_detect(X1,'chr')) %>% 
  mutate(Time = (1-X5/100)/(2*1.3e-8)/1e6,Type = 'Arm',Sub = str_sub(X1,5,5))


LTR_dis <- rbind(cent_LTR,arm_LTR)
ggplot(cent_LTR,aes(x = Time)) + geom_density(aes(color = Sub)) + 
  scale_colour_Publication('Dark2') + theme_Publication()
ggsave('../figure/figure2/AB_cent.LTR.insertTime.density.pdf',width = 3,height = 2.2)
# median(cent_LTR$Time)
# [1] 0.4653846
ggplot(arm_LTR,aes(x = Time)) + geom_density(aes(color = Sub)) + 
  scale_colour_Publication() + theme_Publication()
ggsave('../figure/figure2/AB_arm.LTR.insertTime.density.pdf',width = 3,height = 2.2)
# median(arm_LTR$Time)
# [1] 1.869231

LTR_dis$X4 <- as.numeric(LTR_dis$X4)
library(ggpubr) 
ggplot(LTR_dis,aes(x = Type,y = Time,group = Type,color = Type)) + geom_boxplot(outlier.size = 0.1) + 
  scale_colour_Publication('Paired') + theme_Publication() + ylim(0,6.6) + 
  geom_signif(comparisons = list(c("Arm", "Cent")), y_position = 6,map_signif_level=TRUE)
ggsave('../figure/figure2/arm_cent.insertTime.box.pdf',height = 3,width = 2.6)

# cent_LTR_AB -------------------------------------------------------------

A_cent_LTR <- read_delim('../data/fl-LTR/LTR.cent_A.insertTime.list',delim = '\t',col_names = F) %>% 
  mutate(Type = 'A',Time = (1-X1/100)/(2*1.3e-8)/1e6)
B_cent_LTR <- read_delim('../data/fl-LTR/LTR.cent_B.insertTime.list',delim = '\t',col_names = F) %>% 
  mutate(Type = 'B',Time = (1-X1/100)/(2*1.3e-8)/1e6)
AB_cent_LTR <- rbind(A_cent_LTR,B_cent_LTR)

ggplot(AB_cent_LTR,aes(y = Time,  x = Type)) +
  geom_violin(aes(color = Type,),width = 0.3) + 
  geom_boxplot(aes(color = Type),outlier.size = 0.1,width = 0.05) + 
  theme_Publication() + scale_fill_Publication('Dark2') + scale_colour_Publication('Dark2') + 
  geom_signif(comparisons = list(c("A", "B")), y_position = 6,
              map_signif_level=TRUE) + ylim(c(0,6.5)) + 
  labs(y = 'Time(Mya)')
ggsave('../figure/Sf3/AB_insertTime.box.pdf',width =2.6 ,height = 3)
median(A_cent_LTR$Time)
# [1] 0.04692308
median(B_cent_LTR$Time)
# [1] 0.04615385

median(cent_LTR$Time)
# [1] 0.04653846

non_cent_LTR <- read_delim('../data/fl-LTR/LTR.non_cent.insertTime.list',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Non-cen',Time = (1-X1/100)/(2*1.3e-8)/1e6)
median(non_cent_LTR$Time)
# [1] 0.1853846
# cenh3_Te ----------------------------------------------------------------

TR_LTR_time <- read_delim('../data/fl-LTR/TR_flLTR.insertionTime',delim = '\t',col_names = F) %>% 
  mutate(Type = 'TR')
peri_LTR_time <- read_delim('../data/fl-LTR/peri40M_flLTR.insertionTime',delim = '\t',col_names = F) %>% 
  mutate(Type = 'Peri')

tmp <- rbind(TR_LTR_time,peri_LTR_time)
ggplot(tmp) + stat_ecdf(aes(x = X1,group = Type,color = Type)) + 
  scale_colour_Publication('Set1') + theme_Publication()+ 
  xlim(c(0,4e6)) + 
  theme(
    legend.position = 'right',
    legend.direction = 'vertical'
  ) + 
  labs(x = 'Insertion time(year)',y = 'Cumulative Proportion')
ggsave('../figure/figure2/TR_insertionTime.ecdf.pdf',width = 4,height = 3)

ggplot(tmp,aes(x = Type, y = X1)) + geom_boxplot(aes(x = Type,y = X1,color = Type)) + 
  scale_colour_Publication('Set1') + theme_Publication()+ 
  ylim(c(0,4e6)) + 
  theme(
    legend.position = 'None'
  ) + labs(x = 'Type', y = 'Insertion time(year)') + 
  geom_signif(comparisons=list(c("Peri", "TR")),map_signif_level = T,size = 0.9,y_position = 3.6e6)
  
ggsave('../figure/figure2/TR_insertionTime.boxplot.pdf',width = 2.5,height = 3) 



# tree --------------------------------------------------------------------

library(ape)      # 基础树操作
library(treeio)   # 增强的树文件读取
library(ggtree)   # 树的可视化
library(tidytree)

flLTR_tree <- read.newick("../data/fl-LTR/fl_LTR.cent_nonCent.sample3k.aligned.fa2.nwk")
flLTR_tree <- as.treedata(flLTR_tree)
tip_labels <- flLTR_tree@phylo$tip.label

# 创建注释数据框
annot_df <- data.frame(
  label = tip_labels,
  # 第一层注释：根据cent/arm分组
  Layer1 = ifelse(grepl("cent", tip_labels), "cent", "arm"),
  # 第二层注释：根据A_/B_分组
  Layer2 = ifelse(grepl("A_", tip_labels), "A", "B"),
  # 第三层注释：提取染色体号
  Layer3 = gsub(".*(chr[1-7]).*", "\\1", tip_labels)
)

# 设置颜色方案
colors_layer1 <- brewer.pal(3, "Set1")[1:2]  # cent/arm 颜色
names(colors_layer1) <- c("cent", "arm")

colors_layer2 <- brewer.pal(3, "Set2")[1:2]  # A/B 颜色
names(colors_layer2) <- c("A", "B")

colors_layer3 <- brewer.pal(7, "Pastel1")  # 7个染色体颜色
names(colors_layer3) <- paste0("chr", 1:7)

# 创建基础树图
p <- ggtree(flLTR_tree, layout = "circular", size = 0.3, color = "gray80") %<+% annot_df

# 添加三层注释
p <- p +
  # 第一层注释 (最内圈)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = Layer1),
    width = 0.2, 
    offset = 0.05
  ) +
  # 第二层注释 (中间圈)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = Layer2),
    width = 0.2, 
    offset = 0.17
  ) +
  # 第三层注释 (最外圈)
  geom_fruit(
    geom = geom_tile,
    mapping = aes(fill = Layer3),
    width = 0.2, 
    offset = 0.17
  ) +
  # 设置颜色标尺
  scale_fill_manual(
    name = "Annotation",
    values = c(colors_layer1, colors_layer2, colors_layer3),
    breaks = c(names(colors_layer1), names(colors_layer2), names(colors_layer3)),
    labels = c("Centromere", "Arm", "A type", "B type", 
               paste("Chromosome", 1:7))
  ) +
  # 添加图例
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)) +
  guides(fill = guide_legend(ncol = 1))

# 显示图形
print(p)
ggsave('../figure/figure2/flLTR.tree.pdf',plot = p,width = 5,height = 3)
