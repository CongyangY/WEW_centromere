library(tidyverse)
library(pastecs)
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


# 1. flLTR insertime -----------------------------------------------------
fl_insert <- read_delim('../data/fl-LTR/asm.np2.pass.list', delim = '\t', skip = 1,
                        col_names = FALSE) %>% 
  select(X1, X10, X11, X12) %>% 
  mutate(X1 = sub(":.*", "", X1),
         Sub = ifelse(grepl("A", X1),'A',ifelse(grepl("B", X1),'B','Other')))

ggplot(fl_insert) + geom_density(aes(x = X12))


density_adj <- density(fl_insert$X12)  
# 寻找显著峰和谷
tp <- turnpoints(density_adj$y)

# 提取所有候选峰及其对应的高度
all_peaks_x <- density_adj$x[which(tp$peaks)]     # 峰对应x值
all_peaks_y <- density_adj$y[which(tp$peaks)]     # 峰对应y值

# 按峰高度降序排列，选择前两个最大的峰
top_peaks <- order(all_peaks_y, decreasing = TRUE)[1:2]
peak_times <- all_peaks_x[top_peaks]
cat("Inferred burst peaks (MYA):", sort(peak_times/1e6), "\n")  # 按时间顺序输出



# 提取所有谷位置
valleys_x <- density_adj$x[which(tp$pits)]
# 寻找分割点：在两个主峰之间的第一个谷
if (length(valleys_x) > 0) {
  # 确定主峰的时间范围
  peak_range <- sort(peak_times)
  
  # 筛选出位于两个主峰之间的谷
  candidates <- valleys_x[valleys_x > peak_range[1] & valleys_x < peak_range[2]]
  
  # 如果存在候选谷，选择最近的（或通过其他逻辑选择）
  if (length(candidates) > 0) {
    split_time <- candidates[which.min(abs(candidates - mean(peak_range)))]
  } else {
    warning("No valley found between peaks, using midpoint")
    split_time <- mean(peak_range)
  }
} else {
  stop("No valleys detected in density distribution")
}

cat("Split time (MYA):", split_time/1e6, "\n")



fl_insert <- fl_insert %>% 
  mutate(Burst_Period = ifelse(X12<= split_time, "Burst2", "Burst1"))



ggplot(fl_insert) + geom_density(aes(x = X12/1e6)) + 
  geom_vline(xintercept = split_time/1e6, color = "red", linetype = "dashed", linewidth = 0.4) + 
  theme_Publication() + 
  labs(x = 'Insertion Time (Myr)', y = 'Density')
ggsave('../figure/Sf1/flLTR.insert.density.pdf',width = 3.5,height = 2.3)





# different burst period X2 -----------------------------------------------


GTest(table(fl_insert$Burst_Period, fl_insert$X10))
# G = 1570.9, X-squared df = 2, p-value < 2.2e-16

data <- table(fl_insert$Burst_Period, fl_insert$X10)
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v
# 0.0961424 

# 计算每个条形的高度（占比）
fl_insert_2 <- fl_insert %>%
  group_by(Burst_Period, X10) %>%
  summarise(count = n()) %>%
  mutate(prop = count / sum(count))

# 绘制 ggplot
ggplot(fl_insert_2, aes(x = X10, fill = X10)) +
  geom_bar(aes(y = prop), position = "dodge", stat = "identity") +  # 显示占比
  geom_text(aes(y = prop, label = scales::percent(prop)), 
            position = position_dodge(width = 0.9), vjust = 0.05, size = 3) +  # 添加占比标签
  facet_wrap(~ Burst_Period, scales = "free_x", 
             labeller = as_labeller(c(Burst1 = paste("Burst1:", round(peak_times[1]/1e6,2), "Mya"), 
                                      Burst2 = paste("Burst2:", round(peak_times[2]/1e6,2), "Mya")))) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "LTR Type", y = "Proportion of Elements", 
       title = "Dominant LTR Types in Two Burst Periods") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = 'None')

ggsave('../figure/Sf1/flLTR.burstPerid.bar.pdf',width = 4.5,height = 4)


# subgenome ---------------------------------------------------------------
GTest(table(filter(fl_insert,Sub!='Other')$Sub, filter(fl_insert,Sub!='Other')$X10))
# G = 338.95, X-squared df = 2, p-value < 2.2e-16

data <- table(filter(fl_insert,Sub!='Other')$Sub, filter(fl_insert,Sub!='Other')$X10)
chi_test <- GTest(data)
chi2 <- chi_test$statistic
n <- sum(data)
k <- min(nrow(data), ncol(data))  # k=2 (since rows=2)
v <- sqrt(chi2 / (n * (k - 1)))
v
# 0.0446614



