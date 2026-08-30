library(tidyverse)
library(ggpubr)
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






# 1.copy number -----------------------------------------------------------
df <- tibble::tibble(
  Chr = paste0("chr", 1:7, "B"),
  WEW = c(870, 151, 70, 6, 591, 1166, 579),
  CS  = c(406,  48, 69, 6, 802,  392, 604)
)


# df <- tibble::tibble(
#   Chr = paste0("chr", 1:7, "B"),
#   WEW = c(985, 172, 100, 8, 710, 1377, 657),
#   CS  = c(478,  56, 103, 8, 931,  454, 656)
# )
# df <- tibble::tibble(
#   Chr = paste0("chr", 1:7, "B"),
#   WEW = c(985,326,106,22,713,1583,746),
#   CS  = c(478,197,104,22,933,664,861)
# )

# 计算 log2 fold change
eps <- 1  # 避免除零
df <- df %>%
  mutate(
    log2FC = log2((WEW + eps) / (CS + eps)),
    dir = case_when(
      log2FC >  0.15 ~ "increase",   # 上升
      log2FC < -0.15 ~ "decrease",   # 下降
      TRUE ~ "no change"             # 变化小视为不变
    )
  ) %>%
  arrange(log2FC) %>%
  mutate(Chr = factor(Chr, levels = Chr))

# 调色
cols <- c("decrease" = "#fdae6b", "increase" = "#9ecae1", "no change" = "grey70")

# 绘图：横轴为 log2 fold change
ggplot(df, aes(x = log2FC, y = Chr, color = dir)) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%.2f", log2FC)),
            hjust = ifelse(df$log2FC > 0, -0.2, 1.2),
            vjust = 0.4, size = 3.2, color = "black") +
  scale_color_manual(values = cols, name = expression(Delta~"direction")) +
  scale_x_continuous(limits = c(-2.5, 2.5), breaks = seq(-2, 2, 1)) +
  theme_classic(base_size = 13) +
  theme(legend.position = "right") +
  labs(x = expression(log[2]*"(WEW / CS)"), y = NULL)
ggsave('../figure/figure5/WEW.tandem.peri10M.dot.pdf',height = 4,width = 4)
# 2. pheatmap ----------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)

# === 输入数据：单位用“百分比”更直观（例如 0.13 表示 0.13%）===
tetra <- 0.13                 # AABB 四倍体的 566 占比（%）
hexa  <- c(0.12,0.11,0.10,    # ← 把这行替换成你的十几个 AABBDD 样本（%）
           0.11,0.09,0.10,
           0.10,0.12,0.11,
           0.10,0.11,0.09)

df <- data.frame(group="AABBDD", value=hexa)
n   <- nrow(df)
below0.13 <- mean(df$value < tetra)

# 统计检验：AABBDD 是否小于 0.13%（单尾）
p_t  <- t.test(df$value, mu=tetra, alternative="less")$p.value
p_w  <- wilcox.test(df$value, mu=tetra, alternative="less", exact=FALSE)$p.value
title_txt <- sprintf("A-subgenome (within AABBDD): mean=%.3f%% (n=%d) vs AABB=%.2f%%; t=%.2g, Wilcoxon=%.2g; %.0f%%<%.2f%%",
                     mean(df$value), n, tetra, p_t, p_w, 100*below0.13, tetra)

# —— 图 A：基线比较（雨云图：小提琴+箱线+散点）——
pA <- ggplot(df, aes(x=group, y=value)) +
  geom_violin(fill="grey92", color="grey40", width=0.6, trim=FALSE) +
  geom_boxplot(width=0.12, outlier.shape = NA, fill="white") +
  geom_jitter(width=0.06, size=2.6, alpha=0.9) +
  geom_hline(yintercept = tetra, linetype="dashed", linewidth=0.7, color="firebrick") +
  annotate("label", x=1.0, y=tetra, label=sprintf("AABB = %.2f%%", tetra),
           vjust=-0.6, size=3.3, fill="white", label.size=0, color="firebrick") +
  scale_y_continuous(expand = expansion(mult=c(0.05,0.10))) +
  labs(x=NULL, y="Cent566 copy proportion in A-subgenome (%)") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size=11))

# —— 图 B：log2 比值分布（相对 AABB=0.13%）——
df$log2FC <- log2(df$value / tetra)
pB <- ggplot(df, aes(x=log2FC)) +
  geom_histogram(binwidth = 0.025, color="white") +
  geom_vline(xintercept = 0, linetype="dashed", color="grey30") +
  labs(x="log2(AABBDD / 0.13%)", y="Count") +
  theme_classic(base_size = 13)

ggpubr::ggarrange(pA, pB, ncol=2, labels=c("A","B")) %>%
  ggpubr::annotate_figure(top = ggpubr::text_grob(title_txt, face="bold", size=13))


# 1. WEW vs. hexploidy ----------------------------------------------------


# === 输入数据：单位用“百分比”更直观（例如 0.13 表示 0.13%）===
tetra <- 0.13501                 # AABB 四倍体的 566 占比（%）
hexa  <- c(0.09843,0.09951,0.12769,    # ← 把这行替换成你的十几个 AABBDD 样本（%）
           0.09624,0.09801,0.11468,
           0.07545,0.08964,0.10863,
           0.11418,0.09522,0.10986,
           0.10863,0.10003,0.11430,
           0.06984,0.07892,0.07275,
           0.08556)

habit <- c("Strong Winter","Semi Winter","Strong Winter",
           "Semi Winter","Spring","Spring",
           "Strong Winter","Winter","Strong Winter",
           "Winter","Strong Winter","Spring",
           "Strong Winter","Strong Winter","Winter",
           "Semi Winter","Winter","Winter",
           "Winter")
# 整理成数据框
df <- tibble(group = "AABBDD",
             value = hexa,
             habit = factor(habit, levels = c("Spring","Semi Winter","Winter","Strong Winter")))

# 统计（可选）：整体是否低于 AABB 基线
p_t <- t.test(df$value, mu = tetra, alternative = "less")$p.value
p_w <- wilcox.test(df$value, mu = tetra, alternative = "less", exact = FALSE)$p.value
title_txt <- sprintf("AABBDD (A-subgenome): n=%d, mean=%.4f%%  vs  AABB=%.4f%%  |  t=%.2g, Wilcoxon=%.2g",
                     nrow(df), mean(df$value), tetra, p_t, p_w)

# —— 图：总体箱线 + 彩色散点（春/冬性） + 基线 —— 
ggplot(df, aes(x = group, y = value)) +
  geom_violin(fill = "grey93", color = "grey60", width = 0.55, trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "grey30") +
  # ← 把散点放在最后一层
  geom_jitter(aes(color = habit), width = 0.07, size = 2.8, alpha = 0.95) +
  geom_hline(yintercept = tetra, linetype = "dashed", linewidth = 0.7, color = "firebrick") +
  annotate("label", x = 1, y = tetra, label = sprintf("AABB = %.3f%%", tetra),
           vjust = -0.6, size = 3.2, label.size = 0, fill = "white", color = "firebrick") +
  scale_color_manual(values = c("Spring"="#1B9E77","Semi Winter"="#D95F02",
                                "Winter"="#7570B3","Strong Winter"="#E7298A")) +
  labs(x = NULL, y = "Cent566 proportion（AB subgenome %）", color = "Ecotype") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(size = 11)) 
ggsave('../figure/figure5/tetra_hexa.box.pdf',height = 4.8,width = 5.2)


df <- tibble(value = hexa,
             diff  = value - tetra,
             log2FC = log2(value / tetra))
ggqqplot(df$log2FC)  + theme_Publication()
ggsave('../figure/figure5/tetra_hexa.qqplot.pdf',width = 3,height = 3)



# 2. similarity of tandem repeat ------------------------------------------


WEW_similar <- read_delim('../data/CENH3/cent566/566_in_WEWnew.slop10M.m8.similar',delim = '\t',col_names = F) %>% 
  mutate(Type = 'WEW')
CSGL_similar <- read_delim('../data/CENH3/cent566/566_in_CSGL.slop10M.m8.similar',delim = '\t',col_names = F) %>% 
  mutate(Type = 'CS')
TR_similar <- rbind(WEW_similar,CSGL_similar)

ggplot(TR_similar,aes(x = Type,y = X1)) + geom_boxplot(aes(fill = Type),outlier.alpha = 0) + 
  geom_signif(comparisons = list(c("WEW", "CS")),map_signif_level=TRUE) + 
  scale_fill_Publication() + theme_Publication()

ggplot(TR_similar,aes(x = Type,y = X1)) + geom_violin(aes(fill = Type),outlier.alpha = 0) + 
  geom_boxplot(aes(fill = Type),outlier.alpha = 0,width = 0.2) + 
  geom_signif(comparisons = list(c("WEW", "CS")),map_signif_level=TRUE) + 
  scale_fill_Publication() + theme_Publication()

ggsave('../figure/figure5/tetra_hexa.similar.violin.pdf',width = 3.6,height = 4)
  
ggplot(TR_similar) + stat_ecdf(aes(color = Type,x = X1)) + 
  theme_Publication() + scale_colour_Publication() + xlim(75,99)
ggsave('../figure/figure5/tetra_hexa.similar.ecdf.pdf',width = 3.6,height = 4)





ggplot(WEW_similar) + geom_density(aes(x = X1))
ggplot(CSGL_similar) + geom_density(aes(x = X1))



# 依赖
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(viridis)
})

# m8 -> 下三角热图矩阵函数
m8_to_tri_heatmap <- function(infile,
                              out_matrix = NULL,
                              out_pdf = NULL,
                              triangle = "lower",
                              palette = "viridis",
                              viridis_option = "C") {
  triangle <- match.arg(triangle)
  palette  <- match.arg(palette)
  
  # 读 m8（无表头），只用 q s pident（第1、2、3列）
  dt <- fread(infile, header = FALSE, sep = "\t", quote = "", data.table = TRUE)
  setnames(dt, paste0("V", seq_len(ncol(dt))))
  dt <- dt[, .(q = V1, s = V2, pident = as.numeric(V3))]
  
  # 同一对取最大 pident
  dt_max <- dt[, .(pident = max(pident, na.rm = TRUE)), by = .(q, s)]
  
  # ID 顺序
  ids <- sort(unique(c(dt_max$q, dt_max$s)))
  n <- length(ids)
  idx <- setNames(seq_along(ids), ids)
  
  # 构造对称矩阵
  mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
  diag(mat) <- 100
  if (nrow(dt_max) > 0) {
    for (k in seq_len(nrow(dt_max))) {
      i <- idx[[ dt_max$q[k] ]]; j <- idx[[ dt_max$s[k] ]]
      if (!is.na(i) && !is.na(j)) {
        v <- dt_max$pident[k]
        if (is.na(mat[i,j]) || v > mat[i,j]) { mat[i,j] <- v; mat[j,i] <- v }
      }
    }
  }
  
  # 可选导出矩阵
  if (!is.null(out_matrix)) {
    write.table(
      cbind(ID = ids, as.data.frame(round(mat, 2), check.names = FALSE)),
      file = out_matrix, sep = "\t", quote = FALSE, row.names = FALSE
    )
  }
  
  # 下/上三角数据
  df <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  colnames(df) <- c("i", "j", "pident")
  if (triangle == "lower") {
    df <- df[ match(df$i, ids) >= match(df$j, ids), ]
  } else {
    df <- df[ match(df$i, ids) <= match(df$j, ids), ]
  }
  df <- df[!is.na(df$pident), ]
  df$i <- factor(df$i, levels = ids)
  df$j <- factor(df$j, levels = ids)
  
  # 颜色范围：最小值到 100
  min_pid <- min(df$pident, na.rm = TRUE)
  if (!is.finite(min_pid)) min_pid <- 100
  
  # 调色
  fill_scale <- if (palette == "viridis") {
    scale_fill_viridis(option = viridis_option, direction = 1,
                       limits = c(min_pid, 100), name = "pident (%)")
  } else {
    scale_fill_gradientn(
      colours = c("#1f77b4", "#a1d99b", "#fee08b", "#f46d43", "#a50026"),
      limits = c(min_pid, 100), name = "pident (%)"
    )
  }
  
  # 画图（去掉坐标轴文字和刻度）
  p <- ggplot(df, aes(x = j, y = i, fill = pident)) +
    geom_tile() +
    fill_scale +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x  = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid   = element_blank()
    )
  
  # 可选导出 PDF
  if (!is.null(out_pdf)) {
    ggsave(out_pdf, p,
           width = max(6, min(12, n * 0.25)),
           height = max(6, min(12, n * 0.25)))
  }
  
  invisible(list(matrix = mat, plot = p, ids = ids, min_pid = min_pid))
}

# ===== 示例 =====
# res <- m8_to_tri_heatmap(
#   infile = "../data/CENH3/cent566/566_in_CSGL.slop10M.chr5A.m8",
#   out_matrix = "chr5A_matrix.tsv",
#   out_pdf = "chr5A_tri_heatmap.pdf",
#   triangle = "lower",
#   palette = "viridis",         # 红绿色盲友好
#   viridis_option = "C"         # 可选 "A"/"B"/"C"/"D"/"E"
# )
# res$plot







