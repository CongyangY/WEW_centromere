library(tidyverse)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(scales)

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






# centromere size ----------------------------------------------------------
cent.size <- data.frame(
  Chr = c(paste('chr',1:7,'A',sep = ''),paste('chr',1:7,'B',sep = '')),
  WEW = c(7.2,7.5,7.3,9.23,7.2,8.4,8.1,8,9.2,6.9,10,8.5,7.5,7.9),
  Langdon = c(7,6.2,7.1,7.5,5.4,7.6,7.9,7.5,8.1,6.6,9.5,7.2,7.4,7.7) 
) %>% melt(id = c('Chr'))


ggplot(cent.size) + 
  geom_bar(aes(x = Chr,y = value,group = variable,fill = variable),
           stat = 'identity',position = 'dodge',width = 0.6) + 
  theme_Publication() + scale_fill_Publication('Set2') + 
  theme(
    axis.text.x = element_text(angle = 60,hjust = 1)
  )
ggsave('../figure/Figure4/cent.size.bar.pdf',width =4,height = 4.8)


chr = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B","chr4A","chr4B",
        "chr5A","chr5B","chr6A","chr6B","chr7A","chr7B")

Langdon_rep1 = c(7.4,6.0,6.8,6.9,4.0,7.2,6.4,7.5,6.8,6.0,8.5,7.5,6.4,6.0)
Langdon_rep2 = c(8.0,6.8,7.2,7.5,7.2,7.9,7.5,8.0,8.3,7.0,9.1,9.0,7.2,8.5)
Langdon_rep3 = c(8.0,7.0,7.0,7.5,7.2,8.1,7.4,8.1,8.2,7.1,9.1,9.0,7.4,8.2)

Zavitan_rep1 = c(7.0,6.3,7.0,7.9,7.1,8.0,8.4,7.6,9.1,6.5,9.1,8.4,7.3,7.6)
Zavitan_rep2 = c(7.2,7.4,7.1,8.7,7.2,8.5,8.7,7.7,9.4,7.2,9.4,8.3,7.6,7.0)

cent.size <- data.frame(
  Chr = chr,
  Langdon_rep1 = Langdon_rep1,
  Langdon_rep2 = Langdon_rep2,
  Langdon_rep3 = Langdon_rep3,
  Zavitan_rep1 = Zavitan_rep1,
  Zavitan_rep2 = Zavitan_rep2
) %>% 
  melt(id = c("Chr"))

cent.size$Chr <- factor(cent.size$Chr, levels = chr)

cent.size$variable <- factor(
  cent.size$variable,
  levels = c("Zavitan_rep1", "Zavitan_rep2",
             "Langdon_rep1", "Langdon_rep2", "Langdon_rep3")
)

p <- ggplot(cent.size) + 
  geom_bar(
    aes(x = Chr, y = value, group = variable, fill = variable),
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.7
  ) + 
  theme_Publication() + 
  scale_fill_Publication("Set2") + 
  labs(
    x = NULL,
    y = "Centromere size (Mb)",
    fill = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),
    legend.position = "top"
  )

p



library(ggpubr)
ggplot(cent.size,aes(x = variable,y = value)) + 
  geom_boxplot(aes(fill = variable)) + 
  geom_signif(comparisons = list(c("WEW", "Langdon")),
              map_signif_level=TRUE,test =  "wilcox.test",test.args = list(paired = TRUE)) +
  theme_Publication() + scale_fill_Publication('Set2')
ggsave('../figure/Figure4/cent.size.box.pdf',height = 4.6,width = 2.6)



# proportion of sv --------------------------------------------------------


# 1️⃣ 合并类别（不考虑 HDR）
syri_m <- syri.clean %>%
  filter(class != "HDR") %>% 
  mutate(
    class6 = case_when(
      class %in% c("INV","INVTR","INVDP") ~ "Inversion",
      class %in% c("TRANS")               ~ "Translocation",
      class %in% c("DUP","TDM")           ~ "Duplication",
      class %in% c("CPL","CPG")           ~ "Complex",
      class %in% c("SYN")                 ~ "Syntenic",
      TRUE                                ~ "Other"
    ),
    class6 = factor(class6,
                    levels = c("Inversion","Translocation","Duplication",
                               "Complex","Syntenic","Other")),
    subg = factor(subg, levels = c("A","B"))
  )

# 🎨 使用 ColorBrewer 的 Dark2 调色板
# 取 6 种颜色并与类别一一对应
cols6 <- setNames(brewer.pal(6, "Set2"),
                  c("Inversion","Translocation","Duplication",
                    "Complex","Syntenic","Other"))

# 2️⃣ 累计长度加权的 100% 堆叠柱图（按 A/B 亚基因组）
df_stack_bp <- syri_m %>%
  group_by(subg, class6) %>%
  summarise(bp = sum(length), .groups = "drop") %>%
  group_by(subg) %>%
  mutate(prop = bp / sum(bp))

p_stack <- ggplot(df_stack_bp, aes(x = subg, y = prop, fill = class6)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = cols6, name = "SV type") +
  labs(x = "Subgenome", y = "Composition (% by affected bp)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text = element_text(color = "black")
  )
p_stack
ggsave('../figure/Figure4/SV.proportion.stack_bar.pdf',width = 4,height = 4)

# 3️⃣ 累计长度加权的 100% 饼图（全基因组）
df_pie_bp <- syri_m %>%
  group_by(class6) %>%
  summarise(bp = sum(length), .groups = "drop") %>%
  mutate(prop = bp / sum(bp),
         label = paste0(class6, " (", scales::percent(prop, accuracy = 0.1), ")"))

p_pie <- ggplot(df_pie_bp, aes(x = "", y = prop, fill = class6)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cols6, name = "SV type") +
  theme_void(base_size = 12)
p_pie
ggsave('../figure/Figure4/SV.proportion.pieplot.pdf',width = 4,height = 4)


# frequency of structure variation ----------------------------------------
syri.clean <- read_delim("../figure/Figure4/syri.out.clean", delim = "\t", col_names = F) %>%
  # ① 合并类别
  mutate(class = case_when(
    X11 %in% c("SYN", "SYNAL") ~ "SYN",
    X11 %in% c("INV", "INVAL") ~ "INV",
    X11 %in% c("INVTR", "INVTRAL", "INVDP", "INVDPAL") ~ "INVTR",
    X11 %in% c("TRANS", "TRANSAL") ~ "TRANS",
    X11 %in% c("DUP", "DUPAL") ~ "DUP",
    X11 == "HDR" ~ "HDR",
    X11 == "TDM" ~ "TDM",
    X11 == "CPL" ~ "CPL",
    X11 == "CPG" ~ "CPG",
    TRUE ~ "OTHER"
  )) %>%
  # ② 再过滤掉镜像记录（防止重复）
  filter(!str_detect(X11, "AL$")) %>%
  # ③ 再过滤掉 NOTAL（未比对）
  filter(X11 != "NOTAL") %>%  mutate(subg = str_extract(X1, "[AB]$")) %>% 
  mutate(length = as.numeric(X3) - as.numeric(X2) + 1)

syri.clean.sv <- filter(syri.clean,X11 %in% c('INV','TRANS','DUP'))
## 设定类别顺序（可按需要调整）
# 1) 汇总
df_sum <- syri.clean.sv %>%
  group_by(subg) %>%
  summarise(
    count = n(),
    median_len = median(length, na.rm = TRUE),
    .groups = "drop"
  )
# 缩放系数（右轴）
s <- max(df_sum$count, na.rm=TRUE) / max(df_sum$median_len, na.rm=TRUE)

# 红绿色盲友好色系（Set2风格）
cols_bar  <- c("A"="#1B9E77", "B"="#D95F02")   # bar 主色
cols_point <- c("A"="#0C7C59", "B"="#B44E03")  # 稍深一点（同色系更饱和）用于点线



ggplot(df_sum, aes(x=subg)) +
  # 柱子（主色 + 黑边）
  geom_col(aes(y=count, fill=subg), width=0.6, color="black", alpha=0.9) +
  # 点与连线（更深的同系色）
  geom_point(aes(y=median_len * s, color=subg), size=3) +
  geom_line(aes(y=median_len * s, group=1), color="gray20", linewidth=0.9, linetype="solid") +
  # 手动配色
  scale_fill_manual(values=cols_bar) +
  scale_color_manual(values=cols_point) +
  scale_y_continuous(
    name = "Number of structural variants",
    sec.axis = sec_axis(
      ~ . / s,
      name = "Median size (bp)",
      labels = label_number(scale_cut = cut_si("bp"))
    )
  ) +
  labs(title = " ", x = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    axis.title.y.left  = element_text(color="gray20"),
    axis.title.y.right = element_text(color="gray20"),
    plot.title = element_text(face="bold", hjust=0.5)
  ) + theme_Publication()
ggsave('../figure/Figure4/SV.bar.dot.pdf',height = 4.5,width = 3.5)


ggplot(syri.clean.sv) + geom_boxplot(aes(x = subg,y = log10(length)))



syri.clean.sv <- syri.clean.sv %>%
  mutate(
    # 先挑出一个“像染色体名”的列作为来源：优先 X1（参考端），否则退到 X6（查询端）
    chr_src = case_when(
      str_detect(X1, "^[Cc]hr\\d+[AB]$") ~ X1,
      str_detect(X6, "^[Cc]hr\\d+[AB]$") ~ X6,
      TRUE ~ X1  # 实在不匹配就用 X1，避免 NA
    ),
    # 去掉前缀 chr/Chr，得到 “1A/2B” 这样的短名
    chr_short = str_replace(chr_src, "^[Cc]hr", ""),
    # 提取亚基因组与编号
    subg = str_extract(chr_short, "[AB]$") %>% toupper(),
    chr_id = suppressWarnings(as.integer(str_extract(chr_short, "\\d+")))
  )

# 现在这步就能跑
syri.clean.sv.count <- syri.clean.sv %>%
  count(chr_short, name = "SV") %>%
  mutate(
    subg = str_extract(chr_short, "[AB]$") %>% toupper(),
    chr_id = suppressWarnings(as.integer(str_extract(chr_short, "\\d+")))
  ) %>%
  filter(!is.na(subg), !is.na(chr_id), chr_id %in% 1:7) %>%
  mutate(subg = factor(subg, levels = c("A","B")))

ggplot(syri.clean.sv.count,aes(x = subg,y = SV)) + geom_boxplot(aes(fill = subg)) + 
  theme_Publication() + scale_fill_Publication() + 
  geom_signif(comparisons = list(c("A", "B")), test =  "wilcox.test",test.args = list(paired = TRUE,alternative = 'less'))
ggsave('../figure/Figure4/SV.subgenome.frequency.box.pdf',width = 2.5,height = 4.5)

ggplot(syri.clean.sv,aes(x = subg,y = log10(length))) + geom_boxplot(aes(fill = subg)) + 
  theme_Publication() + scale_fill_Publication() + 
  geom_signif(comparisons = list(c("A", "B")),test =  "wilcox.test")
ggsave('../figure/Figure4/SV.subgenome.median_length.box.pdf',width = 2.5,height = 4.5)



# frequency of different distance -----------------------------------------------------------------
centro <- read_delim("../figure/Figure4/WEW_new.cent.relative_pos.bed", delim="\t",
                     col_names=c("chr", "cen_start", "cen_end")) %>%
  mutate(cen_mid = (cen_start + cen_end)/2)
sv <- syri.clean %>%
  filter(!str_detect(X11, "AL$"), X11 != "NOTAL") %>%
  mutate(
    chr = X1,
    start = as.numeric(X2),
    end   = as.numeric(X3),
    mid   = (start + end)/2,
    chr_short = str_replace(chr, regex("^chr", ignore_case = TRUE), ""),
    subg = str_extract(chr, "[AB]") %>% toupper()
  ) %>%
  filter(!is.na(subg), is.finite(mid))

## 若只看 INV，解注释下一行；或改为 c("INV","INVTR") 等
# sv <- sv %>% filter(class == "INV")

## 合并着丝粒中心并算距离(Mb)
sv <- sv %>%
  left_join(centro %>% select(chr, cen_mid), by = "chr") %>%
  mutate(dist_to_cen_Mb = abs(mid - cen_mid) / 1e6)

##---------------------------
## 分箱：0–2 / 2–5 / 5–10 / 10–25 Mb
## 箱宽：2, 3, 5, 15（Mb）
##---------------------------
breaks <- c(0, 2, 5, 10, 25)
labels <- c("0–2", "2–5", "5–10", "10–25")
bin_widths <- c(`0–2`=2, `2–5`=3, `5–10`=5, `10–25`=15)

sv <- sv %>%
  filter(dist_to_cen_Mb >= 0, dist_to_cen_Mb <= 25) %>%  # 只分析 0–25 Mb 范围
  mutate(dist_bin = cut(dist_to_cen_Mb, breaks = breaks, labels = labels, include.lowest = TRUE))

##---------------------------
## 逐“染色体×距离箱”的密度：count / bin_width
## 然后在 A/B 上对每条染色体的密度取均值 ± SE
##---------------------------
per_chr_density <- sv %>%
  count(subg, chr_short, dist_bin, name = "count") %>%
  complete(subg, chr_short, dist_bin = factor(labels, levels = labels), fill = list(count = 0)) %>% 
  mutate(bin_w = as.numeric(bin_widths[as.character(dist_bin)]),
         density = count / bin_w)  # SV per Mb

## 聚合到亚基因组（A/B）：每个距离箱内，对“每条染色体密度”取均值 ± SE
by_subg <- per_chr_density %>%
  group_by(subg, dist_bin) %>%
  summarise(
    n_chr = n_distinct(chr_short),
    mean_density = mean(density, na.rm = TRUE),
    sd_density = sd(density, na.rm = TRUE),
    se_density = ifelse(n_chr > 0, sd_density / sqrt(n_chr), NA_real_),
    .groups = "drop"
  )

## 颜色（Okabe & Ito）
cols_bar <- c("A"="#1B9E77", "B"="#D95F02")

## 图1：四个箱的密度对比（均值±SE；A/B 并排）
p_bins <- ggplot(by_subg, aes(x = dist_bin, y = mean_density, fill = subg)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = mean_density - se_density,
                    ymax = mean_density + se_density),
                position = position_dodge(width = 0.7), width = 0.25) +
  scale_fill_manual(values = cols_bar) +
  labs(x = "Distance to centromere (Mb)",
       y = "SV density (events per Mb)",
       title = "SV density vs. distance to centromere (per chromosome, mean ± SE)") +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank())

##---------------------------
## 可选：把四个箱合并成三组：
## 近心：0–2；中间：2–10；远心：10–25
## 先对每条染色体在组内“加权合并”：总count / 总宽度 = (sum count)/(sum bin_width)
## 然后再在 A/B 上取均值 ± SE
##---------------------------
group_map <- tibble(
  dist_bin = factor(labels, levels = labels),
  group3 = c("proximal (0–2)", "intermediate (2–10)", "intermediate (2–10)", "distal (10–25)")
)

per_chr_group3 <- per_chr_density %>%
  left_join(group_map, by = "dist_bin") %>%
  group_by(subg, chr_short, group3) %>%
  summarise(
    count_sum = sum(count),
    width_sum = sum(bin_w),
    density = ifelse(width_sum > 0, count_sum / width_sum, NA_real_),
    .groups = "drop"
  )

by_subg_group3 <- per_chr_group3 %>%
  group_by(subg, group3) %>%
  summarise(
    n_chr = n_distinct(chr_short),
    mean_density = mean(density, na.rm = TRUE),
    sd_density = sd(density, na.rm = TRUE),
    se_density = ifelse(n_chr > 0, sd_density / sqrt(n_chr), NA_real_),
    .groups = "drop"
  )

p_group3 <- ggplot(by_subg_group3, aes(x = group3, y = mean_density, fill = subg)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = mean_density - se_density,
                    ymax = mean_density + se_density),
                position = position_dodge(width = 0.7), width = 0.25) +
  scale_fill_manual(values = cols_bar) +
  labs(x = "Distance group",
       y = "SV density (events per Mb)") +
  theme_bw(base_size = 13) +
  theme(legend.title = element_blank(),
        panel.grid.minor = element_blank()) + theme_Publication() + scale_fill_Publication()

# 打印图
p_bins
p_group3
ggsave('../figure/Figure4/SV.distance.bar.pdf',width = 5.2,height = 5.2)

# per_chr_group3: 每条染色体在 group3 内的加权密度（总count/总宽度）
library(dplyr)
library(stringr)
library(tidyr)

## per_chr_group3: 每条染色体×组（proximal/intermediate/distal）的密度
## 列：subg, chr_short, group3, density

# 1) 先去重汇总：同一 subg+chr+group3 若有多行，先求平均密度
per_chr_group3_dedup <- per_chr_group3 %>%
  mutate(chr_id = as.integer(str_extract(chr_short, "\\d+"))) %>%
  filter(!is.na(chr_id)) %>%
  group_by(subg, chr_id, group3) %>%
  summarise(density = mean(density, na.rm = TRUE), .groups = "drop")

# 2) A/B 成对整理
wide_pair3 <- per_chr_group3_dedup %>%
  select(subg, chr_id, group3, density) %>%
  pivot_wider(names_from = subg, values_from = density) %>%
  filter(!is.na(A), !is.na(B))

# 3) 逐 group3 做配对 Wilcoxon（A<B） + 双侧 + 效应量 + FDR
tests_group3 <- wide_pair3 %>%
  group_by(group3) %>%
  summarise(
    n_pairs = n(),
    A_mean = mean(A),  B_mean = mean(B),
    A_median = median(A), B_median = median(B),
    diff_median = median(A - B),
    p_wilcox_less = suppressWarnings(
      wilcox.test(B, A, paired = TRUE, alternative = "less")$p.value
    ),
    p_wilcox_two  = suppressWarnings(
      wilcox.test(A, B, paired = TRUE)$p.value
    ),
    # rank-biserial 效应量（负值表示 A<B 倾向）
    r_rb = {
      d <- B - A
      ranks <- rank(abs(d))
      Wpos <- sum(ranks[d > 0])
      n <- length(d)
      1 - (2 * Wpos) / (n * (n + 1))
    },
    .groups = "drop"
  ) %>%
  mutate(
    p_wilcox_less_FDR = p.adjust(p_wilcox_less, method = "BH"),
    p_wilcox_two_FDR  = p.adjust(p_wilcox_two,  method = "BH")
  )

tests_group3




