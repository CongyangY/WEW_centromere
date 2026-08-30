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






# 1. insertion time -------------------------------------------------------

# 1.1 WEW -----------------------------------------------------------------



Cereba_new <- read_delim('../data/Cereba_new/Cereba.new_peaks.flLTR.bed',
                         delim = '\t',col_names = F) %>% 
  mutate(Type = 'Distal',Time = (1-X5/100)/(2*1.3e-8)/1e6)
Cereba_old <- read_delim('../data/fl-LTR/fl_LTR.cent.bed',
                         delim = '\t',col_names = F) %>% 
  filter(X1 %in% c('chr1A', 'chr2A', 'chr3A')) %>% 
  mutate(Type = 'Cent',Time = (1-X5/100)/(2*1.3e-8)/1e6)

Cereba <- rbind(Cereba_new,Cereba_old)

ggplot(Cereba,aes(x = Type, y = Time,fill = Type)) + geom_boxplot(outlier.size = 0.2) +
  theme_Publication() +
  facet_wrap(~ X1,ncol = 1) + 
  scale_fill_Publication('Paired')+ 
  theme(strip.background = element_blank())


ggplot(filter(Cereba,X1 == 'chr1A'),aes(x = Type, y = Time,fill = Type)) + geom_boxplot(outlier.size = 0.2) +
  theme_Publication() +
  scale_fill_Publication('Paired')+
  geom_signif(comparisons = list(c("Distal", "Cent")),
              map_signif_level=TRUE,y_position = 5.5) 
ggsave('../figure/figure3/Ceraba_insertTime.box.chr1A.pdf',width = 2.6,height = 3.5)






# 1.2 Lachesis ------------------------------------------------------------

Lachesis_newPeaks.1A <- read_delim('../data/Cereba_new/Lachesis_1A.second_trans.flLTR.bed',
                         delim = '\t',col_names = F) %>% 
  mutate(Type = 'Distal',Time = (1-X5/100)/(2*1.3e-8)/1e6)
Lachesis_cent.1A <- read_delim('../data/Cereba_new/Lachesis_1A.cent.flLTR.bed',
                         delim = '\t',col_names = F) %>% 
  mutate(Type = 'Cent',Time = (1-X5/100)/(2*1.3e-8)/1e6)
Lachesis.1A <- rbind(Lachesis_newPeaks.1A,Lachesis_cent.1A)
ggplot(Lachesis.1A,aes(x = Type, y = Time,fill = Type)) + geom_boxplot(outlier.size = 0.2) +
  theme_Publication() +
  scale_fill_Publication('Paired')+
  geom_signif(comparisons = list(c("Distal", "Cent")),
              map_signif_level=TRUE,y_position = 5.5) 
ggsave('../figure/Sf3/Lachesis_CRW.box.chr1A.pdf',width = 2.6,height = 3.5)



# 2. Cereba -------------------------------------------------------


library(ape)
library(stringr)
library(dplyr)
library(tidyr)


tree_file <- "../data/Cereba_new/Lachesis_1A.3_region.flLTR.aligned.fa.nwk"  # <-- 改成你的 nwk 文件
n_perm    <- 2000                     # 置换次数（可调）


suppressPackageStartupMessages({
  library(ape)
  library(stringr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

#  读树 & 三类标签
tr   <- read.tree(tree_file)
tips <- tr$tip.label

class_of <- function(x){
  # 先 second，再 cent，再 arm
  if (str_detect(x, regex("second", ignore_case = TRUE))) return("second")
  if (str_detect(x, regex("cent", ignore_case = TRUE))) return("cent")
  if (str_detect(x, regex("arm",      ignore_case = TRUE))) return("arm")
  return(NA_character_)
}
classes <- vapply(tips, class_of, character(1))
if (anyNA(classes)) {
  stop("有 tip 未匹配到 arm/cent/second，请检查命名。")
}
names(classes) <- tips

second_tips <- names(classes)[classes == "second"]
cent_tips   <- names(classes)[classes == "cent"]
arm_tips    <- names(classes)[classes == "arm"]
stopifnot(length(second_tips) > 0, length(cent_tips) > 0, length(arm_tips) > 0)

# ================== 距离矩阵
D <- cophenetic.phylo(tr)
diag(D) <- Inf

# ================== 最近的“非-second”最近邻 
nonsecond_tips <- setdiff(tips, second_tips)
nearest_nonsecond_tip <- vapply(second_tips, function(ti){
  cand <- nonsecond_tips
  cand[which.min(D[ti, cand])]
}, character(1))
nearest_nonsecond_class_obs <- classes[nearest_nonsecond_tip]

# 观测比例：second 的最近“非-second”最近邻为 cent 的比例
n_cent <- sum(nearest_nonsecond_class_obs == "cent")
n_arm  <- sum(nearest_nonsecond_class_obs == "arm")
n_tot  <- length(nearest_nonsecond_class_obs)
prop_cent_obs2 <- n_cent / n_tot

# Binomial 95% CI（Clopper-Pearson）
bt <- binom.test(n_cent, n_tot, p = 0.5, alternative = "greater")  # 也给出 P(>0.5)
ci95 <- bt$conf.int  # 注意：在 one-sided 下会返回单侧区间，示意这里也可用

# 置换检验：在非-second（arm+cent）间打乱标签
set.seed(42)
pool_labels <- classes[nonsecond_tips]
prop_cent_perm2 <- replicate(n_perm, {
  perm_labels <- classes
  perm_labels[nonsecond_tips] <- sample(pool_labels, length(pool_labels), FALSE)
  mean(perm_labels[nearest_nonsecond_tip] == "cent")
})
p_perm2 <- (sum(prop_cent_perm2 >= prop_cent_obs2) + 1) / (n_perm + 1)

# ================== second -> cent / arm 的最小距离 & 配对检验
min_to_cent <- apply(D[second_tips, cent_tips, drop = FALSE], 1, min)
min_to_arm  <- apply(D[second_tips, arm_tips,  drop = FALSE], 1, min)
wil <- suppressWarnings(wilcox.test(min_to_cent, min_to_arm,
                                    paired = TRUE, alternative = "less"))
delta_med <- median(min_to_arm - min_to_cent)  # >0 表示 second 更接近 cent

# 汇总表（可导出补充表）
summary_df <- tibble(
  second_tip = second_tips,
  nn_tip     = nearest_nonsecond_tip,
  nn_class   = nearest_nonsecond_class_obs,
  d_min_cent = min_to_cent[second_tips],
  d_min_arm  = min_to_arm[second_tips]
)

# ================== 图 A：最近非-second 邻居（柱状
bar_df <- tibble(
  class = c("cent", "arm"),
  count = c(n_cent, n_arm)
) %>%
  mutate(class = factor(class, levels = c("cent","arm")))

annot_line1 <- sprintf("Proportion(second→cent) = %.3f", prop_cent_obs2)
# 用双侧 CI 展示更直观；如需单侧可改 bt$conf.int
ci_bt <- binom.test(n_cent, n_tot)$conf.int
annot_line2 <- sprintf("95%% CI = [%.2f, %.2f]; n = %d", ci_bt[1], ci_bt[2], n_tot)
annot_line3 <- sprintf("Permutation P = %.3g", p_perm2)

pA <- ggplot(bar_df, aes(x = class, y = count,fill = class)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = count), vjust = -0.5, size = 4) +
  labs(x = "Nearest non-second neighbor", y = "Count (second tips)",
       title = "Nearest non-second neighbor of 'second' tips") +
  annotate("text", x = 1.5, y = max(bar_df$count) * 1.15, label = annot_line1, size = 4) +
  annotate("text", x = 1.5, y = max(bar_df$count) * 1.08, label = annot_line2, size = 4) +
  annotate("text", x = 1.5, y = max(bar_df$count) * 1.01, label = annot_line3, size = 4) +
  coord_cartesian(ylim = c(0, max(bar_df$count) * 1.25)) +
  theme_classic(base_size = 13) + 
  scale_fill_Publication('Set2')
print(pA)
ggsave("../figure/figure3/supp_second_nearest_nonsecond_bar.pdf", pA, width = 5, height = 4)

# ================== 图 B：成对箱线（second→cent vs arm 的最小距离）
long_df <- summary_df %>%
  select(second_tip, d_min_cent, d_min_arm) %>%
  pivot_longer(cols = c(d_min_cent, d_min_arm),
               names_to = "target", values_to = "distance") %>%
  mutate(target = recode(target, d_min_cent = "cent", d_min_arm = "arm"),
         target = factor(target, levels = c("cent","arm")))

y_max <- max(long_df$distance) * 1.15
annotB <- sprintf("Paired Wilcoxon P = %.3g; median(arm − cent) = %.3f",
                  wil$p.value, delta_med)

pB <- ggplot(long_df, aes(x = target, y = distance)) +
  geom_line(alpha = 0.15,aes(group = second_tip)) +
  geom_boxplot(aes(fill=target),width = 0.4, outlier.alpha = 0.25) +
  # 如想要小提琴+箱线，替换为：
  # geom_violin(trim = TRUE, alpha = 0.2) + geom_boxplot(width = 0.15, outlier.alpha = 0.25)
  labs(x = "Target lineage", y = "Minimum patristic distance from 'second'",
       title = "Minimum distances: 'second' → cent vs arm") +
  annotate("text", x = 1.5, y = y_max, label = annotB, size = 4) +
  coord_cartesian(ylim = c(0, y_max * 1.02)) +
  theme_classic(base_size = 13) + 
  scale_fill_Publication('Set2')
print(pB)
ggsave("../figure/figure3/supp_second_min_distance_boxpaired.pdf", pB, width = 6, height = 4.2)

# ================== 控制台小结（便于复制到补充材料）
cat("\n=== Summary for Supplementary text ===\n")
cat(sprintf("[Nearest non-second] cent = %d, arm = %d, n = %d; proportion = %.3f; 95%% CI = [%.3f, %.3f]; Permutation P = %.3g\n",
            n_cent, n_arm, n_tot, prop_cent_obs2, ci_bt[1], ci_bt[2], p_perm2))
cat(sprintf("[Paired Wilcoxon] P = %.3g; median(arm - cent) = %.3f\n",
            wil$p.value, delta_med))
# 如需导出明细表：
# write.csv(summary_df, "supp_second_phylo_summary.csv", row.names = FALSE)





# 3. solo/intact  ratio ----------------------------------------------------------

CSGL.1A.intact.bg <- read_delim('../data/Cereba_new/CSGL.1A.intact.bg',delim = '\t',col_names = F)
CSGL.1A.solo.bg <- read_delim('../data/Cereba_new/CSGL.1A.solo.bg',delim = '\t',col_names = F)
CSGL.1A.ratio <- data.frame(
  Ratio = CSGL.1A.solo.bg$X4/(CSGL.1A.solo.bg$X4+CSGL.1A.intact.bg$X4),
  Type = 'Hexaploid'
)


WEW.1A.intact.bg <- read_delim('../data/Cereba_new/WEW.1A.intact.bg',delim = '\t',col_names = F)
WEW.1A.solo.bg <- read_delim('../data/Cereba_new/WEW.1A.solo.bg',delim = '\t',col_names = F)
WEW.1A.ratio <- data.frame(
  Ratio = WEW.1A.solo.bg$X4/(WEW.1A.solo.bg$X4+WEW.1A.intact.bg$X4),
  Type = 'Tetraploid'
)

Lachesis.1A.intact.bg <- read_delim('../data/Cereba_new/Lachesis.1A.intact.bg',delim = '\t',col_names = F)
Lachesis.1A.solo.bg <- read_delim('../data/Cereba_new/Lachesis.1A.solo.bg',delim = '\t',col_names = F)
Lachesis.1A.ratio <- data.frame(
  Ratio = Lachesis.1A.solo.bg$X4/(Lachesis.1A.solo.bg$X4+Lachesis.1A.intact.bg$X4),
  Type = 'Diploid'
)


all.1A.ratio <- rbind(CSGL.1A.ratio,WEW.1A.ratio,Lachesis.1A.ratio)
all.1A.ratio$Type <- factor(all.1A.ratio$Type,levels = c('Diploid','Tetraploid','Hexaploid'))

ggplot(all.1A.ratio,aes(x = Type,y = Ratio,fill = Type)) + geom_boxplot() + 
  theme_Publication() + scale_fill_Publication('Set2') + 
  labs(x = 'Ploidy level',y = 'Proportion of solo-LTRs') + 
  theme(
    legend.position = 'None'
  ) + 
  geom_signif(comparisons = list(c("Tetraploid", "Hexaploid")), y_position = 0.4,
                  map_signif_level=TRUE) +
  geom_signif(comparisons = list(c("Diploid", "Tetraploid")), y_position = 0.5,
              map_signif_level=TRUE) +
  geom_signif(comparisons = list(c("Diploid", "Hexaploid")), y_position = 0.6,
              map_signif_level=TRUE)
ggsave('../figure/figure3/Cereba.soloLTRs.box.pdf',width = 3.65,height = 4.35)
  


# 4. 1A second LTR age -------------------------------------------------------
# —— 依赖包 —— #
library(tidyverse)
library(binom)      # Wilson 置信区间
library(FSA)        # Dunn 检验（如需）
# install.packages("logistf")
library(logistf)    # Firth logistic（小样本稳健）


################################### 1) 读入、计算 Age（Myr）、清洗
# —— 参数：碱基替换率（与正文一致）—— #
mu <- 1.3e-8  # 每位点每年
to_age <- function(id_percent, mu = 1.3e-8) {
  # id_percent 为百分数（如 98.7）
  (1 - id_percent/100) / (2*mu) / 1e6  # Myr
}

# —— 读入 & 计算 Age —— #
Lachesis.1A.age <- readr::read_delim('../data/Cereba_new/Lachesis.1A.Cereba.flLTR.list',
                                     col_names = FALSE, delim = '\t') %>%
  mutate(Type = 'Diploid', Age = to_age(X1, mu))

WEW.1A.age <- readr::read_delim('../data/Cereba_new/WEW.1A.Cereba.flLTR.list',
                                col_names = FALSE, delim = '\t') %>%
  mutate(Type = 'Tetraploid', Age = to_age(X1, mu))

CSGL.1A.age <- readr::read_delim('../data/Cereba_new/CSGL.1A.Cereba.flLTR.list',
                                 col_names = FALSE, delim = '\t') %>%
  mutate(Type = 'Hexaploid', Age = to_age(X1, mu))

all.1A <- bind_rows(Lachesis.1A.age, WEW.1A.age, CSGL.1A.age) %>%
  filter(!is.na(X1), X1 > 0, X1 <= 100, !is.na(Age), is.finite(Age)) %>%
  mutate(Type = factor(Type, levels = c('Diploid','Tetraploid','Hexaploid')))


########################## “按年龄分层”的主图（柱图：各倍性在不同年龄层的占比 + Wilson CI）
# —— 按年龄分层 —— #
all.1A <- all.1A %>%
  mutate(AgeBin = cut(Age,
                      breaks = c(-Inf,  0.50, 2.00, Inf),
                      labels = c('Post-4x young (0.00–0.50)',
                                 'Late young (0.50–2.00)',
                                 'Old (>2.00)'),
                      right = TRUE, ordered_result = TRUE))

# —— 统计每组在各年龄层的占比 + Wilson CI —— #
tab_agebin <- all.1A %>%
  count(Type, AgeBin, name = "k") %>%
  group_by(Type) %>%
  mutate(n_total = sum(k)) %>%
  ungroup() %>%
  mutate(prop = k / n_total) %>%
  group_by(Type) %>%
  group_modify(~{
    ci <- binom::binom.confint(x = .x$k, n = .x$n_total, methods = "wilson")
    bind_cols(.x, select(ci, lower, upper))
  }) %>% ungroup()

# —— 作图：各年龄层的占比（分组柱 + Wilson CI） —— #
p_age_bins <- ggplot(tab_agebin, aes(x = AgeBin, y = prop, fill = Type)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                position = position_dodge(width = 0.72), width = 0.18) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult=c(0,0.05))) +
  labs(x = "Insertion age bins (Myr)", y = "Fraction within ploidy",
       title = "Age-stratified composition of full-length LTRs by ploidy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 18, hjust = 1)) + 
  scale_fill_Publication('Set2')
print(p_age_bins)
ggsave('../figure/figure3/Cereba.new_insertion.CI.pdf',p_age_bins,width = 6,height = 4.2)

######################## 3) 连续年龄的分布图（小提琴+盒须）+ 整体/两两秩和检验（可选）
# 小提琴+盒须（连续 Age）
p_violin <- ggplot(all.1A, aes(Type, Age, fill = Type)) +
  geom_violin(trim = FALSE, alpha = 0.7, color = NA) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, alpha = 0.9) +
  labs(x = "", y = "Insertion age (Myr)",
       title = "Age distribution of full-length LTRs") +
  theme_classic() +
  theme(legend.position = "none") +
  geom_signif(comparisons = list(c("Diploid", "Tetraploid")), map_signif_level=TRUE)+
  geom_signif(comparisons = list(c("Hexaploid", "Tetraploid")), map_signif_level=TRUE)
print(p_violin)

# Kruskal–Wallis（整体）
kw <- kruskal.test(Age ~ Type, data = all.1A); kw

# Dunn 事后（如需）
# library(FSA)
# dunn <- FSA::dunnTest(Age ~ Type, data = all.1A, method = "bh"); dunn

############################## 4) “post-4x 年轻插入”主阈值（Age ≤ 0.5 Myr ≡ X1 ≥ 0.987）：Fisher 检验
thr_age <- 0.50  # Myr

# —— 各组命中/总数 —— #
tab_hit <- all.1A %>%
  mutate(hit = Age <= thr_age) %>%
  count(Type, hit, name="n") %>%
  pivot_wider(names_from = hit, values_from = n, values_fill = 0) %>%
  rename(no = `FALSE`, yes = `TRUE`) %>%
  mutate(total = yes + no, prop = yes/total) %>%
  rowwise() %>%
  mutate(ci = list(binom::binom.confint(yes, total, methods = "wilson")[,c("lower","upper")])) %>%
  ungroup() %>% tidyr::unnest_wider(ci, names_sep = "_")

tab_hit   # 每组：命中 yes、未命中 no、总数、比例、Wilson CI

# —— 成对 Fisher + OR（0.5 校正） —— #
pair_fisher <- function(df, g1, g2) {
  a <- df %>% filter(Type == g1) %>% select(yes, no) %>% as.numeric()
  b <- df %>% filter(Type == g2) %>% select(yes, no) %>% as.numeric()
  m <- rbind(a, b)
  ft <- fisher.test(m)  # 双侧
  OR <- ((m[1,1]+0.5)*(m[2,2]+0.5))/((m[1,2]+0.5)*(m[2,1]+0.5))  # HA 0.5 校正
  tibble(comp = paste(g2, "vs", g1),
         g1_yes = m[1,1], g1_no = m[1,2],
         g2_yes = m[2,1], g2_no = m[2,2],
         g1_prop = m[1,1]/sum(m[1,]),
         g2_prop = m[2,1]/sum(m[2,]),
         OR = OR,
         p_two_sided = ft$p.value)
}

res24 <- pair_fisher(tab_hit, "Diploid",   "Tetraploid")
res46 <- pair_fisher(tab_hit, "Tetraploid","Hexaploid")
res26 <- pair_fisher(tab_hit, "Diploid",   "Hexaploid")
bind_rows(res24, res46, res26)

# —— 3×2 总体 Fisher —— #
mat_all <- tab_hit %>% select(Type, yes, no) %>%
  column_to_rownames("Type") %>% as.matrix()
fisher.test(mat_all)



################## 5) 阈值敏感性扫描（展示稳健性）
scan_thr <- function(thrs = seq(0.45, 0.60, by=0.01)) {
  map_dfr(thrs, function(t0) {
    tmp <- all.1A %>%
      mutate(hit = Age <= t0) %>%
      count(Type, hit, name="n") %>%
      pivot_wider(names_from = hit, values_from = n, values_fill=0) %>%
      rename(no = `FALSE`, yes = `TRUE`)
    getm <- function(g) tmp %>% filter(Type==g) %>% select(yes,no) %>% as.matrix()
    m24 <- rbind(getm("Diploid"),    getm("Tetraploid"))
    m46 <- rbind(getm("Tetraploid"), getm("Hexaploid"))
    tibble(thr=t0,
           p_2x4x = fisher.test(m24)$p.value,
           p_4x6x = fisher.test(m46)$p.value)
  })
}
sens <- scan_thr()

# 可视化（可选）
p_sens <- sens %>%
  pivot_longer(cols = starts_with("p_"), names_to = "contrast", values_to = "p") %>%
  ggplot(aes(thr, p, color=contrast)) +
  geom_line(size=1) + geom_point() +
  geom_hline(yintercept=0.05, linetype=2) +
  labs(x="Age threshold (Myr)", y="Fisher p-value",
       title="Sensitivity of post-4x young enrichment to the age threshold") +
  theme_classic()
print(p_sens)
ggsave('../figure/Sf3/sensitive_of_age_threshold.line.pdf',width = 6,height = 4)

################## 6) Firth logistic（小样本更稳；可作为主分析的佐证）
df_bin <- all.1A %>%
  mutate(hit = Age <= thr_age,
         Type = factor(Type, levels=c("Diploid","Tetraploid","Hexaploid")))

fit <- logistf(hit ~ Type, data = df_bin)
summary(fit)         # OR（相对 Diploid）、Wald z、p 值、置信区间

# 若需要 Tetraploid 为参照来比较 6x：
df_bin2 <- df_bin %>% mutate(Type = relevel(Type, ref="Tetraploid"))
fit2 <- logistf(hit ~ Type, data = df_bin2)
summary(fit2)


############### 7) 连续年龄主图（替代你原来的 KDE，更直观）
# 已在上面 p_violin 给出；如需密度曲线保持一致带宽：
p_density <- ggplot(all.1A, aes(x=Age, color=Type)) +
  geom_density(adjust=1, size=1) +   # adjust 可统一带宽; 可调试
  labs(x="Insertion age (Myr)", y="Density",
       title="Kernel density of flLTR age by ploidy") +
  theme_classic()
print(p_density)


ggplot(all.1A, aes(x=Age, color=Type)) +
  stat_ecdf() +   # adjust 可统一带宽; 可调试
  labs(x="Insertion age (Myr)", y="Cumulative Fraction",
       title="Kernel density of flLTR age by ploidy") +
  theme_classic()



# 5. 2A and 3A cluster CRW ------------------------------------------------


