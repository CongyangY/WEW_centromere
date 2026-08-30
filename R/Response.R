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








# 1. CENH3 signal in contigs ----------------------------------------------

contigs.cenh3 <- read_delim('../2026/data/centromere/zavitan.contigs.CENH3.q20.bg',
                            delim = '\t',col_names = F) %>% mutate(Type = 'Contigs')
cent.cenh3 <- read_delim('../2026/data/centromere/zavitan.cent.CENH3.q20.bg',
                            delim = '\t',col_names = F) %>% mutate(Type = 'Cent')
arm_random.cenh3 <- read_delim('../2026/data/centromere/zavitan.arm_random.CENH3.q20.bg',
                            delim = '\t',col_names = F)%>% mutate(Type = 'Random')

cen3.compare <- rbind(contigs.cenh3,cent.cenh3,arm_random.cenh3)
cen3.compare$Type <- factor(cen3.compare$Type, levels = c('Cent', 'Random', 'Contigs'))
ggplot(cen3.compare,aes(x = Type, y = log10(X4),fill = Type)) + geom_boxplot() + 
  theme_Publication() + scale_fill_Publication() + 
  geom_signif(comparisons = list(c("Cent", "Random")), y_position = 1.4,
              map_signif_level=TRUE) + 
  geom_signif(comparisons = list(c("Random", "Contigs")), y_position = 1.4,
              map_signif_level=TRUE) + ylim(-2,2) + 
  theme(
    legend.position = 'None'
  ) + labs(x = 'Region',y = 'Log10(relative signal of CENH3)')
ggsave('../2026/figure/Sf2/CENH3.signal.compare.boxplot.pdf',width = 3,height = 4)


# 2 unscaffold region stat-----------------------------------------------------------------------
all_stat <- read_delim('../2026/data/Response/syn/best_All_type.tsv',delim = '\t',col_names = F)

# 重命名
colnames(all_stat) <- c(
  "Region_ID", "Chr", "Start", "End", "Length_bp",
  "TE_Type", "BestOv", "TotalOv", "Frac"
)
all_stat2 <- all_stat %>%
  mutate(
    Subgenome = case_when(
      str_detect(Chr, "A$") ~ "A",
      str_detect(Chr, "B$") ~ "B",
      TRUE ~ "Other"
    )
  ) %>%
  filter(Subgenome %in% c("A", "B"))

# 如果类型太多，可以保留前10类，其余合并为 Other
top_types <- all_stat2 %>%
  group_by(TE_Type) %>%
  summarise(TotalLength = sum(Length_bp, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalLength)) %>%
  slice(1:10) %>%
  pull(TE_Type)

all_stat2 <- all_stat2 %>%
  mutate(
    TE_Type2 = ifelse(TE_Type %in% top_types, TE_Type, "Other")
  )

# 统计 count
stat_count <- all_stat2 %>%
  count(Subgenome, TE_Type2, name = "Value") %>%
  mutate(Metric = "Count")

# 统计 length
stat_length <- all_stat2 %>%
  group_by(Subgenome, TE_Type2) %>%
  summarise(Value = sum(Length_bp, na.rm = TRUE), .groups = "drop") %>%
  mutate(Metric = "Length")

# 合并
plot_df <- bind_rows(stat_count, stat_length)

# 统一排序：按总长度排序
type_order <- all_stat2 %>%
  group_by(TE_Type2) %>%
  summarise(TotalLength = sum(Length_bp, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(TotalLength)) %>%
  pull(TE_Type2)

plot_df <- plot_df %>%
  mutate(
    TE_Type2 = factor(TE_Type2, levels = type_order),
    Metric = factor(Metric, levels = c("Count", "Length"))
  )

p <- ggplot(plot_df, aes(x = TE_Type2, y = Value, fill = Subgenome)) +
  geom_col(position = "dodge") +
  facet_wrap(~Metric, scales = "free_y", ncol = 1) +
  theme_bw() +
  labs(
    title = "Counts/length of dominant repeat types in previously unresolved regions",
    x = "Repeat type",
    y = NULL
  ) + scale_fill_Publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")
  ) 

print(p)

ggsave("../2026/figure/Sf_syn/TE_type_count_length_AB.pdf", p, width = 10, height = 8)
ggsave("../2026/figure/Sf_syn/TE_type_count_length_AB.png", p, width = 10, height = 8, dpi = 300)



# chromosome density ------------------------------------------------------
library(CMplot)

# 你的染色体长度
chr_len <- data.frame(
  Chromosome = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B",
                 "chr4A","chr4B","chr5A","chr5B","chr6A","chr6B",
                 "chr7A","chr7B"),
  Position = c(611773303,735254030,790741867,842231976,766661948,872156323,
               751530861,696107520,722676616,738030801,634596046,759895744,
               752998143,773842621)
)

# 给每条染色体添加一个末端“伪 SNP”
chr_len$SNP <- paste0(chr_len$Chromosome, "_end")

# 调整列顺序，确保前3列是 SNP / Chromosome / Position
chr_len <- chr_len[, c("SNP", "Chromosome", "Position")]
all_stat3 <- all_stat %>% mutate(
  Position = (Start + End)/2,
  SNP = Region_ID,
  Chromosome = Chr
) %>% select(SNP,Chromosome,Position)
# 假设原始数据 dat 也只有前三列，或者你只取前三列作密度图
dat2 <- rbind(all_stat3, chr_len)

# 染色体顺序建议固定，不然容易乱
dat2$Chromosome <- factor(
  dat2$Chromosome,
  levels = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B",
             "chr4A","chr4B","chr5A","chr5B","chr6A","chr6B",
             "chr7A","chr7B")
)

dat2 <- dat2[order(dat2$Chromosome, dat2$Position), ]

CMplot(
  dat2,
  plot.type = "d",
  bin.size = 1e6,
  chr.pos.max = TRUE,
  chr.den.col = c("darkgreen", "yellow", "red"),
  file = "pdf",
  file.name = "wheat_density",
  dpi = 300,
  file.output = TRUE,
  verbose = TRUE,
  width = 12,
  height = 6
)
file.rename('./Marker_Density.wheat_density.pdf','../2026/figure/Sf_syn/unscaffold_distribution.cmplot.pdf')





# centromere size ---------------------------------------------------------

chr_order <- c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B",
               "chr4A","chr4B","chr5A","chr5B","chr6A","chr6B",
               "chr7A","chr7B")

zavitan_cent <- data.frame(
  Chr = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B","chr4A","chr4B",
          "chr5A","chr5B","chr6A","chr6B","chr7A","chr7B"),
  Cent_start = c(224022000,259542000,376728000,374156000,324530000,357246000,
                 309020000,285180000,242936000,218638000,290828000,324594000,
                 360664000,317216000),
  Cent_end = c(230984000,265824000,383684000,382080000,331610000,365226000,
               317410000,292732000,252080000,225106000,299964000,332962000,
               367956000,324834000),
  Chr_len = c(611773303,735254030,790741867,842231976,766661948,872156323,
              751530861,696107520,722676616,738030801,634596046,759895744,
              752998143,773842621),
  Genome = "Zavitan"
)

langdon_cent <- data.frame(
  Chr = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B","chr4A","chr4B",
          "chr5A","chr5B","chr6A","chr6B","chr7A","chr7B"),
  Cent_start = c(214204000,260794000,347258000,371326000,326956000,361898000,
                 289512000,282202000,247902000,230736000,291122000,334914000,
                 370862000,310448000),
  Cent_end = c(221564000,266838000,354012000,378250000,330956000,369094000,
               295948000,289688000,254700000,236700000,299590000,342432000,
               377232000,316448000),
  Chr_len = c(602454030,712931786,803473472,823385595,763004926,863844817,
              767942507,697541477,722581089,738974146,631877482,744936831,
              754097208,751300603),
  Genome = "Langdon"
)

cent <- bind_rows(zavitan_cent, langdon_cent)

# 2. 数据整理
chr_y <- data.frame(
  Chr = factor(chr_order, levels = chr_order),
  base_y = rev(seq_along(chr_order)) * 1.15
)

cent <- cent %>%
  mutate(
    Chr = factor(Chr, levels = chr_order),
    Genome = factor(Genome, levels = c("Zavitan", "Langdon"))
  ) %>%
  left_join(chr_y, by = "Chr") %>%
  mutate(
    y = ifelse(Genome == "Zavitan", base_y + 0.18, base_y - 0.18),
    x_start = Cent_start / Chr_len * 100,
    x_end = Cent_end / Chr_len * 100
  )

backbone <- cent %>%
  select(Genome, Chr, Chr_len, y) %>%
  distinct() %>%
  mutate(
    x0 = 0,
    x1 = 100
  )

chr_labels <- cent %>%
  group_by(Chr, base_y) %>%
  summarise(y = mean(y), .groups = "drop")

genome_labels <- cent %>%
  mutate(label = as.character(Genome))
# 3. 作图

p <- ggplot() +
  # 染色体骨架
  geom_segment(
    data = backbone,
    aes(x = x0, xend = x1, y = y, yend = y),
    linewidth = 3,
    color = "grey82",
    lineend = "round"
  ) +
  # 着丝粒区间
  geom_rect(
    data = cent,
    aes(
      xmin = x_start, xmax = x_end,
      ymin = y - 0.11, ymax = y + 0.11,
      fill = Genome
    ),
    color = "black",
    linewidth = 0.25
  ) +
  # 左侧染色体标签
  geom_text(
    data = chr_labels,
    aes(x = -6, y = y, label = Chr),
    hjust = 1,
    size = 3.2
  ) +
  # 右侧 genome 标签
  geom_text(
    data = genome_labels,
    aes(x = 103, y = y, label = label),
    hjust = 0,
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(-10, 118),
    breaks = c(0, 25, 50, 75, 100),
    labels = c("0", "25", "50", "75", "100"),
    expand = c(0, 0)
  ) +
  labs(
    x = "Relative chromosome position (%)",
    y = NULL
  ) +
  scale_fill_manual(values = c("Zavitan" = "#4C78A8", "Langdon" = "#E45756")) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    axis.title.x = element_text(size = 10),
    axis.text.x = element_text(size = 9, color = "black"),
    plot.margin = margin(8, 25, 8, 25)
  )

# 显示图
print(p)

ggsave("../2026/figure/Figure3/Fig3A_paired_centromere_loc.pdf", p, width = 5, height = 4)




# centemre compare --------------------------------------------------------


cent <- data.frame(
  Chr = c("chr1A","chr1B","chr2A","chr2B","chr3A","chr3B","chr4A","chr4B",
          "chr5A","chr5B","chr6A","chr6B","chr7A","chr7B"),
  Langdon_rep1 = c(7.4,6.0,6.8,6.9,4.0,7.2,6.4,7.5,6.8,6.0,8.5,7.5,6.4,6.0),
  Langdon_rep2 = c(8.0,6.8,7.2,7.5,7.2,7.9,7.5,8.0,8.3,7.0,9.1,9.0,7.2,8.5),
  Langdon_rep3 = c(8.0,7.0,7.0,7.5,7.2,8.1,7.4,8.1,8.2,7.1,9.1,9.0,7.4,8.2),
  Zavitan_rep1 = c(7.0,6.3,7.0,7.9,7.1,8.0,8.4,7.6,9.1,6.5,9.1,8.4,7.3,7.6),
  Zavitan_rep2 = c(7.2,7.4,7.1,8.7,7.2,8.5,8.7,7.7,9.4,7.2,9.4,8.3,7.6,7.0)
)

cent_mean <- cent %>%
  mutate(
    Langdon = rowMeans(select(., starts_with("Langdon"))),
    Zavitan = rowMeans(select(., starts_with("Zavitan"))),
    Subgenome = ifelse(str_detect(Chr, "A$"), "A", "B")
  ) %>%
  select(Chr, Subgenome, Zavitan, Langdon)

cent_long <- cent_mean %>%
  pivot_longer(cols = c(Zavitan, Langdon),
               names_to = "Accession",
               values_to = "Centromere_size")

p1 <- ggplot(cent_long, aes(x = Accession, y = Centromere_size, group = Chr)) +
  geom_line(color = "grey60", linewidth = 0.5) +
  geom_point(aes(fill = Subgenome), shape = 21, size = 3, color = "black", stroke = 0.3) +
  theme_Publication() +
  labs(x = NULL, y = "Functional centromere span (Mb)") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(size = 11)
  )

ggsave("../2026/figure/Figure3/Fig3B_paired_centromere_size.pdf", p1, width = 2.5, height = 4)


cent_delta <- cent_mean %>%
  mutate(
    Delta = Langdon - Zavitan,
    Chr = factor(Chr, levels = Chr)
  )

p2 <- ggplot(cent_delta, aes(x = Chr, y = Delta, fill = Subgenome)) +
  geom_hline(yintercept = 0, linewidth = 0.4, linetype = "dashed") +
  geom_col(width = 0.65, color = "black", linewidth = 0.2) +
  theme_Publication() +
  labs(x = NULL, y = "Δ span (Langdon − Zavitan, Mb)") +
  theme(
    legend.title = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 9)
  )

ggsave("../2026/figure/Figure3/Fig3B_delta_centromere_size.pdf", p2, width = 2.5, height = 4)


t.test(cent_mean$Langdon, cent_mean$Zavitan, paired = TRUE) # p-value = 0.04677
wilcox.test(cent_mean$Langdon, cent_mean$Zavitan, paired = TRUE) # p-value = 0.05798



# SV ----------------------------------------------------------------------

library(tidyverse)
library(patchwork)

# 1. Input data

# Summary of all inversions >100 kb
inv_context <- data.frame(
  Category = c("Within single contig", "Across contigs/gaps"),
  Count = c(51, 8)
) %>%
  mutate(
    Category = factor(Category, levels = c("Within single contig", "Across contigs/gaps")),
    Percent = Count / sum(Count) * 100,
    Label = paste0(Count, " (", sprintf("%.1f", Percent), "%)")
  )

# Breakpoint-spanning read support for the 8 inversions crossing contig joins or gaps
inv_support <- data.frame(
  Event = c("INV9_chr1B", "INV58_chr2B", "INV25_chr4B", "INV10_chr5B",
            "INV4_chr6B", "INV5_chr6B", "INV10_chr6B", "INV40_chr7B"),
  Chr = c("chr1B", "chr2B", "chr4B", "chr5B", "chr6B", "chr6B", "chr6B", "chr7B"),
  Size_Mb = c(14.172404, 12.616582, 12.971629, 46.768496,
              17.976417, 17.673433, 13.649203, 9.331101),
  left_ONT_ratio = c(0.8333, 0.8182, 0.8462, 0.7742, 0.8276, 0.6667, 0.9062, 0.9565),
  right_ONT_ratio = c(0.9000, 0.9032, 0.7879, 0.8571, 0.7586, 0.9091, 0.7619, 0.9000),
  left_HiFi_ratio = c(0.7879, 0.7812, 0.7143, 0.7536, 0.6875, 0.5676, 0.6596, 0.6977),
  right_HiFi_ratio = c(0.6905, 0.7308, 0.7308, 0.5914, 0.7541, 0.6702, 0.7030, 0.7243)
) %>%
  mutate(
    ONT_min_ratio = pmin(left_ONT_ratio, right_ONT_ratio),
    HiFi_min_ratio = pmin(left_HiFi_ratio, right_HiFi_ratio),
    ONT_min_percent = ONT_min_ratio * 100,
    HiFi_min_percent = HiFi_min_ratio * 100
  )

# Save processed data
write.table(inv_support, "../2026/figure/Figure3/Fig3E_inversion_validation_data.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Convert to long format for plotting
inv_support_long <- inv_support %>%
  select(Event, Chr, Size_Mb, ONT_min_percent, HiFi_min_percent) %>%
  pivot_longer(
    cols = c(ONT_min_percent, HiFi_min_percent),
    names_to = "Data_type",
    values_to = "Support_percent"
  ) %>%
  mutate(
    Data_type = recode(Data_type,
                       "ONT_min_percent" = "ONT",
                       "HiFi_min_percent" = "HiFi"),
    Event = factor(Event, levels = rev(inv_support$Event))
  )

# 2. Plot style

base_theme <- theme_classic(base_size = 10) +
  theme(
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    plot.margin = margin(5, 5, 5, 5)
  )

# 3. Left panel: assembly context
p_context <- ggplot(inv_context, aes(x = Category, y = Count, fill = Category)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.25) +
  geom_text(
    aes(label = Label),
    hjust = -0.05,
    size = 3.2
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 60),
    breaks = seq(0, 60, 20),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    title = "Assembly context",
    x = NULL,
    y = "Number of inversions"
  ) +
  base_theme +
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 9, color = "black"),
    axis.ticks.y = element_blank()
  )

# 4. Right panel: conservative breakpoint-spanning support

p_support <- ggplot(inv_support_long,
                    aes(x = Support_percent, y = Event, shape = Data_type)) +
  geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.35) +
  geom_point(size = 2.7, stroke = 0.8) +
  scale_x_continuous(limits = c(50, 100), breaks = seq(50, 100, 10)) +
  labs(
    title = "Breakpoint-spanning support",
    x = "Conservative support ratio (%)",
    y = NULL
  ) +
  base_theme +
  theme(
    legend.position = "bottom"
  )

# 5. Combine and save

p_final <- p_context + p_support +
  plot_layout(widths = c(1.05, 1.55), guides = "collect") &
  theme(legend.position = "bottom")

ggsave("../2026/figure/Figure3/Fig3E_inversion_validation_left.pdf", p_context,
       width = 5, height = 1.8, useDingbats = FALSE)
ggsave("../2026/figure/Figure3/Fig3E_inversion_validation_right.pdf", p_support,
       width = 3.5, height = 5, useDingbats = FALSE)

