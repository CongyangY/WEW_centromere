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







# 1. coverage -------------------------------------------------------------
hifi.cov <- read_delim('../data/QC/hifi.500k.bg',delim = '\t',col_names = F) %>% 
  filter(grepl("chr", X1))
ggplot(hifi.cov) + geom_line(aes(x = (X2+X3)/2,y = X4)) +  
  facet_wrap(~ X1,nrow = 7) + theme_Publication() + 
  ylim(c(0,0.4))


ont.cov <- read_delim('../data/QC/ont.500k.bg',delim = '\t',col_names = F) %>% 
  filter(grepl("chr", X1))
ggplot(ont.cov) + geom_line(aes(x = (X2+X3)/2,y = X4)) +  
  facet_wrap(~ X1,nrow = 7) + theme_Publication() + 
  ylim(c(0,0.2))


# 2. RNAfold --------------------------------------------------------------

RNAfold <- read_delim('../data/QC/RNAfold/WEW_new.RNAfold2.500k.cent30Mb.bg',delim = '\t',
                      col_names = F)

ggplot(RNAfold) +
  geom_line(aes(x = (X2 + X3) / 2, y = X4)) +
  facet_wrap(~X1, ncol = 2) + 
  theme_Publication() + ylim(c(-35,-25))


# centromere size ---------------------------------------------------------

cent.size <- data.frame(
  sub = c(rep('A',7),rep('B',7)),
  chr = c(1:7,1:7),
  size = c(7.2,7.5,7.3,9.23,7.2,8.4,8.1,8,9.2,6.9,10,8.5,7.5,7.9),
  len = c(611.773,735.254,790.742,842.232,766.662,872.156,751.531,696.108,722.677,738.031,634.596,759.896,752.998,773.843)
) %>% mutate(id = paste('chr',chr,sub,sep = ''))

ggplot(cent.size,aes(x = sub,y = size)) + geom_boxplot() + 
  geom_signif(comparisons = list(c("A", "B")), test =  "wilcox.test",test.args = list(paired = TRUE,alternative = 'less')) + 
  theme_Publication() 
ggsave('../figure/Sf4/cent.size.box.pdf',width = 2.5,height = 3.5)


library(ggplot2)
library(ggpubr)
ggplot(cent.size, aes(x = size, y = len)) +
  geom_point() +
  theme_Publication()
ggsave('../figure/Sf4/cent.size.point.cor.pdf',width = 3.4,height = 3.5)


# Langdon slop 1mb --------------------------------------------------------
Langdon.cent_core <- read_delim('../data/fl-LTR/Langdon/Langdon.cent.flLTR.insertTime.bg',
                           delim = '\t',col_names = F) %>% 
  mutate(Type = 'Cent')
Langdon.cent_flank <- read_delim('../data/fl-LTR/Langdon/Langdon.cent.flank_1mb.flLTR.insertTime.bg',
                           delim = '\t',col_names = F) %>% 
  mutate(Type = 'Flank')

Langdon.cent <- rbind(Langdon.cent_flank,Langdon.cent_core)

ggplot(Langdon.cent,aes(x = Type, y = X4,fill = Type)) + geom_boxplot(outlier.alpha = 0) + 
  theme_Publication() + scale_fill_Publication() + 
  geom_signif(comparisons = list(c("Cent", "Flank")),y_position = 5.3)
ggsave('../figure/Sf7/Langdon.cent.LTR_time.box.pdf',width = 2.5,height = 3.5)
  

ggplot(Langdon.cent) + stat_ecdf(aes(x = X4, color = Type)) + 
  theme_Publication() + scale_colour_Publication()
ggsave('../figure/Sf7/Langdon.cent.LTR_time.ecdf.pdf',width = 2.5,height = 3.5)




# fish imageJ -------------------------------------------------------------

chr2A.fish <- read_delim(pipe("pbpaste"),delim = '\t')
write_csv(chr2A.fish,'../data/Cereba_new/Cereba.chr2A.fish.csv')

chr1A.fish <- read_delim(pipe("pbpaste"),delim = '\t')
write_csv(chr1A.fish,'../data/Cereba_new/Cereba.chr1A.fish.csv')

chr3A.fish <- read_delim(pipe("pbpaste"),delim = '\t')
write_csv(chr3A.fish,'../data/Cereba_new/Cereba.chr3A.fish.csv')


chr1A.fish <- read_delim('../data/Cereba_new/Cereba.chr1A.fish.csv',delim = ',',
                         col_names = T)
colnames(chr1A.fish) <- c('Distance','Value')
ggplot(chr1A.fish) + geom_line(aes(x = Distance,y = Value)) + theme_Publication() 
ggsave('../figure/Sf10/chr1A.fish.line.pdf',width = 2.5,height = 1.8)

colnames(chr2A.fish) <- c('Distance','Value')
ggplot(chr2A.fish) + geom_line(aes(x = Distance,y = Value)) + theme_Publication() 
ggsave('../figure/Sf10/chr2A.fish.line.pdf',width = 2.5,height = 1.8)


colnames(chr3A.fish) <- c('Distance','Value')
ggplot(chr3A.fish) + geom_line(aes(x = Distance,y = Value)) + theme_Publication() 
ggsave('../figure/Sf10/chr3A.fish.line.pdf',width = 2.5,height = 1.8)



# 6 cent566 similarity ----------------------------------------------------
cent566.zavitan <- read_delim('../data/CENH3/cent566/566_in_WEWnew.slop10M.m8.similar',
                              delim = '\t',col_names = F) %>% mutate(Type = 'Zavitan')
cent566.CS <- read_delim('../data/CENH3/cent566/566_in_CSGL.slop10M.m8.similar',
                              delim = '\t',col_names = F) %>% mutate(Type = 'CS')
cent566 <- rbind(cent566.zavitan,cent566.CS)

ggplot(cent566) + geom_density(aes(x = X1, color = Type)) + theme_Publication() + 
  scale_colour_Publication()
ggsave('../figure/Sf12/cent566.metaplot.pdf',height = 3,width = 3)
wilcox.test(cent566.zavitan$X1,cent566.CS$X1)
# p-value < 2.2e-16


