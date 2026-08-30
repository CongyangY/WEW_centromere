library(tidyverse)
library(ggsignif)

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






# A_AABB ------------------------------------------------------------------

invade_A <- read_delim('../data/fl-LTR/polyploid/A_AB.invade.flTR_in_WEW.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'A_invade')
all_A <- read_delim('../data/fl-LTR/polyploid/Lachesis_Tu.cent.fa.finder.combine.bed',delim = '\t',col_names = F) %>% 
  mutate(Type = 'A')

A_flLTR <- rbind(invade_A,all_A)
ggplot(A_flLTR,aes(x = Type,y = X5,fill = Type)) + geom_boxplot()+
  theme_Publication() + scale_fill_Publication() + ylim(c(0,8.4e7)) + 
  geom_signif(comparisons=list(c("A", "A_invade")),map_signif_level = T,size = 0.9,y_position = 7.9e7) + 
  labs(x = 'Type', y = 'Insertion time (Myr)')
ggsave('../figure/figure2/A_B_insert_time.boxplot.pdf',width =3.6 ,height = 3.8)

ggplot(A_flLTR,aes(x = Type,y = X5,fill = Type)) + geom_boxplot()+
  theme_Publication() + scale_fill_Publication() + ylim(c(0,1.2e7)) + 
  geom_signif(comparisons=list(c("A", "A_invade")),map_signif_level = T,size = 0.9,y_position = 1e7) + 
  labs(x = 'Type', y = 'Insertion time (Myr)')
ggsave('../figure/figure2/A_B_insert_time.enlarged.boxplot.pdf',width =3.6 ,height = 3.8)
