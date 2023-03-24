#rm(list = ls(all.names = TRUE))
library(ggplot2)
library(dplyr)
#arg <- commandArgs(trailingOnly = TRUE);
#file_name = arg[1]
#gene_name = arg[2]
setwd("E:/CisCross/R_graph/light")
file_name = "EIN3.csv"
gene_name = "EIN3 AT3G20770"
#file_name = "AT1G46480.csv"
#df <- rbind(arr10, ebs_2, ebs_1, ebs_gc, ebs_teil)
#df <- read.table( file_name, stringsAsFactors = F, header = T)
df <- read.csv(file_name, header = F, sep = "\t", stringsAsFactors = F)

len = length(df[,2])
TF <- df[1:len,1] 
rel_beg <- df[1:len,2]
rel_end <- df[1:len,3]

input_gene <- df
input_gene%>% 
  mutate(height = 1+(1:nrow(df)))%>%
  ggplot(aes(y = height))+#, group = TF, color = TF, fill = TF))+
  geom_rect(mapping=aes(xmin=rel_beg, xmax=rel_end, ymin=height, ymax=height-0.2), alpha=0.5)+
  
  geom_text(aes(rel_beg-120, height-0.08, label = rel_beg), size = 2.5, show.legend = FALSE)+
  geom_text(aes(rel_end+100, height-0.08, label = rel_end), size = 2.5, show.legend = FALSE)+
  geom_text(aes((rel_end+rel_beg)/2, height+0.4, label = TF), color = 'black', 
            size = 3, show.legend = FALSE)+
  scale_x_continuous(name = 'promoter', limits = c(-2000,550), breaks = c(seq(-2000, 400, 300), 1))+
  scale_y_continuous(limits = c(1,(nrow(df)+1+0.7)), breaks = seq(1,(nrow(df)+1),1), expand = c(0,0))+
  theme_classic()+
  ggtitle(gene_name)+
  theme(plot.title = element_text(hjust = 0.5, face = 'italic'))+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  geom_segment(aes(x = 1, y = 0.8, xend = 1, yend = 1.8), color = 'black')+
  geom_segment(aes(x = 1, y = 1.8, xend = 100, yend = 1.8), arrow = arrow(length = unit(0.02,units = "npc")),
               size = 0.8, color = 'black')+
  geom_text(aes(-50, 1.8, label = 'TSS'), color = 'black', size = 3)

ggsave(paste0(gene_name, '.png'), device = 'png',dpi = 300)#, height = 2100, width = 5600, units = 'px')

