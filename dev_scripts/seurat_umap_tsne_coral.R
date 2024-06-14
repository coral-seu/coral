#' @title Visualization of Coral Cell Atlas
#' @description This tool is designed to carry out graphical rendering after 
#' obtaining t-SNE or UMAP results using the Seurat tool.
#' 
#' @param path should be the sample of umap.csv or tsne.csv
#' @param cluster should be the sample of cluster.csv
#' @param name refers to the possible cell types that can be identified before 
#' and after clustering analysis using the Seurat tool in scRNA-seq data analysis
#' @param outdir should be the target path for saving
#' @param projection choose to use UMAP or t-SNE as the input
#' @param sample should be the name of sample

library(argparse)
library(ggplot2)
library(tidyverse)
library(pals)
library(ggtext)

parser=ArgumentParser()
parser$add_argument("--path", help="the sample of umap.csv or tsne.csv")
parser$add_argument("--cluster", help="the sample of cluster.csv")
parser$add_argument("--name", help="the name of cluster")
parser$add_argument("--outdir", help="outdir of project")
parser$add_argument("--projection", help="umap or tsne")
parser$add_argument("--sample", help="the name of sample")
args <- parser$parse_args()
str(args)

path=args$path
cluster=args$cluster
name=args$name
projection=args$projection
sample=args$sample
outdir=args$outdir

define_theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(), 
  axis.title = element_text(color='black', size=15),
  axis.ticks.length = unit(0.4, "lines"),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text = element_blank(),
  legend.title = element_blank(),
  legend.text = element_text(size=18,face = 'bold'),
  legend.key = element_blank(),
  legend.key.size = unit(1, 'cm')
) 
define_color1 <- c("#4A8360","#A14F1E","#BD9826","#6376B8","#A4256D","#5D5595","#A1B9D1","#4F8730","#C46D26","#5A086B","#6376B8","#181311","#859C5E","#A12120","#A9CA78","#CD8483","#D9AE64","#B19DC4","#668F2E","#535353","#7F6324","#49277A","#804825","#363D6C","#B24F3C","#65A7AC","#CB89B1","#6981AD","#E9ED6F","#94694F","#DC050C","#FB8072","#1965B0","#7BAFDE","#882E72",
                   "#B17BA6","#FF7F00","#FDB462","#E7298A","#E78AC3",
                   "#33A02C","#B2DF8A","#55A1B1","#8DD3C7","#A6761D",
                   "#E6AB02","#7570B3","#BEAED4","#666666","#999999",
                   "#aa8282","#d4b7b7","#8600bf","#ba5ce3","#808000",
                   "#aeae5c","#1e90ff","#00bfff","#56ff0d","#ffff00")



data_2<-read_csv(cluster)
rename<-read_csv(name)
rename$newname<-as.character(rename$newname)

if(projection=="umap"){
  data_1<- read_csv(path)
  data<-merge.data.frame(data_1,data_2,by="Barcode")
  for(line in rename$oldname){
    clu<-match(line,rename$oldname)
    data$Lable[data$Cluster==line]<-as.character(rename[clu,2])
  }
  Lable<-as.factor(data$Lable)
  
  class_avg <- data %>%group_by(Lable) %>%summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))
  
  p<-ggplot(data, aes(x = UMAP_1, y = UMAP_2))+ geom_point(aes(color=Lable))+scale_color_manual(values = define_color1) +geom_label(aes(label = Lable),data=class_avg,nudge_x=0, label.r = unit(0.2, 'lines'),alpha = .5, fontface = 'bold',size = 6) +theme_bw() +define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$UMAP_1), 
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1) + 10, 
                     yend = min(data$UMAP_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$UMAP_1),
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1), 
                     yend = min(data$UMAP_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, file = paste0(outdir,"/",sample,"_Cluster_UMAP_1.png"), dpi = 300, width = 12, height = 8)
  ggsave(p,file=paste0(outdir,"/",sample,"_Cluster_UMAP_1.pdf"), width = 12, height = 8)
  
  
  ####彩色字白色边框
  p<-ggplot(data, aes(x = UMAP_1, y = UMAP_2))+ geom_point(aes(color=Lable))+ggrepel::geom_text_repel(data = class_avg,
                                                                                                      mapping = aes(label = Lable, col = Lable),bg.color = "white",bg.r = .25,
                                                                                                      size = 6,
                                                                                                      fontface = 'bold',force=0,
                                                                                                      show.legend = FALSE) +
    scale_color_manual(values = define_color1) +
    theme_bw() +define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$UMAP_1), 
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1) + 10, 
                     yend = min(data$UMAP_2)),
                 colour = "black", size = 1,arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$UMAP_1),
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1), 
                     yend = min(data$UMAP_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, filename = paste0(outdir,"/",sample,"_Cluster_UMAP_2.png"), dpi = 300, width = 12, height = 8)
  ggsave(file=paste0(outdir,"/",sample,"_Cluster_UMAP_2.pdf"),plot = p,  width = 12, height = 8)
  ####黑色字白色边框
  p<-ggplot(data, aes(x = UMAP_1, y = UMAP_2))+ geom_point(aes(color=Lable))+ggrepel::geom_text_repel(data = class_avg,
                                                                                                      mapping = aes(label = Lable),bg.color = "white",####文字边框
                                                                                                      bg.r = .25,
                                                                                                      size = 6,
                                                                                                      fontface = 'bold',force=0,
                                                                                                      show.legend = FALSE) +
    scale_color_manual(values = define_color1) +
    theme_bw() +
    define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$UMAP_1), 
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1) + 10, 
                     yend = min(data$UMAP_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$UMAP_1),
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1), 
                     yend = min(data$UMAP_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, file = paste0(outdir,"/",sample,"_Cluster_UMAP_3.png"), dpi = 300, width = 12, height = 8)
  ggsave(paste0(outdir,"/",sample,"_Cluster_UMAP_3.pdf"),plot = p,  width = 12, height = 8)
  ####无字版
  p<-ggplot(data, aes(x = UMAP_1, y = UMAP_2))+ geom_point(aes(color=Lable))+scale_color_manual(values = define_color1) +
    theme_bw() +
    define_theme +guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$UMAP_1), 
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1) + 10, 
                     yend = min(data$UMAP_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$UMAP_1),
                     y = min(data$UMAP_2),
                     xend = min(data$UMAP_1), 
                     yend = min(data$UMAP_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(plot = p, filename = paste0(outidr,"/",sample,"_Cluster_UMAP_4.png"), dpi = 300, width = 12, height = 8)
  ggsave(paste0(outdir,"/",sample,"_Cluster_UMAP_4.pdf"),plot = p,  width = 12, height = 8)
  
}

if(projection=="tsne"){ 
  data_1<- read_csv(path)
  data<-merge.data.frame(data_1,data_2,by="Barcode")
  for(line in rename$oldname){
    clu<-match(line,rename$oldname)
    data$Lable[data$Cluster==line]<-as.character(rename[clu,2])
  }
  Lable<-as.factor(data$Lable)
  
  class_avg <- data %>%group_by(Lable) %>%summarise(tSNE_1 = median(tSNE_1),tSNE_2 = median(tSNE_2))
  
  p<-ggplot(data, aes(x = tSNE_1, y = tSNE_2))+ geom_point(aes(color=Lable))+scale_color_manual(values = define_color1) +geom_label(aes(label = Lable),data=class_avg,nudge_x=0, label.r = unit(0.2, 'lines'),alpha = .5, fontface = 'bold',size = 6) +theme_bw() +define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$tSNE_1), 
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1) + 10, 
                     yend = min(data$tSNE_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$tSNE_1),
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1), 
                     yend = min(data$tSNE_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, file = paste0(outdir,"/",sample,"_Cluster_tSNE_1.png"), dpi = 300, width = 12, height = 8)
  ggsave(p,file=paste0(outdir,"/",sample,"_Cluster_tSNE_1.pdf"), width = 12, height = 8)
  
  
  ####彩色字白色边框
  p<-ggplot(data, aes(x = tSNE_1, y = tSNE_2))+ geom_point(aes(color=Lable))+ggrepel::geom_text_repel(data = class_avg,
                                                                                                      mapping = aes(label = Lable, col = Lable),bg.color = "white",bg.r = .25,
                                                                                                      size = 6,
                                                                                                      fontface = 'bold',force=0,
                                                                                                      show.legend = FALSE) +
    scale_color_manual(values = define_color1) +
    theme_bw() +define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$tSNE_1), 
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1) + 10, 
                     yend = min(data$tSNE_2)),
                 colour = "black", size = 1,arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$tSNE_1),
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1), 
                     yend = min(data$tSNE_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, filename = paste0(outdir,"/",sample,"_Cluster_tSNE_2.png"), dpi = 300, width = 12, height = 8)
  ggsave(file=paste0(outdir,"/",sample,"_Cluster_tSNE_2.pdf"),plot = p,  width = 12, height = 8)
  ####黑色字白色边框
  p<-ggplot(data, aes(x = tSNE_1, y = tSNE_2))+ geom_point(aes(color=Lable))+ggrepel::geom_text_repel(data = class_avg,
                                                                                                      mapping = aes(label = Lable),bg.color = "white",####文字边框
                                                                                                      bg.r = .25,
                                                                                                      size = 6,
                                                                                                      fontface = 'bold',force=0,
                                                                                                      show.legend = FALSE) +
    scale_color_manual(values = define_color1) +
    theme_bw() +
    define_theme +
    guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$tSNE_1), 
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1) + 10, 
                     yend = min(data$tSNE_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$tSNE_1),
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1), 
                     yend = min(data$tSNE_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(p, file = paste0(outdir,"/",sample,"_Cluster_tSNE_3.png"), dpi = 300, width = 12, height = 8)
  ggsave(paste0(outdir,"/",sample,"_Cluster_tSNE_3.pdf"),plot = p,  width = 12, height = 8)
  ####无字版
  p<-ggplot(data, aes(x = tSNE_1, y = tSNE_2))+ geom_point(aes(color=Lable))+scale_color_manual(values = define_color1) +
    theme_bw() +
    define_theme +guides(colour = guide_legend(override.aes = list(size = 5))) +
    geom_segment(aes(x = min(data$tSNE_1), 
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1) + 10, 
                     yend = min(data$tSNE_2)),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm"))) + 
    geom_segment(aes(x = min(data$tSNE_1),
                     y = min(data$tSNE_2),
                     xend = min(data$tSNE_1), 
                     yend = min(data$tSNE_2) + 10),
                 colour = "black", 
                 size = 1,
                 arrow = arrow(length = unit(0.3,"cm")))+
    theme(axis.title.x=element_text(hjust=0.05,vjust=8),axis.title.y=element_text(hjust=0.05,vjust=-8))
  ggsave(plot = p, filename = paste0(outdir,"/",sample,"_Cluster_tSNE_4.png"), dpi = 300, width = 12, height = 8)
  ggsave(paste0(outdir,"/",sample,"_Cluster_tSNE_4.pdf"),plot = p,  width = 12, height = 8)
}
