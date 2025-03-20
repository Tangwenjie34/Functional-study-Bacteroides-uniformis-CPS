####Data Processing and Visualization
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(readxl)
library(tidyverse)
library(plotrix)
library(ggh4x)
library(pheatmap)
library(ggplot2)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggprism)
library(ggbeeswarm)
library(patchwork)

a <- read.csv("./data.csv", row.names=1)
ggplot(a,aes(Cage.group,IgA))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=Cage.group),pch=21,stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  scale_fill_manual(values = c("grey","grey","grey","grey","grey",
                                     "#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
                                     "#e95f5c", "#e95f5c", "#e95f5c", "#e95f5c", "#e95f5c"))+
                                       scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 16,  
              base_line_size = 0.25, 
              axis_text_angle = 45)+ 
  facet_wrap( ~group, scales = 'free_x', ncol = 3)

group_names <- c( "Con","IgA-High",  "IgA-Low")
comparisons <- combn(group_names, 2, simplify = FALSE)
ggplot(a,aes(group,IgA))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),pch=21,stroke=0.1,alpha=1,size=3,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  scale_fill_manual(values = c( "grey","#5cc3e8","#e95f5c"))+
  theme(axis.text.x=element_text(angle=30,vjust=1, hjust=1))+
  theme_prism(base_fontface = "plain", 
              base_family = "serif", 
              base_size = 12,
              base_line_size = 0.25, 
              axis_text_angle = 0) 

a <- read.csv("./Elisa-data.csv", row.names=1)
data<-subset(a,a$IgA.Coating!="NA") 
group_names <- c( "Con","IgA-High",  "IgA-Low")
comparisons <- combn(group_names, 2, simplify = FALSE) 
ggplot(data,aes(group,IgA.Coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),stroke=0,alpha=1,size=3,
              width = 0.1, height = 0,pch=21)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c( "grey","#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_y_continuous(limits = c(0, 40), breaks = seq(5,30, 5))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", 
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.25, 
              axis_text_angle = 45) 

a <- read.csv("Elisa-data1.csv", row.names=1)
data<-subset(a,a$IgA.cohusing!="NA") 
group_names <- c( "IgA-High-Initial","IgA-Low-Initial","IgA-High-post","IgA-Low-post")
data$group.cohousing <- factor(data$group.cohousing,levels = group_names)
ggplot(data,aes(group.cohousing,IgA.cohusing,fill=group.cohousing))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ 
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.5, 
              axis_text_angle = 45)+
  scale_fill_manual(values = c("#dedede", "grey","#dedede","grey"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  geom_signif(comparisons = list(c("IgA-Low-Initial","IgA-Low-post"),
                                 c("IgA-High-post","IgA-Low-post")),
              map_signif_level = T,
              test = "wilcox.test", 
              y_position = c(500,520,550,570),
              size=0.8,color="black")+
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 100)) +
  facet_wrap( ~regroup.cohousing, scales = 'free_x')

a <- read.csv("./SI-Colon-Elisa-data.csv", row.names=1)
group_names <- c( "IgA-High","IgA-Low")
a$GI.tract <- factor(a$GI.tract,levels =c("Pro-SI","Distal-SI","Cecum","Pro-colon","Distal-colon"))
ggplot(a,aes(group,IgA,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+
  theme_classic()+ 
  scale_fill_manual(values = c("#dedede", "grey","#dedede","grey"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  geom_signif(comparisons = list(c("IgA-Low","IgA-High")),
              map_signif_level = T,
              test = "wilcox.test", 
              y_position = c(450),
              size=0.25,color="black")+
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  facet_wrap( ~GI.tract, scales = 'free_x',ncol = 5)

a <- read.csv("./Elisa-allbacterialtest.csv", row.names=1)
group_names <- c( "Con","B.uniformis",  "B.ovatus", "B.xylanisolvens","B.faecis", "B.intestinalis", "B.caccae",
                  "B.vulgatus", "P.copri", "P.intermedia", "P.mirabilis", "F.varium")
comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA))+
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),pch=21,alpha=1,size=2.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  scale_fill_manual(values = c(
    "grey","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
          "#e95f5c","#e95f5c","#ffdb00","#79ceb8"))+
     theme_classic()+ scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 100))

a <- read.csv("Iga-coating.csv", row.names=1)
group_names <- c( "Con","B.uniformis",  "B.ovatus", "B.xylanisolvens","B.faecis", "B.intestinalis", "B.caccae",
                  "B.vulgatus", "P.copri", "P.intermedia", "P.mirabilis", "F.varium")
a$group <- factor( a$group, level=group_names )
comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
ggplot(a,aes(group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),shape=21,stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
    "grey","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
          "#e95f5c","#e95f5c","#ffdb00","#79ceb8"))+
            theme_classic()+
  scale_y_continuous(limits = c(0, 40), breaks = seq(0,30, 5))

ggplot(data = a ,mapping = aes( x = IgA, y= IgA.Coating)) + 
  geom_point(size= 2) +
  stat_smooth(method = "lm")+ 
  ggpubr::stat_cor(aes(), label.x =200, label.y = 28)+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30, 4))+
  scale_x_continuous(limits = c(200, 480), breaks = seq(200,480, 100))

c <- read.csv("./Elisa-data3.csv", row.names=1)
c <- subset(c, c$Group1 %in% c("Con_1","PEC+BU_1" ,  "NOPEC+BU_1"))
group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU")
c$Group <- factor(c$Group, level=group_names )
comparisons <- combn(group_names, 2, simplify = FALSE)

ggplot(c,aes(Group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=Group),shape=21,stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
    "grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_y_continuous(limits = c(0, 40), breaks = seq(0,30, 5))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain",
              base_family = "serif", 
              base_size = 16, 
              base_line_size = 0.25, 
              axis_text_angle = 45)

b <- read.csv("./Elisa-data5.csv", row.names=1)
group_names <- c( "CON",  "PEC", "Starch","XYL","Inulin")
comparisons <- lapply(group_names[-12], function(group) c(group, "CON"))
b$Group <- factor(b$Group, level=group_names )

topbar <- function(x){
  return(mean(x)+sd(x)/sqrt(length(x)))
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}

ggplot(b, aes(week, IgA ,group = Group, colour = Group)) +
  geom_point() +
  geom_smooth(method = 'loess')+
  facet_wrap( ~Group, scales = 'free_x', ncol = 5)+
  scale_y_continuous(limits = c(50, 400), breaks = seq(0, 400, 80))+
  scale_color_manual(values = c("grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_fill_manual(values = c("grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  theme_classic(base_size = 16)+
  theme(legend.position = 'none')+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.62) +
  scale_x_continuous(limits = c(0, 5), breaks = seq(0, 5, 1))

c <- read.csv("./four-weeks-Elisa.csv", row.names=1)
group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU"  ,     "XYL+BU"   ,  "NOXYL+BU" , 
                  "Inulin+BU" ,   "NO_Inulin+BU" )
c$Group <- factor(c$Group, level=group_names )
topbar <- function(x){
  return(mean(x)+sd(x)/sqrt(length(x))) 
}
bottombar <- function(x){
  return(mean(x)-sd(x)/sqrt(length(x)))
}
ggplot(c,aes(week,IgA,color=Group))+
  stat_summary(geom = 'line',fun='mean',cex=1)+
  stat_summary(geom = 'errorbar',
               fun.min = bottombar,fun.max = topbar,
               width=0.12,cex=0.8,aes(color=Group))+
  stat_summary(geom = 'point',fun='mean',aes(fill=Group),
               size=2.5,pch=21,color='black')+
  scale_y_continuous(limits = c(230, 400), breaks = seq(0, 400, 50))+
  scale_x_continuous(limits = c(1, 7.5), breaks = seq(1, 7, 1)) +
  scale_color_manual(values = c("grey", "#5cc3e8","#3b8000","#0151b8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_fill_manual(values = c("grey", "#5cc3e8", "#3b8000","#0151b8","#e95f5c","#ffdb00","#79ceb8"))+  
  theme() +
  annotate(geom = 'text',label='***',x=6.5,y=342,size=4,angle=90)+
  annotate(geom = 'text',label='****',x=6.8,y=312,size=4,angle=90)
