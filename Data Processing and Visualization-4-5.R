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

a <- read.csv("./Elisa-ABUNDANCE.csv", row.names=1)
ggplot(data = a ,mapping = aes( x = IgA, y= logCFUg)) + 
  geom_point(size= 2) +
  stat_smooth(method = "lm")+ 
  ggpubr::stat_cor(aes(), label.x =20, label.y = 9)+
  scale_y_continuous(limits = c(0, 10), breaks = seq(0,10, 2))+
  scale_x_continuous(limits = c(20, 220), breaks = seq(20,220, 50))+theme_classic2()

a <- read.csv("CFU.csv", row.names=1)
a$position <- factor(a$position ,levels =c("Lumen","Epithelium","Tissue"))
a$group <- factor(a$group ,levels =c("WT","rcsC"))
p1 <-ggplot(a,aes(group,logCFUg,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ 
  theme_classic()+ 
  scale_fill_manual(values = c( "#5f5f60","#1754a5"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=3.5,
              width = 0.1, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  facet_wrap( ~position, scales = 'free_x',ncol = 5)
p2 <-ggplot(a,aes(group,Ration,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ 
  theme_classic()+ 
  scale_fill_manual(values = c( "#5f5f60","#1754a5"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=5,
              width = 0.1, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.2, 0.2))
combined_plot <- p1 + p2
print(combined_plot)

a <- read.csv("./ADHERENCE.csv", row.names=1)
a$group <- factor(a$group ,levels =c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA"))
group <- c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA","CrcsC-C57","CrcsC-IgA")
comparisons <- combn(group, 2, simplify = FALSE) 
ggplot(a,aes(group,adhereRatio,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge", width = 0.65)+ 
  theme_classic()+
  scale_fill_manual(values = c( "#9fa0a0","#9fa0a0","#1d519a","#1d519a"))+
  geom_jitter(shape=21,stroke=0.01,alpha=1,size=4.5,
              width = 0.1, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 20)) 

a <- read.csv("./DSSHISTOLOGYSCROE.csv", row.names=1)
a$group <- factor(a$group ,levels =c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA"))
group <- c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA")
ggplot(a,aes(group,score,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge", width = 0.65)+ 
  theme_classic()+ 
  scale_fill_manual(values = c( "#9fa0a0","#1d519a","#8b212e","#ef7e1b"))+
  geom_jitter(shape=21,stroke=0.01,alpha=1,size=4,
              width = 0.2, height = 0,color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 2)) 

a <- read.csv("./Elisa.csv", row.names=1)
group_names <- c( "CON-PBS","CON-WT",  "CON-rcsC", "LY6G-PBS","LY6G-WT", "LY6G-rcsC")
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA))+
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),pch=21,alpha=1,size=2.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  scale_fill_manual(values = c(
  "grey","grey","grey","#5cc3e8","#5cc3e8","#5cc3e8"))+
   theme_classic()+
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 100))

