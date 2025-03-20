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
c <- read.csv("./BU-rcsCCAPSULETHICHNESS.csv", row.names=1)
group_names <- c( "WT","rcsC","C-rcsC")
c$group <- factor( c$group, level=group_names )
ggplot(c,aes(group,capsule,fill=group))+
  geom_beeswarm(cex = 2.5, show.legend = FALSE,pch=21,size=3)+
  stat_summary(fun=mean,fun.ymin=mean,fun.ymax=mean,geom='crossbar',width=0.3,size=0.5,color='black')+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.2)+
  scale_fill_manual(values = c("grey","#5cc3e8","#e95f5c"))+theme_classic()+
  scale_y_continuous(limits = c(0,200))+
  stat_compare_means(comparisons = list(c("WT","rcsC"), c("WT","C-rcsC"), c("C-rcsC","rcsC")) ,
                     method = "wilcox.test",label = "p.signif")

m3h.g2 <-read.delim("./ODVALUE.txt",header = T,row.names= 1)
m3h.g2 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[8]<-"OD600"
group_names <- c( "WT","rcsC","CrcsC")
new.dat$group <- factor(new.dat$group, level=group_names )
ggplot(new.dat,aes(x=Weeks,y=OD600))+
  geom_point(aes(color=group),size=0.91)+
  geom_line(aes(color=group),size=0.91)+
  geom_errorbar(aes(ymin=OD600-std_value,
                    ymax=OD600+std_value, color=group),width=0.42,size=0.45)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,72), breaks = seq(0,72,12))+
  scale_color_manual(values = c("WT"="grey","rcsC"="#5cc3e8","CrcsC"="#e95f5c"))+
  theme_classic()

m3h.g1 <-read.delim("./ODVALUE.txt",header = T,row.names= 1)
m3h.g1 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[9]<-"OD600"
group_names <- c("Pectin", "Glucose","HG","HG-oligos","RG-I","RG-I-oligos")
new.dat$group <- factor(new.dat$group, level=group_names )
group_names <- c( "WT","rcsC","CrcsC")
new.dat$strain <- factor(new.dat$strain, level=group_names )
ggplot(new.dat,aes(x=Weeks,y=OD600))+
  geom_point(aes(color=strain),size=0.51)+
  geom_line(aes(color=strain),size=0.51)+
  geom_errorbar(aes(ymin=OD600-std_value,
                    ymax=OD600+std_value,
                    color=strain),
                width=0.242,
                size=0.5)+
  facet_wrap( ~ group, scales = 'free_x', ncol = 6)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,72),
                     breaks = seq(0,72,12))+
  scale_color_manual(values = c("WT"="grey","rcsC"="#5cc3e8","CrcsC"="#e95f5c"))

c <- read.csv("./IgA-Elisa.csv", row.names=1)
group_names <- c( "Con","WT", "rcsC",  "CrcsC" )
c$Group <- factor(c$Group, level=group_names )

topbar <- function(x){
  return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
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
  scale_y_continuous(limits = c(25, 250), breaks = seq(0, 250, 50))+
  scale_x_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)) +
  scale_color_manual(values = c("#ffdb00","grey", "#5cc3e8", "#e95f5c"))+
  scale_fill_manual(values = c("#ffdb00","grey", "#5cc3e8", "#e95f5c"))+ 
  theme() +
  annotate(geom = 'text',label='***',x=22.5,y=134.2,size=4,angle=90)+
  annotate(geom = 'text',label='****',x=23.8,y=125.12,size=4,angle=90)
