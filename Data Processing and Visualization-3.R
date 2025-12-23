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



m3h.g2 <-read.delim("./ODVALUE.txt",header = T,row.names= 1)
m3h.g2 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[8]<-"OD600"
group_names <- c( "WT","rcsC")
new.dat$group <- factor(new.dat$group, level=group_names )
ggplot(new.dat,aes(x=Weeks,y=OD600))+
  geom_point(aes(color=group),size=0.91)+
  geom_line(aes(color=group),size=0.91)+
  geom_errorbar(aes(ymin=OD600-std_value,
                    ymax=OD600+std_value, color=group),width=0.42,size=0.45)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,72), breaks = seq(0,72,12))+
  scale_color_manual(values = c("WT"="grey","rcsC"="#5cc3e8")+
  theme_classic()

m3h.g1 <-read.delim("./ODVALUE.txt",header = T,row.names= 1)
m3h.g1 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[9]<-"OD600"
group_names <- c("Pectin", "Glucose","HG","HG-oligos","RG-I","RG-I-oligos")
new.dat$group <- factor(new.dat$group, level=group_names )
group_names <- c( "WT","rcsC")
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
  scale_color_manual(values = c("WT"="grey","rcsC"="#5cc3e8"))

c <- read.csv("./IgA-Elisa.csv", row.names=1)
group_names <- c( "Con","WT", "rcsC" )
c$Group <- factor(c$Group, level=group_names )

topbar <- function(x){
  return(mean(x)+sd(x)/sqrt(length(x))) #mean+-sem
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
  scale_color_manual(values = c("#ffdb00","grey", "#5cc3e8"))+
  scale_fill_manual(values = c("#ffdb00","grey", "#5cc3e8"))+ 
  theme() 
                     
 a <- read.csv("./Iga-coating.csv", row.names=1)
group_names <- c( "WT","rcsC")
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),shape=21,stroke=0.1,alpha=1,size=2,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  scale_fill_manual(values = c(
    "grey","#5cc3e8"))+
  theme_classic()+
  scale_y_continuous(limits = c(0, 22), breaks = seq(0,20, 5))

a <- read.csv("./Iga-coating_1.csv", row.names=1)
group_names <- c("HP-WT","HP-rcsC" , "LP-WT", "LP-rcsC")
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),shape=21,stroke=0.1,alpha=1,size=1.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  scale_fill_manual(values = c("grey","#5cc3e8","grey","#5cc3e8"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,40), breaks = seq(0,25, 5))

a <- read.csv("IgA-LPMC.csv", row.names=1)
group_names <- c("LPMC-SPF-WT","LPMC-SPF-rcsC" , "LPMC-GF-WT", "LPMC-GF-rcsC")
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),shape=21,stroke=0.1,alpha=1,size=1.5,
              width = 0.1, height = 0)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  scale_fill_manual(values = c("grey","#5cc3e8","grey","#5cc3e8"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,15), breaks = seq(0,15, 3))                   
