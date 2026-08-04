####Data Processing and Visualization 
##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
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
library(ggrepel)
library(RColorBrewer)
library(viridis)

a <- read.csv("./data1.csv", row.names=1)
df_subset <- a[,c("X47", "X54","X61","X98","X119", "taxonomy") ]
df_subset$taxonomy[is.na(df_subset$taxonomy) | df_subset$taxonomy == ""] <- "other"
df_summary <- df_subset %>%
  group_by(taxonomy) %>%
  summarise(
    X47 = sum(X47, na.rm = TRUE),
    X54 = sum(X54, na.rm = TRUE),
    X61 = sum(X61, na.rm = TRUE),
    X98 = sum(X98, na.rm = TRUE),
    X119 = sum(X119, na.rm = TRUE)
  ) %>%
  ungroup()
df_percent <- df_summary %>%
  mutate(across(c(X47, X54, X61, X98, X119), 
                ~ . / sum(., na.rm = TRUE) * 100,
                .names = "{.col}_percent"))
write.csv(df_percent, "df_percent.csv")
b <- read.csv("./df_percent.csv", row.names=1)
 pie(b$FMTP_A, labels = b$taxonomy, main = "FMTP_A") 
 pie(b$FMT_A,  labels = b$taxonomy,
          main = "FMT_A") 

taxa_data <- read.csv("./data.csv", header = TRUE, row.names = 1, check.names = FALSE)
taxa_datat <- data.frame(t(taxa_data))
taxa_datat$Sample <- rownames(taxa_datat)
taxa_long <- taxa_datat %>%pivot_longer( cols = -Sample, names_to = "Genus",values_to = "Abundance")
N <- 5
top_genera <- taxa_long %>%
  group_by(Genus) %>%
  summarise(MeanAbundance = mean(Abundance, na.rm = TRUE)) %>%
  arrange(desc(MeanAbundance)) %>%
  slice_head(n = N) %>%
  pull(Genus)
taxa_plot <- taxa_long %>%
  mutate(
    Genus_plot = ifelse(Genus %in% top_genera, Genus, "Other")
  ) %>%
  group_by(Sample, Genus_plot) %>%
  summarise(
    Abundance = sum(Abundance, na.rm = TRUE),
    .groups = "drop" )
ggplot(taxa_plot, aes(x = Sample, y = Abundance, fill = Genus_plot)) +
  geom_bar(stat = "identity", position = "fill", width = 0.85) +
  scale_fill_manual(values = c("g__Bacteroides" = "#b2182b",
                               "g__Allobaculum" = "#d9d9d9",
                               "g__Ligilactobacillus" = "#d9d9d9",
                               "g__norank_f__Muribaculaceae" = "#d9d9d9",
                               "g__Candidatus_Saccharimonas" = "#d9d9d9",
                               "Other"="#d9d9d9"), name = "Genus") +
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust = 1,size = 10))

a <- read.csv("./data1.csv", row.names=1)
 ggplot(data = a ,mapping = aes( x = IgA, y= IgA.Coating)) + 
   geom_point(size= 2) +
   stat_smooth(method = "lm")+ 
   ggpubr::stat_cor(aes(), label.x =200, label.y = 28)+
  scale_x_continuous(limits = c(200, 480), breaks = seq(200,480, 100))

a <- read.csv("./data2.csv", row.names=1)
 ggplot(a,aes(Cage.group,IgA))+ 
   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
   geom_jitter(aes(fill=Cage.group),pch=21,stroke=0.1,alpha=1,size=2.5,
               width = 0.1, height = 0)+
   stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
   scale_fill_manual(values = c("grey","grey","grey","grey","grey",
                                "#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
                                "#e95f5c", "#e95f5c", "#e95f5c", "#e95f5c", "#e95f5c"))+
   scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100))+
   theme_prism(palette = "candy_bright",
               base_fontface = "plain", # 字体样式，可选 bold, plain, italic
               base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
               base_size = 16,  # 图形的字体大小
               base_line_size = 0.25, # 坐标轴的粗细
               axis_text_angle = 45)+ # 可选值有 0，45，90，270
  facet_wrap( ~group, scales = 'free_x', ncol = 3)
group_names <- c( "Con","IgA-High",  "IgA-Low")
comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较

ggplot(a,aes(group,IgA))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),pch=21,stroke=0.1,alpha=1,size=3,
              width = 0.1, height = 0)+
    stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  scale_fill_manual(values = c( "grey","#5cc3e8","#e95f5c"))+
            theme(axis.text.x=element_text(angle=30,vjust=1, hjust=1))+
  theme_prism(base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 12,  # 图形的字体大小
              base_line_size = 0.25, # 坐标轴的粗细
              axis_text_angle = 0) 


gene_expression_matrix <- read.csv("./data3.csv", row.names=1)
bk <- c(seq(0,1,by=0.01))
pheatmap(gene_expression_matrix, scale = "none", cluster_rows=T, cluster_cols=T,
         show_rownames=T, show_colnames=TRUE, 
         color=colorRampPalette(c("black", "pink", "firebrick3"))(length(bk)),
         breaks = bk)


m3h.g2 <-read.delim("./Bu.txt",header = T,row.names= 1)
m3h.g2 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[8]<-"OD600"
ggplot(new.dat,aes(x=Weeks,y=OD600))+
  geom_point(aes(color=group),size=0.71)+
  geom_line(aes(color=group),size=0.71)+
  geom_errorbar(aes(ymin=OD600-std_value,
                    ymax=OD600+std_value,
                    color=group),
                width=0.32,
                size=0.35)+
  scale_y_continuous(limits = c(0,1))+
  scale_x_continuous(limits = c(0,72),
                     breaks = seq(0,72,12))+
  scale_color_manual(values = c("Glucose"="#f68a15",
                                "Pectin"="#008c00",
                                "Xyloglucan"="#094cc3",
                                "Pectin_galactan"="#784cc3",
                                "Arabinogalactan"="#784c23",
                                "Type II mucin"="#22cac3",
                                "Inulin"="#ba3a15",
                                "Laminarin"="#AABC15",
                                "A_Carbon(NO)"="#C19995",
                                "Agarose"="#D97775",
                                "Amylopectin"="#D33775","AStarch"="#344775","Dextran"="#127775"))

b <- read.csv("./Elisa检测.csv", row.names=1)
 group_names <- c( "CON",  "PEC", "Starch","XYL","Inulin")
 comparisons <- lapply(group_names[-12], function(group) c(group, "CON"))
 # comparisons <- combn(group_names, 2, simplify = FALSE) 
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
