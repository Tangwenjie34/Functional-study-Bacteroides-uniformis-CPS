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

####Figure 1 单形拟杆菌多糖利用能力OD值绘制####
m3h.g2 <-read.delim("E:/AAAAA-博后/单形拟杆菌生长曲线//BU多糖利用能力OD值测定原始数据 .txt",header = T,row.names= 1)
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

ggsave("E:/AAAAA-博后/单形拟杆菌生长曲线/BU多糖利用能力OD值测定原始数据 .pdf",height = 5, width =6.5)

chars_to_extract <- c("0.801308535",
                      "0.8174416854",
                      "0.7730670854",
                      "0.8830872854",
                      "0.59351165",
                      "0.346612148",
                      "0.115813774",
                      "0.121031032",
                      "0.1399394554",
                      "0.1500604232",
                      "0.154451544",
                      "0.1517277294",
                      "0.1393348984")
colnames(new.dat)
new.dat_extracted <- subset(new.dat, OD600 %in% chars_to_extract)
# write.csv(new.dat_extracted, "gene_expression_matrix.csv")
gene_expression_matrix <- read.csv("gene_expression_matrix_多糖.csv", row.names=1)
# gene_expression_matrix_scaled <- t(scale(t(gene_expression_matrix)))


# 创建热图
bk <- c(seq(0,1,by=0.01))

pheatmap(gene_expression_matrix, scale = "none", cluster_rows=T, cluster_cols=T,
         show_rownames=T, show_colnames=TRUE, 
         color=colorRampPalette(c("black", "pink", "firebrick3"))(length(bk)),
         breaks = bk)
# 热图没法用ggsave，直接图片输出


####Figure 2果胶降解能力和引出rcsC基因####
m3h.g1 <-read.delim("E:/AAAAA-博后/单形拟杆菌生长曲线//各类聚糖OD值测定原始数据.txt",header = T,row.names= 1)
m3h.g1 %>% rowwise() %>% 
  mutate(mean_value=mean(c_across(rep1:rep5)),std_value=sd(c_across(rep1:rep5))) %>% 
  filter(Weeks>=1)-> new.dat
colnames(new.dat)[8]<-"OD600"
ggplot(new.dat,aes(x=Weeks,y=OD600))+
  geom_point(aes(color=group),size=1)+
  geom_line(aes(color=group),size=1)+
  geom_errorbar(aes(ymin=OD600-std_value,
                    ymax=OD600+std_value,
                    color=group),
                width=0.42,
                size=1)+
facet_wrap( ~ group, scales = 'free_x', ncol = 5)+
scale_y_continuous(limits = c(0,1))+
scale_x_continuous(limits = c(0,72),
                     breaks = seq(0,72,12))+
scale_color_manual(values = c("Glucose"="#f68a15",
                              "Pectin"="#008c00",
                              "RG-I-oligos"="#094cc3",
                              "RG-II"="#784cc3",
                              "X_Ara_dSB"="#784c23",
                              "HG-oligos"="#22cac3",
                              "Galactan"="#ba3a15",
                              "RG-I"="#AABC15",
                              "HG"="#C19995",
                              "X_Ara_SB"="#D97775"))
ggsave("E:/AAAAA-博后/单形拟杆菌生长曲线/各类聚糖绘制生长曲线.pdf",height = 6, width =12)


chars_to_extract <- c("0.805366",
                      "0.8855918",
                      "0.757035972",
                      "0.437993028",
                      "0.1399394554",
                      "0.1500604232",
                      "0.6626195202",
                      "0.154451544",
                      "0.1517277294",
                      "0.8161498326",
                      "0.1393348984",
                      "0.4375367172",
                      "0.6493279532",
                      "0.5696508",
                      "0.601028356")
colnames(new.dat)
new.dat_extracted <- subset(new.dat, OD600 %in% chars_to_extract)
# write.csv(new.dat_extracted, "gene_expression_matrix.csv")
gene_expression_matrix <- read.csv("gene_expression_matrix.csv", row.names=1)
# gene_expression_matrix_scaled <- t(scale(t(gene_expression_matrix)))

# 创建热图
bk <- c(seq(0,1,by=0.01))

 pheatmap(gene_expression_matrix, scale = "none", cluster_rows=T, cluster_cols=T,
         show_rownames=T, show_colnames=TRUE, 
         color=colorRampPalette(c("black", "pink", "firebrick3"))(length(bk)),
         breaks = bk)
 # 热图没法用ggsave，直接图片输出
 
 
####BU五株不同菌果胶利用能力OD值测定原始数据####
 
 m3h.g3 <-read.delim("E:/AAAAA-博后/单形拟杆菌生长曲线//BU五株不同菌果胶利用能力OD值测定原始数据.txt",header = T,row.names= 1)
 m3h.g3 %>% rowwise() %>% 
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
   facet_wrap( ~ group, scales = 'free_x', ncol = 5)+
   scale_y_continuous(limits = c(0,1))+
   scale_x_continuous(limits = c(0,72),
                      breaks = seq(0,72,12))+
   scale_color_manual(values = c("BU-1"="#f68a15",
                                 "BU-2"="#008c00",
                                 "BU-3"="#094cc3",
                                 "BU-4"="#784cc3",
                                 "BU-5"="#784c23",
                                 "Type II mucin"="#22cac3",
                                 "Inulin"="#ba3a15",
                                 "Laminarin"="#AABC15",
                                 "A_Carbon(NO)"="#C19995",
                                 "Agarose"="#D97775",
                                 "Amylopectin"="#D33775","AStarch"="#344775","Dextran"="#127775"))
 
 ggsave("E:/AAAAA-博后/单形拟杆菌生长曲线/BU五株不同菌果胶利用能力OD值测定.pdf",height = 1.8, width =8)
 

 
 
 ####果胶日粮引起IgA反应性差异化Elisa检测####
 
 a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/果胶日粮引起IgA反应性差异Elisa检测.csv", row.names=1)
 ggplot(a,aes(Cage.group,IgA))+ 
   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
   geom_jitter(aes(fill=Cage.group),pch=21,stroke=0.1,alpha=1,size=2.5,
               position = position_jitter(0.18))+
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
ggsave("果胶日粮引起IgA反应性差异化1.pdf",height =4.5, width =9) 

group_names <- c( "Con","IgA-High",  "IgA-Low")
comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较

ggplot(a,aes(group,IgA))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),pch=21,stroke=0.1,alpha=1,size=3,
              position = position_jitter(0.18))+
    stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  scale_fill_manual(values = c( "grey","#5cc3e8","#e95f5c"))+
            theme(axis.text.x=element_text(angle=30,vjust=1, hjust=1))+
  theme_prism(base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 12,  # 图形的字体大小
              base_line_size = 0.25, # 坐标轴的粗细
              axis_text_angle = 0) 

ggsave("果胶日粮引起IgA反应性差异化分组组合1.pdf",height = 3, width =4.5) 
  
####果胶日粮引起IgA-coating反应性差异化####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/果胶日粮引起IgA反应性差异Elisa检测.csv", row.names=1)
data<-subset(a,a$IgA.Coating!="NA") 
group_names <- c( "Con","IgA-High",  "IgA-Low")
comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较

ggplot(data,aes(group,IgA.Coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),stroke=0,alpha=1,size=3,
              position = position_jitter(0.18),pch=21)+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c( "grey","#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_y_continuous(limits = c(0, 40), breaks = seq(5,30, 5))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.25, # 坐标轴的粗细
              axis_text_angle = 45) 
ggsave("果胶日粮引起IgA-coating反应性差异化1.pdf",height = 4, width =4.5)  

####果胶日粮引起IgA反应性差异化后的cohousing实验ELISA结果####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/果胶日粮引起IgA反应性差异Elisa检测.csv", row.names=1)
data<-subset(a,a$IgA.cohusing!="NA") 
group_names <- c( "IgA-High-Initial","IgA-Low-Initial","IgA-High-post","IgA-Low-post")

data$group.cohousing <- factor(data$group.cohousing,levels = group_names)

ggplot(data,aes(group.cohousing,IgA.cohusing,fill=group.cohousing))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.5, # 坐标轴的粗细
              axis_text_angle = 45)+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c("#dedede", "grey","#dedede","grey"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=2.5,
              position = position_jitter(0.28),color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  geom_signif(comparisons = list(c("IgA-Low-Initial","IgA-Low-post"),
                                 c("IgA-High-post","IgA-Low-post")),# 设置需要比较的组
            map_signif_level = T, #是否使用星号显示
            test = "wilcox.test", ##计算方法
            y_position = c(500,520,550,570),#横线下方的竖线设置
            size=0.8,color="black")+
  scale_y_continuous(limits = c(0, 600), breaks = seq(0, 600, 100)) +
  facet_wrap( ~regroup.cohousing, scales = 'free_x')
ggsave("果胶日粮引起IgA反应性差异化后的cohousing实验ELISA结果.pdf",height = 4, width =5) 


####果胶日粮引起IgA反应性ELISA差异化-各肠段分段比较####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/果胶日粮引起IgA反应性差异各肠段SI-Colon-Elisa检测.csv", row.names=1)
group_names <- c( "IgA-High","IgA-Low")
a$GI.tract <- factor(a$GI.tract,levels =c("Pro-SI","Distal-SI","Cecum","Pro-colon","Distal-colon"))
ggplot(a,aes(group,IgA,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  theme_classic()+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c("#dedede", "grey","#dedede","grey"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=2.5,
              position = position_jitter(0.28),color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  geom_signif(comparisons = list(c("IgA-Low","IgA-High")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = "wilcox.test", ##计算方法
              y_position = c(450),#横线下方的竖线设置
              size=0.25,color="black")+
  scale_y_continuous(limits = c(0, 500), breaks = seq(0, 500, 100)) +
  facet_wrap( ~GI.tract, scales = 'free_x',ncol = 5)
ggsave("果胶日粮引起IgA反应性ELISA差异化各肠段分段比较.pdf",height = 3.4, width =6) 








####所有菌株IgA反应性Elisa检测####
 
 # a <- subset( a, BOR %in% c("cPD","PD","PR","SD"))
 # write.csv(a, "a.csv")
 a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/菌株IgA反应性Elisa检测.csv", row.names=1)
 group_names <- c( "Con","B.uniformis",  "B.ovatus", "B.xylanisolvens","B.faecis", "B.intestinalis", "B.caccae",
                  "B.vulgatus", "P.copri", "P.intermedia", "P.mirabilis", "F.varium")
 comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
 # comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较
a$group <- factor( a$group, level=group_names )

 ggplot(a,aes(group,IgA))+
   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
   geom_jitter(aes(fill=group),pch=21,alpha=1,size=2.5,
               position = position_jitter(0.18))+
     stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
   stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
   scale_fill_manual(values = c(
     "grey","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
           "#e95f5c","#e95f5c","#ffdb00","#79ceb8"))+
  theme_classic()+
   scale_y_continuous(limits = c(0, 800), breaks = seq(0, 800, 100))
 
 ggsave("菌株IgA反应性Elisa检测.pdf",height = 5.5, width =7) 
 # method:比较的方法，
 # 默认为”wilcox.test”, 其他可选方法为：”t.test”、”anova”、”kruskal.test”
 
 
 ####菌株诱导粪便IgA-coating反应性流式检测####
 a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/菌株IgA反应性Iga-coating值.csv", row.names=1)
 group_names <- c( "Con","B.uniformis",  "B.ovatus", "B.xylanisolvens","B.faecis", "B.intestinalis", "B.caccae",
                   "B.vulgatus", "P.copri", "P.intermedia", "P.mirabilis", "F.varium")
 a$group <- factor( a$group, level=group_names )
 comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
 
 ggplot(a,aes(group,IgA_coating))+ 
   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
   geom_jitter(aes(fill=group),shape=21,stroke=0.1,alpha=1,size=2.5,
               position = position_jitter(0.18))+
   stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
   stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
   stat_compare_means(method="anova")+
   scale_fill_manual(values = c(
     "grey","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8","#5cc3e8",
     "#e95f5c","#e95f5c","#ffdb00","#79ceb8"))+
   theme_classic()+
   scale_y_continuous(limits = c(0, 40), breaks = seq(0,30, 5))
 
ggsave("菌株诱导粪便IgA-coating反应性流式检测.pdf",height = 5.5, width =8.5) 
 
 ####粪便中Iga_coating与Iga含量的相关性#### 
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/果胶日粮引起IgA反应性差异Elisa检测.csv", row.names=1)
 ggplot(data = a ,mapping = aes( x = IgA, y= IgA.Coating)) + 
   geom_point(size= 2) +
   geom_jitter(position = position_jitter(1), size=1, alpha=1)+
   stat_smooth(method = "lm")+ 
   ggpubr::stat_cor(aes(), label.x =200, label.y = 28)+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30, 4))+
  scale_x_continuous(limits = c(200, 480), breaks = seq(200,480, 100))
 
ggsave("粪便中Iga_coating与Iga含量的相关性.pdf",height = 5, width =5)  


####联合果胶BU中对Iga_coating含量的流式测定#### 
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/植物多糖-BU联合对IgA反应性的影响-四周Elisa检测.csv", row.names=1)
c <- subset(c, c$Group1 %in% c("Con_1","PEC+BU_1" ,  "NOPEC+BU_1"))
group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU")
c$Group <- factor(c$Group, level=group_names )
comparisons <- combn(group_names, 2, simplify = FALSE)

ggplot(c,aes(Group,IgA_coating))+ 
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=Group),shape=21,stroke=0.1,alpha=1,size=2.5,
              position = position_jitter(0.18))+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
  "grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
            # theme(axis.text.x=element_text(angle=0,vjust=1, hjust=1))+
  scale_y_continuous(limits = c(0, 40), breaks = seq(0,30, 5))+
  theme_prism(palette = "candy_bright",
              base_fontface = "plain", # 字体样式，可选 bold, plain, italic
              base_family = "serif", # 字体格式，可选 serif, sans, mono, Arial等
              base_size = 16,  # 图形的字体大小
              base_line_size = 0.25, # 坐标轴的粗细
              axis_text_angle = 45)

ggsave("联合果胶BU中对Iga_coating含量的流式测.pdf",height = 4, width =4)  


####联合果胶BU中对B220-Iga细胞亚群含量的流式测定#### 
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/植物多糖-BU联合对IgA反应性的影响-四周Elisa检测.csv", row.names=1)
c <- subset(c, c$Group1 %in% c("Con_1","PEC+BU_1" ,  "NOPEC+BU_1"))

group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU")
c$Group <- factor(c$Group, level=group_names )
comparisons <- combn(group_names, 2, simplify = FALSE)
ggplot(c,aes(Group,B220_IgA_Subset))+ geom_boxplot(width = 0.46)+
  geom_jitter(aes(fill=Group),shape=21,stroke=0.1,alpha=1,size=4,
              position = position_jitter(0.18))+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
    "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
               # theme(axis.text.x=element_text(angle=0,vjust=1, hjust=1))+
  scale_y_continuous(limits = c(3, 30), breaks = seq(3,30, 6))

ggsave("联合果胶BU中对B220-Iga细胞亚群含量的流式测定.pdf",height = 3.5, width =4.5)  


####联合果胶BU中对B220阴性-Iga细胞亚群含量的流式测定#### 
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/植物多糖-BU联合对IgA反应性的影响-四周Elisa检测.csv", row.names=1)
c <- subset(c, c$Group1 %in% c("Con_1","PEC+BU_1" ,  "NOPEC+BU_1"))

group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU")
c$Group <- factor(c$Group, level=group_names )
comparisons <- combn(group_names, 2, simplify = FALSE)
ggplot(c,aes(Group,NOB220_IgA_Subset))+ geom_boxplot(width = 0.46)+
  geom_jitter(aes(fill=Group),shape=21,stroke=0.1,alpha=1,size=4,
              position = position_jitter(0.18))+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
    "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
               # theme(axis.text.x=element_text(angle=0,vjust=1, hjust=1))+
  scale_y_continuous(limits = c(0, 12), breaks = seq(0,12, 4))

ggsave("联合果胶BU中对B220阴性-Iga细胞亚群含量的流式测定.pdf",height = 3.5, width =4.5)  

 
 ####单独饲喂植物多糖IgA反应性处理四周Elisa检测####
b <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/植物多糖IgA反应性处理四周Elisa检测.csv", row.names=1)
 group_names <- c( "CON",  "PEC", "Starch","XYL","Inulin")
 comparisons <- lapply(group_names[-12], function(group) c(group, "CON"))
 # comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较
   b$Group <- factor(b$Group, level=group_names )
 
   topbar <- function(x){
     return(mean(x)+sd(x)/sqrt(length(x))) #误差采用了mean+-sem
   }
   bottombar <- function(x){
     return(mean(x)-sd(x)/sqrt(length(x)))
   }
   
ggplot(b,aes(week,IgA,color=Group))+
     stat_summary(geom = 'line',fun='mean',cex=1)+
     stat_summary(geom = 'errorbar',
                  fun.min = bottombar,fun.max = topbar,
                  width=0.2,cex=0.8,aes(color=Group))+
     stat_summary(geom = 'point',fun='mean',aes(fill=Group),
                  size=2.5,pch=21,color='black')+
     theme_classic(base_size = 16)+
     theme(legend.position = 'none')+
  facet_wrap( ~Group, scales = 'free_x', ncol = 5)+
  scale_y_continuous(limits = c(100, 400), breaks = seq(0, 400, 100))+
  scale_color_manual(values = c("grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_fill_manual(values = c("grey", "#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))
ggsave("单独饲喂多糖IgA反应性Elisa检测.pdf",height = 3.5, width =8)   
  
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
ggsave("单独饲喂多糖IgA反应性Elisa检测w.pdf", height = 2.7, width =8)  
 
   
####植物多糖+BU之后IgA反应性处理四周Elisa检测####
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/植物多糖-BU联合对IgA反应性的影响-四周Elisa检测.csv", row.names=1)
group_names <- c( "Con","PEC+BU" ,  "NOPEC+BU"  ,     "XYL+BU"   ,  "NOXYL+BU" , 
                 "Inulin+BU" ,   "NO_Inulin+BU" )
# comparisons <- lapply(group_names[-12], function(group) c(group, "CON"))
# comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较
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
  # facet_wrap( ~Group, scales = 'free_x', ncol = 5)+
  scale_y_continuous(limits = c(230, 400), breaks = seq(0, 400, 50))+
  scale_x_continuous(limits = c(1, 7.5), breaks = seq(1, 7, 1)) +
  scale_color_manual(values = c("grey", "#5cc3e8","#3b8000","#0151b8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_fill_manual(values = c("grey", "#5cc3e8", "#3b8000","#0151b8","#e95f5c","#ffdb00","#79ceb8"))+  # 添加图例标签
  theme() +
  annotate(geom = 'text',label='***',x=6.5,y=342,size=4,angle=90)+
  annotate(geom = 'text',label='****',x=6.8,y=312,size=4,angle=90)

# 显著性在graphprim里面统计的
ggsave("植物多糖_BU之后IgA反应性处理四周Elisa检测.pdf",height = 3.25, width =4.5)   

   
####Figure 2_BU1-5诱导IgA-coating反应性的同一菌株之间的功能差异####   
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/Figure 2_BU1-5诱导IgA-coating反应性的同一菌株之间的功能差异.csv", row.names=1)
group_names <- c( "PECCON","PECBU1",  "PECBU2", "PECBU3","PECBU4", "PECBU5",
                  "NOPECCON", "NOPECBU1", "NOPECBU2", "NOPECBU3", "NOPECBU4","NOPECBU5")
a$group <- factor( a$group, level=group_names )

ggplot(a,aes(group,IgA_coating))+
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),stroke=0.1,alpha=1,size=2.5,
              position = position_jitter(0.18),pch=21)+
  stat_compare_means(comparisons = list(c("PECCON","PECBU1"), c("PECBU1","PECBU4")) ,
                     method = "wilcox.test",label = "p.signif")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
  stat_compare_means(method="anova")+
  scale_fill_manual(values = c(
    "grey","#5cc3e8","#5cc3e8","#5cc3e8","#e95f5c","#e95f5c","grey","#5cc3e8","#5cc3e8","#5cc3e8",
          "#e95f5c","#e95f5c"))+
            theme(axis.text.x=element_text(angle=30,vjust=1, hjust=1))+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30, 5))+
   facet_wrap( ~PEC, scales = 'free_x', ncol = 5)

ggsave("Figure 2_BU1-5诱导IgA-coating反应性的同一菌株之间的功能差异.pdf",height = 3, width =6.5) 



####Figure 3_BU1-5荚膜厚度统计图####   
b <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/Figure 3_BU1-5荚膜厚度统计.csv", row.names=1)
group_names <- c( "CON","BU1",  "BU2", "BU3","BU4", "BU5")
b$group <- factor( b$group, level=group_names )

# ggplot(b,aes(group,capsule,fill=group))+
#   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+ #绘制柱状图 #误差棒
#   geom_jitter(stroke=0.1,alpha=1,size=2.5, position = position_jitter(0.18),pch=21)+
#   stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
#   stat_compare_means(comparisons = list(c("CON","BU1"), c("BU1","BU4"), c("BU5","BU4")) ,
#                      method = "wilcox.test",label = "p.signif")+
#   scale_fill_manual(values = c("grey","#5cc3e8","#5cc3e8","#5cc3e8","#e95f5c","#e95f5c"))+theme_classic()
# ggsave("Figure 3_BU1-5荚膜厚度统计.pdf",height = 3.2, width =5)

ggplot(b,aes(group,capsule,fill=group))+
  geom_beeswarm(cex = 1.82, show.legend = FALSE,pch=21,size=1.85)+
  stat_summary(fun=mean,fun.ymin=mean,fun.ymax=mean,geom='crossbar',width=0.3,size=0.3,color='black')+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.25,size=0.25)+
  scale_fill_manual(values = c("grey","#5cc3e8","#5cc3e8","#5cc3e8","#e95f5c","#e95f5c"))+theme_classic()+
  scale_y_continuous(limits = c(0,250))+
  stat_compare_means(comparisons = list(c("CON","BU1"), c("BU1","BU4"), c("BU5","BU4")) ,
                     method = "wilcox.test",label = "p.signif")
ggsave("Figure 3_BU1-5荚膜厚度统计1.pdf",height = 3.2, width =5)

####Figure 3_BU-rcsC基因敲除对荚膜厚度影响统计图####   
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/Figure 3_BU-rcsC突变株荚膜厚度统计.csv", row.names=1)
group_names <- c( "WT","rcsC","C-rcsC")
c$group <- factor( c$group, level=group_names )

# ggplot(c,aes(group,capsule,fill=group))+
#   geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+ #绘制柱状图 #误差棒
#   geom_jitter(stroke=0.1,alpha=1,size=3.5, position = position_jitter(0.18),pch=21)+
#   stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.35)+
#   stat_compare_means(comparisons = list(c("WT","rcsC"), c("WT","C-rcsC"), c("C-rcsC","rcsC")) ,
#                      method = "wilcox.test",label = "p.signif")+
#   scale_fill_manual(values = c("grey","#5cc3e8","#e95f5c"))+theme_classic()
# 
# ggsave("Figure 3_BU-rcsC基因敲除对荚膜厚度影响统计图.pdf",height = 3.5, width =5)

ggplot(c,aes(group,capsule,fill=group))+
  geom_beeswarm(cex = 2.5, show.legend = FALSE,pch=21,size=3)+
  stat_summary(fun=mean,fun.ymin=mean,fun.ymax=mean,geom='crossbar',width=0.3,size=0.5,color='black')+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.2)+
  scale_fill_manual(values = c("grey","#5cc3e8","#e95f5c"))+theme_classic()+
  scale_y_continuous(limits = c(0,200))+
  stat_compare_means(comparisons = list(c("WT","rcsC"), c("WT","C-rcsC"), c("C-rcsC","rcsC")) ,
                     method = "wilcox.test",label = "p.signif")

ggsave("Figure 3_BU-rcsC基因敲除对荚膜厚度影响统计图1.pdf",height = 3.5, width =5)


####Figure 3_BU-rcsC基因敲除果胶利用能力OD值绘制####
m3h.g2 <-read.delim("E:/AAAAA-博后/单形拟杆菌生长曲线//Figure 3_BU-rcsC基因敲除果胶利用能力OD值绘制.txt",header = T,row.names= 1)
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

ggsave("E:/AAAAA-博后/单形拟杆菌生长曲线//Figure 3_BU-rcsC基因敲除果胶利用能力OD值绘制.pdf",
       height = 2.5, width =4.5)

####Figure S8C_BU-rcsC基因敲除对寡糖和单糖利用能力OD值绘制####
m3h.g1 <-read.delim("E:/AAAAA-博后/单形拟杆菌生长曲线/Figure S8C_BU-rcsC基因敲除对寡糖和单糖利用能力OD值绘制.txt",header = T,row.names= 1)
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
ggsave("E:/AAAAA-博后/单形拟杆菌生长曲线/Figure S8C_BU-rcsC基因敲除对寡糖和单糖利用能力OD值绘制.pdf",
       height = 2, width =11)



####散点加误差棒图代码####
ggplot(c,aes(Group,B220_IgA_Subset)) + 
  geom_jitter(aes(fill=Group),shape=21,stroke=0.1,alpha=1,size=4,
              position = position_jitter(0.18))+
  stat_summary(fun = mean, fun.ymin = mean, fun.ymax = mean,
               geom = 'crossbar', width = 0.3, size = 0.3, color = 'black')+  #计算各组均值并添加在图中
  stat_summary(fun.data = function(x) mean_se(x),
               geom = 'errorbar', width = 0.25, size = 0.2, color = 'black') +
  scale_fill_manual(values = c("#5cc3e8", "#e95f5c","#ffdb00","#79ceb8"))+
  scale_y_continuous(limits = c(0, 30), breaks = seq(0,30, 6))+  theme_classic()



####无菌小鼠BU之后IgA反应性处理三周Elisa检测####
c <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/无菌小鼠BU之后IgA反应性处理三周Elisa检测.csv", row.names=1)
group_names <- c( "Con","WT", "rcsC",  "CrcsC" )
# comparisons <- lapply(group_names[-12], function(group) c(group, "CON"))
# comparisons <- combn(group_names, 2, simplify = FALSE) #任意两组比较
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
  # facet_wrap( ~Group, scales = 'free_x', ncol = 5)+
  scale_y_continuous(limits = c(25, 250), breaks = seq(0, 250, 50))+
  scale_x_continuous(limits = c(5, 25), breaks = seq(5, 25, 5)) +
  scale_color_manual(values = c("#ffdb00","grey", "#5cc3e8", "#e95f5c"))+
  scale_fill_manual(values = c("#ffdb00","grey", "#5cc3e8", "#e95f5c"))+  # 添加图例标签
  theme() +
  annotate(geom = 'text',label='***',x=22.5,y=134.2,size=4,angle=90)+
  annotate(geom = 'text',label='****',x=23.8,y=125.12,size=4,angle=90)

# 显著性在graphprim里面统计的
ggsave("无菌小鼠BU之后IgA反应性处理三周Elisa检测.pdf",height = 3.25, width =4.5)   

####无菌小鼠BU之后IgA反应性处理三周IgA和rcsC敲除珠CFU相关性#### 
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/无菌小鼠BU之后IgA反应性处理三周Elisa检测.csv", row.names=1)
ggplot(data = a ,mapping = aes( x = IgA, y= logCFUg)) + 
  geom_point(size= 2) +
  geom_jitter(position = position_jitter(1), size=1, alpha=1)+
  stat_smooth(method = "lm")+ 
  ggpubr::stat_cor(aes(), label.x =20, label.y = 9)+
  scale_y_continuous(limits = c(0, 10), breaks = seq(0,10, 2))+
  scale_x_continuous(limits = c(20, 220), breaks = seq(20,220, 50))+theme_classic2()

ggsave("无菌小鼠BU之后IgA反应性处理三周IgA和rcsC敲除珠CFU相关性.pdf",height = 5, width =5)  


####无菌小鼠BU之后肠腔-上皮-黏膜组织CFU比较####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/无菌小鼠BU之后肠腔-上皮-黏膜组织CFU比较.csv", row.names=1)
a$position <- factor(a$position ,levels =c("Lumen","Epithelium","Tissue"))
a$group <- factor(a$group ,levels =c("WT","rcsC"))
p1 <-ggplot(a,aes(group,logCFUg,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  theme_classic()+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c( "#5f5f60","#1754a5"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=3.5,
              position = position_jitter(0.28),color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  geom_signif(comparisons = list(c("WT","rcsC")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = "wilcox.test", ##计算方法
              y_position = c(9),#横线下方的竖线设置
              size=0.25,color="black")+
  scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2)) +
  facet_wrap( ~position, scales = 'free_x',ncol = 5)

p2 <-ggplot(a,aes(group,Ration,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge")+ #绘制柱状图
  theme_classic()+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c( "#5f5f60","#1754a5"))+
  geom_jitter(shape=21,stroke=0.1,alpha=1,size=5,
              position = position_jitter(0.28),color="black")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  geom_signif(comparisons = list(c("WT","rcsC")),# 设置需要比较的组
              map_signif_level = T, #是否使用星号显示
              test = "wilcox.test", ##计算方法
              y_position = c(1.2),#横线下方的竖线设置
              size=0.25,color="black")+
  scale_y_continuous(limits = c(0, 1.5), breaks = seq(0, 1.2, 0.2))
combined_plot <- p1 + p2
print(combined_plot)
ggsave("无菌小鼠BU之后肠腔-上皮-黏膜组织CFU比较.pdf",height = 3.4, width =8.6) 

####C57BL6和IgA敲除鼠BU之后粘附细胞比例计算####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/C57BL6和IgA敲除鼠BU之后粘附细胞比例计算.csv", row.names=1)
a$group <- factor(a$group ,levels =c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA","CrcsC-C57","CrcsC-IgA"))

group <- c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA","CrcsC-C57","CrcsC-IgA")
# comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
comparisons <- combn(group, 2, simplify = FALSE) #任意两组比较

ggplot(a,aes(group,adhereRatio,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge", width = 0.65)+ #绘制柱状图
  theme_classic()+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c( "#9fa0a0","#9fa0a0","#1d519a","#1d519a","#8b212e","#8b212e"))+
  geom_jitter(shape=21,stroke=0.01,alpha=1,size=4.5,
              position = position_jitter(0.15),color="black")+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  scale_y_continuous(limits = c(0, 160), breaks = seq(0, 160, 20)) 

ggsave("C57BL6和IgA敲除鼠BU之后粘附细胞比例计算.pdf",height = 5.5, width =7) 


####IgA敲除小鼠DSS之后炎症评分统计####
a <- read.csv("E:/AAAAA-博后/单形拟杆菌生长曲线/IgA敲除小鼠DSS之后炎症评分统计.csv", row.names=1)
a$group <- factor(a$group ,levels =c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA"))

group <- c("WT-C57","WT-IgA","rcsC-C57","rcsC-IgA")
# comparisons <- lapply(group_names[-13], function(group) c(group, "B.uniformis"))
comparisons <- combn(group, 2, simplify = FALSE) #任意两组比较

ggplot(a,aes(group,score,fill=group))+
  geom_bar(stat="summary",fun=mean,position="dodge", width = 0.65)+ #绘制柱状图
  theme_classic()+ # 可选值有 0，45，90，270
  scale_fill_manual(values = c( "#9fa0a0","#1d519a","#8b212e","#ef7e1b"))+
  geom_jitter(shape=21,stroke=0.01,alpha=1,size=4,
              position = position_jitter(0.15),color="black")+
  stat_compare_means(comparisons = comparisons, method = "wilcox.test",label = "p.signif")+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.3)+#误差棒
  scale_y_continuous(limits = c(0, 15), breaks = seq(0, 15, 2)) 

ggsave("IgA敲除小鼠DSS之后炎症评分统计.pdf",height = 5.5, width =7) 



