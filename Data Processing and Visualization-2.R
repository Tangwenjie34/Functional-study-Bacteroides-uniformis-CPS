####Data Processing and Visualization 
Sys.setenv(LANGUAGE = "en")
options(stringsAsFactors = FALSE)
rm(list=ls())
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(patchwork)
library(ggsankey)
library(gggenes)

setwd("/report")
diff_expr_genes <- read.csv("diff_expr_genes.csv")
all_expr_genes <- read.delim("total_gene_information.xls")
log2fc_threshold <- 1
pvalue_threshold <- 0.05
genes_to_label <- c("rcsC")
a <- ggplot(diff_expr_genes, aes(x=log2FoldChange, y=-log10(p.value), color= group)) +
  geom_point(aes(color=significant), alpha=1, size=2.5) +
  scale_color_manual(values=c("Not Significant"="blue", "Significant"="red")) +
  geom_vline(xintercept=c(-log2fc_threshold, log2fc_threshold), col="black", linetype="dashed") +
  geom_hline(yintercept=-log10(pvalue_threshold), col="black", linetype="dashed") +
  geom_text_repel(data=filter(diff_expr_genes, gene_name %in% genes_to_label), 
                  aes(label=gene_name), size=5, box.padding=0.3) +
  labs(x="Log2 Fold Change", y="-Log10 p-value") +
  ylim(0, 90)+
  xlim(-5, 12)+
  theme_pubr()
gene_expression_matrix <- read.csv("gene_expression_matrix.csv", row.names=1)
gene_expression_matrix_scaled <- t(scale(t(gene_expression_matrix)))
pheatmap(gene_expression_matrix, scale = "row", cluster_rows=TRUE, cluster_cols=TRUE,
         show_rownames=T, show_colnames=TRUE, 
         color=colorRampPalette(c("navy", "white", "firebrick3"))(50))

b <- pheatmap(gene_expression_matrix_scaled,  cluster_rows=TRUE, cluster_cols=TRUE,
              show_rownames=T, show_colnames=TRUE, 
              color=colorRampPalette(c("navy", "white", "red"))(50))

cowplot::plot_grid( a, b$gtable,ncol= 2, labels=LETTERS[1:2])

sankeyGoPlot <- function(goData = NULL,
                         topGenes = 5,
                         sankeyExpand = c(0.5,1),
                         nodeSize = 2.5,
                         nodeColor = NULL,
                         goCol = NULL,
                         xShift = 0.05,
                         downShift = 4.5,
                         upShift = 0.25){
  ego_df <- goData %>%
    dplyr::group_by(Description) %>%
    dplyr::mutate(gene_ratio = eval(parse(text = GeneRatio))) %>%
    dplyr::arrange(pvalue)
  ego_df$Description <- factor(ego_df$Description,levels = rev(ego_df$Description))
  ego_df$term_gene <- sapply(strsplit(ego_df$BgRatio,split = "\\/"),"[",1) %>%
    as.numeric()
  ego_df$term_ratio <- ego_df$term_gene
  sankey_df <- ego_df %>%
    dplyr::select(Description,geneID) %>%
    tidyr::separate_longer_delim(geneID,delim = "/") %>%
    dplyr::group_by(Description) %>%
    dplyr::slice_head(n = topGenes) %>%
    as.data.frame() %>%
    dplyr::mutate(Description = as.character(Description))
  
  sankey_df_long <- sankey_df %>% ggsankey::make_long(geneID,Description)
  sankey_df_long$node <- factor(sankey_df_long$node,
                                levels = c(unique(sankey_df$Description),
                                           unique(sankey_df$geneID)))
  if(is.null(nodeColor)){
    mycol <- cols4all::c4a('rainbow_wh_rd',length(unique(sankey_df_long$node)))
  }else{
    mycol <- grDevices::colorRampPalette(colors = nodeColor)(length(unique(sankey_df_long$node)))
  }
  ps <-
    ggplot(data = sankey_df_long,
           mapping = aes(x = x,
                         next_x = next_x,
                         node = node,
                         next_node = next_node,
                         fill = factor(node),
                         label = node)) +
    geom_sankey(flow.alpha = 0.5,
                flow.fill = 'grey',
                flow.color = 'grey',
                width = 0.1,
                node.fill = mycol) +
    theme_void() +
    scale_x_discrete(expand = expansion(mult = sankeyExpand)) +
    theme(legend.position = 'none')
  ps_data <- ggplot_build(ps)
  ps_data_info <- ps_data$data[[2]]
  ps2 <-
    ps +
    geom_text(data = ps_data_info,
              mapping = aes(next_x = 0,next_node = 0,
                            x = xmin,y = (ymin + ymax)/2,
                            label = label,hjust = 1),
              fontface = "bold",size = nodeSize)

  if(is.null(goCol)){
    pcol <- scale_fill_viridis_c(option = "plasma",direction = -1)
  }else{
    pcol <- scale_fill_gradient(low = goCol[1],high = goCol[2])
  }
  pp <-
    ggplot(ego_df) +
    geom_point(aes(x = -log10(pvalue),y = Description,
                   size = term_ratio,fill = gene_ratio),
               color = "black",shape = 21) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text = element_text(colour = "black"),
          plot.background = element_blank()) +
    ylab("") + pcol
  xmin <- max(ps_data_info$xmax)
  y_range <- subset(ps_data_info,x == "2")
  ymin <- min(c(y_range$ymin,y_range$ymax))
  ymax <- max(c(y_range$ymin,y_range$ymax))
  cb <-
    ps2 +
    annotation_custom(grob = ggplotGrob(pp),
                      xmin = xmin - xShift,xmax = 2 + sankeyExpand[2],
                      ymin = ymin - downShift,ymax = ymax + upShift)
  return(cb)
}
ego_df1 <-read.delim("./GO-enrichment.txt",row.names = 1)
sankeyGoPlot(goData = ego_df1)
c <- sankeyGoPlot(goData = ego_df1,
                  downShift = 6.5,upShift = 1,
                  sankeyExpand = c(0.5,1))
ego_df2 <-read.delim("./KEGG-enrichment.txt",row.names = 1)
d <- sankeyGoPlot(goData = ego_df2,
                  downShift = 6.5,upShift = 1,
                  sankeyExpand = c(0.5,1))
cowplot::plot_grid( c, d ,ncol= 2, labels=LETTERS[3:4])

BU_genes <-read.delim("./total_gene_information.xls")

A_CAZy_ara <-read.delim("./A_CAZy_ara.txt",row.names = 1)

colnames(A_CAZy_ara)
ggplot(A_CAZy_ara, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                       fill = CAZy, label =Gene_ID ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 


ggplot(A_CAZy_ara, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                       fill = CAZy, label =gene_name ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 

A_CAZy_gala <-read.delim("./A_CAZy_gala.txt",row.names = 1)  

ggplot(A_CAZy_gala, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                        fill = CAZy, label = Gene_ID ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 
ggplot(A_CAZy_gala, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                        fill = CAZy, label =gene_name ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 
ggsave("果胶分解gala基因结构图_genename.pdf",height = 4.5, width =38) 

A_CAZy_rham <-read.delim("./A_CAZy_rham.txt",row.names = 1)  

ggplot(A_CAZy_rham, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                        fill = CAZy, label = Gene_ID ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 

ggplot(A_CAZy_rham, aes(xmin = Start_NEW, xmax = End_NEW, y = molecule, 
                        fill = CAZy, label =gene_name ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes() 

A_CAZy_rcsC <-read.delim("./rcsC.txt",row.names = 1)  

ggplot(A_CAZy_rcsC, aes(xmin = Start_NEW, xmax = End_NEW, y = strain, 
                        fill = cluster, label = Gene_ID ))+
  geom_gene_arrow(arrowhead_height = unit(3, "mm"), 
                  arrowhead_width = unit(1, "mm")) +
  geom_gene_label(align = "centre",
                  min.size = 0.1,grow = T) +
  theme_genes()+
  scale_fill_brewer(palette = "Set3")

b <- read.csv("./CAPSULE-THICKNESS.csv", row.names=1)
group_names <- c( "CON","BU1",  "BU2", "BU3","BU4", "BU5")
b$group <- factor( b$group, level=group_names )
ggplot(b,aes(group,capsule,fill=group))+
  geom_beeswarm(cex = 1.82, show.legend = FALSE,pch=21,size=1.85)+
  stat_summary(fun=mean,fun.ymin=mean,fun.ymax=mean,geom='crossbar',width=0.3,size=0.3,color='black')+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.25,size=0.25)+
  scale_fill_manual(values = c("grey","#5cc3e8","#5cc3e8","#5cc3e8","#e95f5c","#e95f5c"))+theme_classic()+
  scale_y_continuous(limits = c(0,250))+
  stat_compare_means(comparisons = list(c("CON","BU1"), c("BU1","BU4"), c("BU5","BU4")) ,
                     method = "wilcox.test",label = "p.signif")
a <- read.csv("BU1-5-IgA-coating.csv", row.names=1)
group_names <- c( "PECCON","PECBU1",  "PECBU2", "PECBU3","PECBU4", "PECBU5",
                  "NOPECCON", "NOPECBU1", "NOPECBU2", "NOPECBU3", "NOPECBU4","NOPECBU5")
a$group <- factor( a$group, level=group_names )
ggplot(a,aes(group,IgA_coating))+
  geom_bar(colour="black",stat="summary",fun=mean,position=position_dodge(0.6),width = 0.65,fill='white')+
  geom_jitter(aes(fill=group),stroke=0.1,alpha=1,size=2.5,
              width = 0.1, height = 0,pch=21)+
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

m3h.g3 <-read.delim("./OD-VALUE.txt",header = T,row.names= 1)
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
