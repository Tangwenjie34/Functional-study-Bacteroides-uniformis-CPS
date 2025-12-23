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
library(rstatix)

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
kruskal_result <- kruskal.test(capsule ~ group, data = b)
print(kruskal_result)
dunn_result <- b %>%
  dunn_test(capsule ~ group, 
            p.adjust.method = "bonferroni")  
print(dunn_result)

c <- read.csv("./CAPSULE-THICKNESS_1.csv", row.names=1)
group_names <- c( "WT","rcsC")
c$group <- factor( c$group, level=group_names )
ggplot(c,aes(group,capsule,fill=group))+
  geom_beeswarm(cex = 2.5, show.legend = FALSE,pch=21,size=3)+
  stat_summary(fun=mean,fun.ymin=mean,fun.ymax=mean,geom='crossbar',width=0.3,size=0.5,color='black')+
  stat_summary(geom = "errorbar",fun.data = 'mean_se', width = 0.2)+
  scale_fill_manual(values = c("grey","#5cc3e8"))+theme_classic()+
  scale_y_continuous(limits = c(0,200))
wilcox_result <- wilcox.test(capsule ~ group, data = c)
print(wilcox_result)
