rm(list=ls())
options(stringsAsFactors = F)

load("../Rdata/deg.Rdata")

###1. 火山图
library(dplyr)
library(ggplot2)
#(1) 火山图
p <- ggplot(data = deg, 
            aes(x = logFC, 
                y = -log10(P.Value))) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
pdf("../figure/vocalno.pdf")
p
dev.off()

##(2) 火山图基因注释
for_label <- deg %>% 
  filter(symbol %in% c("HADHA","LRRFIP1"))

volcano_plot <- p +
  geom_point(size = 3, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black"
  )
pdf("./figures/vocalno.pdf")
volcano_plot
dev.off()

###2. 差异基因热图
if(F){
  #全部差异基因
  cg <-  deg$symbol[deg$change !="stable"]
  length(cg)
}else{
  #取前20上调和前20下调
  library(dplyr)
  dat <-  deg %>%
    filter(change!="stable") %>%
    arrange(logFC)
  cg <-  c(head(dat$symbol,20),
         tail(dat$symbol,20))
}


#差异基因热图
library(pheatmap)
exp_arrange = cbind(exp[,group == "AF"],exp[,group != "AF"])
group_arrange = factor(c(rep("AF",28),rep("control",14)),levels = c("control","AF"))

# 取出差异基因的表达矩阵diff_exp
diff_exp <- exp_arrange[cg,]
dim(diff_exp)

annotation_col=data.frame(group=group_arrange)
rownames(annotation_col)=colnames(diff_exp) 
heatmap_plot <- pheatmap(as.matrix(diff_exp),show_colnames =F,
                         scale = "row",
                         cluster_cols = F, 
                         
                         annotation_col=annotation_col,
                         breaks = seq(-3,3,length.out = 100)
) 
ggsave(heatmap_plot,filename = "../figure/pheatmap.pdf",width = 10,height = 16,dpi=300)

library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("#2fa1dd", "white", "#f87669"))
top_annotation = HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = c("#f87669","#2fa1dd")),
                       labels = c("AF","control"),
                       labels_gp = gpar(col = "white", fontsize = 12)))

m = Heatmap(t(scale(t(as.matrix(diff_exp)))),name = " ",
            col = col_fun,
            top_annotation = top_annotation,
            column_split = group,
            show_heatmap_legend = T,
            border = F,
            show_column_names = F,
            show_row_names = F,
            column_title = NULL)
png("../figure/heatmap.png")
m
dev.off()

###3. 感兴趣基因的相关性----
library(corrplot)
g = c("TIMM50", "GOLGA2", "CLIP1", "DDX6", "CTSZ", "P2RY12", "LTF", 
      "HBD", "F3", "ZFP36L2")
M = cor(t(exp[g,]))
pheatmap(M)

library(paletteer)
my_color = rev(paletteer_d("RColorBrewer::RdYlBu"))
my_color = colorRampPalette(my_color)(10)
corrplot(M, type="upper",
         method="pie",
         order="hclust", 
         col=my_color,
         tl.col="black", 
         tl.srt=45)

