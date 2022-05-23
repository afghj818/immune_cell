rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)


#可视化展示
rm(list=ls())
load("../Rdata/filter_ciber_for_figures.Rdata")
load("../Rdata/exp_group.Rdata")

### 此处生成表达矩阵exp, 免疫细胞丰度矩阵re, 分组信息group

index = as.character(group)
identical(colnames(exp),colnames(re))
exp = cbind(exp[,index == "AF"],exp[,index == "control"])
re = cbind(re[,index == "AF"],re[,index == "control"])
identical(colnames(exp),colnames(re))
group = factor(c(rep("AF",28),rep("control",14)),levels = c("control","AF"))

###1. 直方图
mycol <- ggplot2::alpha(rainbow(nrow(re)), 0.7)
dat <- t(re) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

p = ggplot(dat,aes(Sample,Proportion,fill = Cell_type)) + 
  geom_bar(stat = "identity") +
  labs(fill = "Cell Type",x = "",y = "Estiamted Proportion") + 
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom") + 
  scale_y_continuous(expand = c(0.01,0)) +
  scale_fill_manual(values = mycol)

ggsave(p, filename = "../figure/barplot_ciber.pdf", width = 10, height = 8,dpi = 400)


###2.热图

anno_col = data.frame(
  group = group)
rownames(anno_col) = colnames(exp)

anno_color = list(Type = c(control = "#00dae0", AF = "#ff9289"))  


pdf(file = "../figure/heatmap_ciber.pdf")
pheatmap::pheatmap(re,
         show_colnames = F,
         annotation_col = anno_col, 
         cluster_cols = F,
         color = colorRampPalette(c("#60c049", "red"))(20),
         annotation_colors = anno_color)
dev.off()


###3. PCA图
dat <- as.data.frame(t(re))
library(FactoMineR)
library(factoextra) 
dat.pca <- PCA(dat, graph = FALSE)
pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = group, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups"
)

pdf("../figure/pca_ciber.pdf")
pca_plot
dev.off()

###4. 相关性热图
M <- round(cor(t(re)),2) # 计算相关性矩阵并保留两位小数
library(corrplot)
pdf("../figure/corHeatmap_ciber.pdf",height=13,width=13)              #保存图片的文件名称
corrplot(corr=M,
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()


library(tidyverse)
library(RColorBrewer)
###5. 免疫细胞箱线图
test_name = colnames(exp)

dat <- t(re) %>% as.data.frame() %>%
  rownames_to_column("Sample") %>% 
  gather(key = Cell_type,value = Proportion,-Sample)

mypalette <- colorRampPalette(brewer.pal(8,"Set1"))

dat$Group = ifelse(dat$Sample %in% test_name[1:28],"AF","control")
library(ggpubr)
p = ggplot(dat,aes(Cell_type,Proportion,fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black") + 
  theme_bw() + 
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle=80,vjust = 0.5))+
  scale_fill_manual(values = mypalette(22)[c(1,6)])+ stat_compare_means(aes(group = Group,label = ..p.signif..),method = "wilcox.test")
ggsave(p, filename = "../figure/boxplot_ciber.pdf",width = 10,height = 8,dpi = 300)

diff_immune_cells = rownames(re)[c(1,2,17,21,4,9,8)]
exp_immune = re[diff_immune_cells,]
save(diff_immune_cells,group,exp,re,exp_immune,file = "../Rdata/diff_immune_cells.Rdata")

### 小提琴图
library(vioplot)                                                    #引用包

normal=5                                                            #正常样品数目
tumor=22                                                         #AF样品数目
rt = t(re)
# rt表达矩阵前tumor行为tumor信息，后normal行为normal信息

pdf("../figure/vioplot_ciber.pdf",height=8,width=15)              #保存图片的文件名称
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")

#对每个免疫细胞循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
for(i in 1:ncol(rt)){
  tumorData=rt[1:tumor,i]
  normalData=rt[(tumor+1):(normal+tumor),i]
  vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(tumorData,normalData)
  p=round(wilcoxTest$p.value,3)
  mx=max(c(normalData,tumorData))
  lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
  text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
  text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
}
dev.off()

