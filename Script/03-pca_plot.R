rm(list=ls())
options(stringsAsFactors = F)

# 载入注释后的合并矩阵:exp
load("../Rdata/exp_group.Rdata")

# 载入注释后的4组表达矩阵:dat1,dat2,dat3,dat4
load("../Rdata/after_annoprobe_exp_groups.Rdata")

### 绘制PCA图
library(FactoMineR)
library(factoextra) 

draw_pca = function(pca_exp,pca_group,gse_number){
 dat_pca <- as.data.frame(t(pca_exp))
 dat.pca <- PCA(dat_pca, graph = FALSE)
 pca_plot <- fviz_pca_ind(dat.pca,
                         geom.ind = "point", # show points only (nbut not "text")
                         col.ind = pca_group, # color by groups
                         palette = c("#00AFBB", "#E7B800"),
                         addEllipses = TRUE, # Concentration ellipses
                         legend.title = "Groups",
                         title = gse_number)
 return(pca_plot)
}

draw_pca(dat1,group1,"GSE41177")
draw_pca(dat2,group2,"GSE79768")
draw_pca(dat3,group3,"GSE115574")
draw_pca(dat4,group4,"GSE14975")
draw_pca(exp,group,"merge")

