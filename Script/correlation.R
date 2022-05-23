rm(list = ls())

load("../Rdata/diff_immune_cells.Rdata")
load("../Rdata/choose_hub_genes.Rdata")

exp_genes = exp[choose_hub_genes,]
identical(colnames(exp),colnames(exp_immune))

library(ggpubr)
library(stringr)
exp_genes = t(exp_genes)
exp_immune = t(exp_immune)
re = t(re)

library(Hmisc)#加载包
res2 <- rcorr(exp_genes,re)
res2
res_gene = list()
### 生成绘图数据框
for (i in 1:5){
res_gene[[i]] = data.frame(cor = res2$r[i,6:ncol(res2$r)],
                       pvalue = res2$P[i,6:ncol(res2$P)])
res_gene[[i]]["name"] = rownames(res_gene[[i]])
}
for (i in 1:5){
dat = res_gene[[i]]
p = ggdotchart(dat, x = "name", y = "cor",
           color = "pvalue",                                # Color by groups
            # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           add.params = list(color = "lightgray", size = 2), # Change segment color and size
           # Order by groups
           title = choose_hub_genes[i],
           dot.size = 6,                                 # Large dot size
           label = round(dat$cor,1),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 9,
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = theme_pubr()                        # ggplot2 theme
) + geom_hline(yintercept = 0, linetype = 2, color = "lightgray")
ggsave(p,filename = paste0("../figure/",choose_hub_genes[i],".pdf"))}


