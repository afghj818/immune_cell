rm(list=ls())
options(stringsAsFactors = F)

# 载入注释后的合并矩阵:exp
load("../Rdata/after_annoprobe_exp_groups.Rdata")

load("../Rdata/hub_genes.Rdata")

tr2 = function(gene_d,exp_d,group_d){
g = exp_d[gene_d,]
g = t(g)
g = as.data.frame(g)
g["group"] = group_d
colnames(g) = c("gene","group")
boxplot(gene~group,data = g)
w = wilcox.test(gene~group,data = g)
message("w: ",w[["p.value"]])
#t = t.test(gene~group,data = g)
#message("t: ",t[["p.value"]])
}

for (i in 1:10){
  gene = genes[i]
  
  tr2(gene,dat2,group2)
}



