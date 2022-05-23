rm(list=ls())
options(stringsAsFactors = F)
library(dplyr)

load("../Rdata/exp_group.Rdata")
#检查数据：表达矩阵exp，分组信息group
table(group)
exp[1:4,1:4]

###1. 差异分析：需要表达矩阵exp和分组信息group，不需要改
library(limma)
design <- model.matrix(~group)
fit <- lmFit(exp,design)
fit <- eBayes(fit)
deg <- topTable(fit,coef=2,number = Inf)

###2. 加change列,标记上下调基因
logFC_t=log2(1.5)
P.Value_t = 0.05
k1 = (deg$P.Value < P.Value_t)&(deg$logFC < -logFC_t)
k2 = (deg$P.Value < P.Value_t)&(deg$logFC > logFC_t)
deg <- mutate(deg,change = ifelse(k1,"down",ifelse(k2,"up","stable")))
deg$symbol=rownames(deg)
table(deg$change)

save(group,exp,deg,logFC_t,P.Value_t,gpl_number,file = "../Rdata/deg.Rdata")

