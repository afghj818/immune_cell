#可视化展示
rm(list=ls())
load("../Rdata/cibersort.Rdata")
load("../Rdata/exp_group.Rdata")

res_ciber = res_cibersort
### 此处生成表达矩阵exp_ciber, 免疫细胞丰度矩阵re, 分组信息group_ciber
res_ciber <- res_ciber[, 1:22]   #取前22列为细胞丰度数据
ciber.res <- res_ciber[,colSums(res_ciber) > 0]   #去除丰度全为0的细胞
re = as.data.frame(t(ciber.res))
if (F){
cg_filter = colnames(exp) %in% rownames(res_ciber)
exp_ciber = exp[,cg_filter]}

save(res_ciber,re,group,file = "../Rdata/filter_ciber_for_figures.Rdata")

