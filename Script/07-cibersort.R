rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.13")

#函数下载：https://content.cruk.cam.ac.uk/fmlab/sivakumar2016/Cibersort.R
source("source.R")   #注释文件

#整理基因表达数据
load("../Rdata/exp_group.Rdata")
load("../Rdata/after_annoprobe_exp_groups.Rdata")
ciber_input <- 2^exp - 1
write.table(ciber_input, file = "../data/cibersort_input.txt",
            sep = "\t", row.names = T,col.names = NA,quote = F)

#CIBERSORT计算
sig_matrix <- "LM22.txt"   #注释文件名
mixture_file = '../data/cibersort_input.txt'   #表达数据文件名
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=1000, QN=TRUE)

if (F){
res = as.data.frame(res_cibersort)
res = res %>% mutate(group = sample_group$group)
table(res[res$`P-value`<=0.05,"group"])


# filter for analysis
res_ciber = res[ifelse(res$`P-value`<=0.05,TRUE,FALSE),]
group_ciber = factor(res_ciber$group,levels = c("control","AF"))
}
# res_sample = data.frame(sample = rownames(res_ciber),group = group_ciber)



save(res_cibersort,group, file = "../Rdata/cibersort.Rdata")   #保存结果

