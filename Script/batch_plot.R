rm(list=ls())
options(stringsAsFactors = F)

load("../Rdata/batch_before_after.Rdata")



library(tidyr)
library(tibble)
library(dplyr)
batch <- paste0("batch",rep(c(1,2,4),c(19,13,10)))
colors = c(rep("red",19),rep("blue",13),rep("purple",10))

#### before
dat = t(exprSet) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = batch) 

colnames(dat)[2:(ncol(dat)-1)] <- paste0("gene_",colnames(dat)[2:(ncol(dat)-1)])
# 宽数据变为长数据
pdat = dat %>% 
  pivot_longer(cols = starts_with("gene"),
               names_to = "gene",
               values_to = "exp")

pdf("../figure/before_batch.pdf")
boxplot(exp~rowname,data = pdat,col = colors)
dev.off()


#### after
dat_batch = t(exprSet_batch) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = batch)

colnames(dat_batch)[2:(ncol(dat_batch)-1)] <- paste0("gene_",colnames(dat_batch)[2:(ncol(dat_batch)-1)])
# 宽数据变为长数据
pdat_batch = dat_batch %>% 
  pivot_longer(cols = starts_with("gene"),
               names_to = "gene",
               values_to = "exp")

# 画图
pdf("../figure/after_batch.pdf")
boxplot(exp~rowname,data = pdat_batch,col = colors)
dev.off()
