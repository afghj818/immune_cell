rm(list = ls())

dat <- rio::import("../data/GSE128188_CountsV2_EdgeR.txt",
                   format = "\t")

rownames(dat) = dat[,1]
exp = dat[,-1]

exp = exp[,seq(2,20,2)]

group = c(rep("SR",5),rep("AF",5))
group = factor(group,levels = c("SR","AF"))

exp_cpm = log2(edgeR::cpm(exp)+1)




## exp为表达矩阵，group为分组信息
save(exp,exp_cpm,group,file = "../Rdata/validation_counts.Rdata")

rm(list = ls())
load("../Rdata/validation_counts.Rdata")
load("../Rdata/hub_genes.Rdata")

exprSet = exp_cpm[genes,]
exprSet = as.data.frame(exprSet)

library(tidyverse)

dat = data.frame(t(exprSet))
dat = mutate(dat,group = group)
p = list()
library(pacman)
pacman::p_load(tidyverse,ggpubr,rstatix,ggsci,ggsignif,reshape2)
for(i in 1:(ncol(dat)-1)){
  p[[i]] = ggplot(data = dat,aes_string(x = "group",y=colnames(dat)[i]))+
    geom_boxplot(aes(color = group))+
    geom_jitter(aes(color = group))+
    theme_bw() +
    geom_signif(comparisons = list(c("control", "AF")),
                map_signif_level=T,
                textsize=4,test=wilcox.test,step_increase=0.2)
}
library(patchwork)
wrap_plots(p,nrow = 2,guides = "collect",widths = 4,
           heights = 10)


library(tidyr)
library(tibble)
library(dplyr)
# library(tidyverse)
exp = exprSet
dat = t(exp) %>% 
  as.data.frame() %>% 
  rownames_to_column() %>% 
  mutate(group = group)

colnames(dat)[2:(ncol(dat)-1)] = paste0("gene_",colnames(dat)[2:(ncol(dat)-1)])
# 宽数据变为长数据
pdat = dat %>% 
  pivot_longer(cols = starts_with("gene"),
               names_to = "gene",
               values_to = "exp")

# 画图
library(ggplot2)
p = ggplot(pdat,aes(gene,exp))+
  geom_boxplot(aes(fill = group))+
  geom_jitter()+
  theme_bw()
p
p + facet_wrap(~gene,scales = "free")

