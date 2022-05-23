rm(list=ls())
load("../Rdata/after_annoprobe_exp_groups.Rdata")
save(dat3,group3,file = "../Rdata/validation.Rdata")

rm(list=ls())
load("../Rdata/hub_genes.Rdata")
load("../Rdata/validation.Rdata")

hub_genes[12] = "FYB1"
exp = dat3[hub_genes,]


library(tidyverse)
dat = data.frame(t(exp))
dat = mutate(dat,group = group3)
p = list()
library(pacman)
pacman::p_load(tidyverse,ggpubr,rstatix,ggsci,ggsignif,reshape2)
for(i in 17:26){
  p[[i-16]] = ggplot(data = dat,aes_string(x = "group",y=colnames(dat)[i]))+
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

## 选出验证集中验证的差异基因
choose_hub_genes = hub_genes[c(10,11,16,17,23)]
choose_hub_genes
save(choose_hub_genes,file = "../Rdata/choose_hub_genes.Rdata")
