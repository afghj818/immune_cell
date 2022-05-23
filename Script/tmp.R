rm(list=ls())
dat = read.csv("../data/CIBERSORT.Output_Job5.csv")
load("../Rdata/exp_group.Rdata")
rownames(dat) = dat[,1]
dat =dat[,-1]
library(tidyverse)
identical(rownames(dat),sample_group$sample)
dat = dat %>% mutate(group = sample_group$group[1:32])

table(dat[dat$P.value <=0.05,"group"])
