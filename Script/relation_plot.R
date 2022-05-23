rm(list=ls())

load("../Rdata/exp_group.Rdata")

load("../Rdata/filter_ciber_for_figures.Rdata")


exprSet = exp_ciber["NCF2",]
identical(colnames(exprSet),colnames(re))

dat = rbind(exprSet,re)
plot( as.numeric(dat["NCF2",]), as.numeric(dat["T cells follicular helper",]))
ggplot(data=dat,aes(x = dat["NCF2",],y = dat["B cells naive",])) +
  geom_point()
