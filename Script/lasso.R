rm(list=ls())

load("../Rdata/exp_group.Rdata")
load("../Rdata/deg.Rdata")
gene_up = deg$symbol[deg$change == 'up'] 
gene_down = deg$symbol[deg$change == 'down'] 
gene_diff = c(gene_up,gene_down)

exprSet = exp[gene_diff,]
exprSet = 2^exprSet - 1

x=t(exprSet)
y = ifelse(group == "control",0,1)
library(glmnet)


set.seed(109000)
cv_fit <- cv.glmnet(x=x, y=y)
plot(cv_fit)

fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")

model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)

choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]

length(choose_gene_min)


save(choose_gene_min,file = paste0("../Rdata/lasso_choose_gene_min.Rdata"))