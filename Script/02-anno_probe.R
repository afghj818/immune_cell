rm(list=ls())
options(stringsAsFactors = F)

load("../Rdata/before_annoprobe_exp_groups.Rdata")

library(tidyverse)
together <- function(a,b,c){
  cg = a %>% 
    intersect(b) %>% 
    intersect(c) 
  return(cg)
}

cg=together(rownames(exp1),rownames(exp2),rownames(exp4))
exprSet <- cbind(exp1[cg,],exp2[cg,],exp4[cg,])
boxplot(exprSet)

group = c(group1,group2,group4)

batch <- paste0("batch",rep(c(1,2,4),c(19,13,10)))
mod <- model.matrix(~group)
exprSet_batch <- sva::ComBat(dat=exprSet,batch=batch,mod=mod)
boxplot(exprSet_batch)



### 探针注释
#通过tinyarray包获取探针注释代码,其他方式见find_anno.R
if (F) {
library(tinyarray)
find_anno(gpl_number)
ids <- AnnoProbe::idmap('GPL570')
head(ids)    # ids为探针和基因symbol的对应数据框
save(ids,file = "../Rdata/ids.Rdata")
}

load("../Rdata/ids.Rdata")
if (F){
#探针注释去重(一个基因对应多个探针，取最大表达量探针)，其他方式详见move_duplicated.R
id_dup = function(exprSet_mat,ids_mat){
  colnames(ids_mat)=c('probe_id','symbol')  
  ids_mat <- ids_mat[ids_mat$probe_id %in% rownames(exprSet_mat),]
  exprSet_mat <- exprSet_mat[ids_mat$probe_id,]
  message(identical(rownames(exprSet_mat),ids_mat$probe_id))
  ids_mat$median <- apply(exprSet_mat,1,median) 
  ids_mat <- ids_mat[order(ids_mat$symbol,ids_mat$median,decreasing = T),]
  ids_mat <- ids_mat[!duplicated(ids_mat$symbol),]
  exprSet_mat <- exprSet_mat[ids_mat$probe_id,]
  rownames(exprSet_mat) <- ids_mat$symbol
  return(exprSet_mat)
}


dat1 = id_dup(exp1,ids)
dat2 = id_dup(exp2,ids)
dat3 = id_dup(exp3,ids)
dat4 = id_dup(exp4,ids)


save(dat1,dat2,dat3,dat4,group1,group2,group3,group4,gpl_number,file = "../Rdata/after_annoprobe_exp_groups.Rdata")
}


load("../Rdata/after_annoprobe_exp_groups.Rdata")

exprSet = data.frame(exprSet)
exprSet_batch = data.frame(exprSet_batch)

exprSet = id_dup(exprSet,ids)
exprSet_batch = id_dup(exprSet_batch,ids)


save(exprSet,exprSet_batch,group,gpl_number,file = "../Rdata/batch_before_after.Rdata")

exp = exprSet_batch
sample_group = data.frame(sample = colnames(exp),group = group)

save(exp,group,sample_group,gpl_number,file = "../Rdata/exp_group.Rdata")

