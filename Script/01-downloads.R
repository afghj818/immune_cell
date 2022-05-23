rm(list=ls())
options(stringsAsFactors = F)

library(GEOquery)
library(stringr)
# 在GEO页面确认是表达芯片数据

#1. 下载数据：注意查看下载文件的大小，检查数据完整性，并保存下载数据至本地
if (F){
  gse_number1 <- "GSE41177"
  gse_number2 <- "GSE79768"
  gse_number3 <- "GSE115574"
  gse_number4 <- "GSE14975"
  eSet1 <- getGEO(gse_number1, 
                destdir = '.', 
                getGPL = F, 
                AnnotGPL = F)
  eSet2 <- getGEO(gse_number2, 
                destdir = '.', 
                getGPL = F, 
                AnnotGPL = F)
  eSet3 <- getGEO(gse_number3, 
                destdir = '.', 
                getGPL = F, 
                AnnotGPL = F)
  eSet4 <- getGEO(gse_number4, 
                destdir = '.', 
                getGPL = F, 
                AnnotGPL = F)

 # 载入和检查数据
  class(eSet1)
  length(eSet1)
  eSet1 = eSet1[[1]]

  class(eSet2)
  length(eSet2)
  eSet2 = eSet2[[1]]

  class(eSet3)
  length(eSet3)
  eSet3 = eSet3[[1]]

  class(eSet4)
  length(eSet4)
  eSet4 = eSet4[[1]]

  
  save(eSet1,eSet2,eSet3,eSet4,file = "../Rdata/eSet.Rdata")
}

rm(list = ls())
load("../Rdata/eSet.Rdata")
###1. 提取表达矩阵exp
exp1 <- exprs(eSet1)
dim(exp1)
exp2 <- exprs(eSet2)
dim(exp2)
exp3 <- exprs(eSet3)
dim(exp3)
exp4 <- exprs(eSet4)
dim(exp4)

#检查矩阵是否正常，如果是空的就会报错，空的和有负值的矩阵需要处理原始数据
#如果表达矩阵为空，大多数是转录组数据
exp1[1:4,1:4]
exp2[1:4,1:4]
exp3[1:4,1:4]
exp4[1:4,1:4]

#判断矩阵是否需要log处理
boxplot(exp1)

library(limma)
exp1 = normalizeBetweenArrays(exp1)
boxplot(exp1,las=2)

boxplot(exp2)
boxplot(exp3)
exp4 = log2(exp4 + 1)
boxplot(exp4)

###2. 提取临床信息metadata
pd1 <- pData(eSet1)
pd2 <- pData(eSet2)
pd3 <- pData(eSet3)
pd4 <- pData(eSet4)


## 取子集：取出对应的样本
cg1 <- str_detect(pd1$title,"LAA")
pd1 <- pd1[cg1,]
exp1 <- exp1[,cg1]
exp1 = data.frame(exp1)
# 考虑到样本表达量离群，去除相应的样本做后续分析：c("GSM1006247","GSM1006249")
# del_sample <- c("GSM1006247","GSM1006249")
# pd1 = pd1[!(colnames(exp1) %in% del_sample),]
# exp1 = exp1[,!(colnames(exp1) %in% del_sample)]

p = identical(rownames(pd1),colnames(exp1));p
if(!p) {exp1 = exp1[,match(rownames(pd1),colnames(exp1))]}

cg2 <- str_detect(pd2$title,"left")
pd2 <- pd2[cg2,]
exp2 <- exp2[,cg2]

p = identical(rownames(pd2),colnames(exp2));p
if(!p) {exp2 = exp2[,match(rownames(pd2),colnames(exp2))]}

cg3 <- str_detect(pd3$source_name_ch1,"left")
pd3 <- pd3[cg3,]
exp3 <- exp3[,cg3]

p = identical(rownames(pd3),colnames(exp3));p
if(!p) {exp3 = exp3[,match(rownames(pd3),colnames(exp3))]}

p = identical(rownames(pd4),colnames(exp4));p
if(!p) {exp4 = exp4[,match(rownames(pd4),colnames(exp4))]}


###3. 让exp列名与metadata的行名顺序完全一致,若不一致，调整顺序

###4. 提取芯片平台编号
gpl_number <- eSet4@annotation;gpl_number


###5. 获取分组信息
group1 <- ifelse(str_detect(pd1$title,"AF"),
                 "AF",
                 "control")
# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
group1 <-  factor(group1,levels = c("control","AF"))
group1
table(group1)

group2 <- ifelse(str_detect(pd2$title,"Fibrillation"),
                 "AF",
                 "control")
# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
group2 <-  factor(group2,levels = c("control","AF"))
group2
table(group2)

group3 <- ifelse(str_detect(pd3$source_name_ch1,"SR"),
                 "control","AF"
)
# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
group3 <-  factor(group3,levels = c("control","AF"))
group3
table(group3)

group4 <- ifelse(str_detect(pd4$source_name_ch1,"control"),
                 "control",
                 "AF")
# 需要把Group转换成因子，并设置参考水平，指定levels，对照组在前，处理组在后
group4 <-  factor(group4,levels = c("control","AF"))
group4
table(group4)



exp1 = data.frame(exp1)
exp2 = data.frame(exp2)
exp3 = data.frame(exp3)
exp4 = data.frame(exp4)


save(group1,group2,group3,group4,exp1,exp2,exp3,exp4,gpl_number,file = "../Rdata/before_annoprobe_exp_groups.Rdata")

