rm(list=ls())
options(stringsAsFactors = F)

load(file = "../Rdata/deg.Rdata")

library(stringr)
library(tidyverse)
###. 加ENTREZID列，用于富集分析（symbol转entrezid，然后inner_join）
library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(deg$symbol, 
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = org.Hs.eg.db)
dim(deg)
deg <- inner_join(deg,s2e,by=c("symbol"="SYMBOL"))
dim(deg)
length(unique(deg$symbol))
deg <- deg[!duplicated(deg$symbol),]


###1. GO富集分析
#(1)输入数据
gene_up = deg$ENTREZID[deg$change == 'up'] 
gene_down = deg$ENTREZID[deg$change == 'down'] 
gene_diff = c(gene_up,gene_down)

#(2)富集分析
#以下步骤耗时很长，设置了存在即跳过
f1 <- paste0("../Rdata/",gpl_number,"_GO.Rdata")
if(!file.exists(f1)){
  ego <- enrichGO(gene = gene_diff,
                  OrgDb= org.Hs.eg.db,
                  ont = "ALL",
                  readable = TRUE)
  ego_BP <- enrichGO(gene = gene_diff,
                     OrgDb= org.Hs.eg.db,
                     ont = "BP",
                     readable = TRUE)
  #ont参数：One of "BP", "MF", and "CC" subontologies, or "ALL" for all three.
  save(ego,ego_BP,file = f1)
  save(ego,file = f1)
}

load(f1)

#(3)可视化
#条带图
barplot(ego)
#气泡图
dotplot(ego_BP)

pdf("../figure/go.pdf",width=12, height=10)
dotplot(ego, split = "ONTOLOGY", font.size = 8, 
        showCategory = 15) + 
  facet_grid(ONTOLOGY ~ ., space = "free_y",scales = "free_y") + 
  scale_y_discrete(labels = function(x) str_wrap(x, width = 65))
dev.off()



geneList = deg$logFC
names(geneList)=deg$ENTREZID

#(3)展示top通路的共同基因，要放大看。
#Gene-Concept Network
cnetplot(ego,categorySize="pvalue", foldChange=geneList,colorEdge = TRUE)
pdf("../figure/gene_concept_network.pdf",width=12, height=10)
cnetplot(ego, showCategory = 5,foldChange=geneList, circular = TRUE, colorEdge = TRUE)
dev.off()

###2. KEGG富集分析
#(1) 输入数据
gene_up <-  deg[deg$change == 'up','ENTREZID'] 
gene_down <-  deg[deg$change == 'down','ENTREZID'] 
gene_diff <-  c(gene_up,gene_down)

#(2) kegg富集分析
f2 <- paste0("../Rdata/",gpl_number,"_KEGG.Rdata")
if(!file.exists(f2)){
  kk.up <- enrichKEGG(gene  = gene_up,
                      organism  = 'hsa',
                      pvalueCutoff = 0.05,
                      qvalueCutoff =0.05)

  kk.down <- enrichKEGG(gene   =  gene_down,
                        organism   = 'hsa',
                        pvalueCutoff = 0.05,
                        qvalueCutoff =0.05)
  
  kk.diff <- enrichKEGG(gene  = gene_diff,
                        organism   = 'hsa',
                        pvalueCutoff = 0.05,
                        qvalueCutoff =0.05)

  save(kk.diff,kk.down,kk.up,file = f2)
}
load(f2)
head(kk.up)[,1:6]
dotplot(kk.up)
ggsave('../figure/kk.up.dotplot.pdf')

head(kk.down)[,1:6]
dotplot(kk.down );ggsave('../figure/kk.down.dotplot.pdf')

head(kk.diff)[,1:6]
barplot(kk.diff,showCategory = 12)
ggsave('../figure/kk.diff.barplot.pdf')

#(3) 看看富集到了吗？https://mp.weixin.qq.com/s/NglawJgVgrMJ0QfD-YRBQg
table(kk.diff@result$p.adjust<0.05)
table(kk.up@result$p.adjust<0.05)
table(kk.down@result$p.adjust<0.05)

#(4) 双向富集图
# 富集分析所有图表默认都是用p.adjust,富集不到可以退而求其次用p值，在文中说明即可
down_kegg <- kk.down@result %>%
  filter(pvalue<0.05) %>% #筛选行
  mutate(group=-1) #新增列

up_kegg <- kk.up@result %>%
  filter(pvalue<0.05) %>%
  mutate(group=1)

source("kegg_plot_function.R")
g_kegg <- kegg_plot(up_kegg,down_kegg)
g_kegg

ggsave(g_kegg,filename = '../figure/kegg_up_down.pdf')

