rm(list = ls())

library(tinyarray)
library(dplyr)

### limma差异基因
load("../Rdata/deg.Rdata")
gene_up= deg[deg$change == 'up','symbol'] 
gene_down=deg[deg$change == 'down','symbol']
gene_diff = c(gene_up,gene_down)

### ggcna模块基因
gene_set1 = read.csv("../data/wgcna_genes.csv")
colnames(gene_set1) = c("num","gene")
gene_diff_wgcna = gene_set1$gene

### 取交集
choose_genes = intersect(gene_diff,gene_diff_wgcna)
library (VennDiagram)

venn.diagram(x= list(WGCNA = gene_diff_wgcna,Limma = gene_diff),
             filename = paste0("../figure/choose_genes",".tiff"),
             height = 800, width = 800,
             resolution =300,
             imagetype="tiff",
             col="transparent",
             fill=c("green","yellow"),
             alpha = 0.50,
             cex=0.45,
             cat.cex=0.45)


# 制作string的输入数据
write.table(choose_genes,
            file="../data/diffgene.txt",
            row.names = F,
            col.names = F,
            quote = F)

### 将diffgene.txt上传string网页，获得hub基因文件cyto.csv

cyto_genes = read.csv("../data/cyto.csv")
MCC = rownames(cyto_genes %>% 
                 arrange(desc(MCC)) %>% 
                 head(50))

DMNC = rownames(cyto_genes %>% 
                  arrange(desc(DMNC)) %>% 
                  head(50))

MNC = rownames(cyto_genes %>% 
                 arrange(desc(MNC)) %>% 
                 head(50)
)

Degree = rownames(cyto_genes %>% 
                    arrange(desc(Degree)) %>% 
                    head(50))

EPC = rownames(cyto_genes %>% 
                 arrange(desc(EPC)) %>% 
                 head(50))

hub_genes =intersect_all(MNC, MCC, DMNC, Degree,EPC)
library(UpSetR)
g = list(MCC = MCC,
         DMNC = DMNC,
         MNC = MNC,
         Degree = Degree,
         EPC = EPC)
upset(fromList(g), order.by = "freq",
      queries = list(list(query = intersects,
                                params = list("MCC","DMNC","MNC","Degree","EPC"),
                                      color = "darkred",
                                    active = T)))

venn.diagram(x= g,
             filename = paste0("../figure/genes",".tiff"),
             disable.logging = T,
             height = 1000, width = 1000,
             resolution =400,
             imagetype="tiff",
             col="transparent",
             fill=c("green","yellow","pink","blue","red"),
             alpha = 0.50,
             cex=0.45,
             cat.cex=0.5,cat.pos = c(-40,-30,120,100,20))

save(hub_genes,file = "../Rdata/hub_genes.Rdata")

