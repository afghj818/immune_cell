### ---------------
###
### Create: Jianming Zeng
### Date: 2018-07-09 20:11:07
### Email: jmzeng1314@163.com
### Blog: http://www.bio-info-trainee.com/
### Forum:  http://www.biotrainee.com/thread-1376-1-1.html
### CAFS/SUSTC/Eli Lilly/University of Macau
### Update Log: 2018-07-09  First version
###

### Update: Xiaojie Liang
### Date: 2021-07-12 23:38:07
### Email: liangxiaojiecom@163.com
### ---------------
library(ggthemes)
library(ggplot2)
kegg_plot <- function(up.data,down.data){
  dat=rbind(up.data,down.data)
  colnames(dat)
  dat$pvalue = -log10(dat$pvalue)
  dat$pvalue=dat$pvalue*dat$group 
  
  dat=dat[order(dat$pvalue,decreasing = F),]
  
  gk_plot <- ggplot(dat,aes(reorder(Description, pvalue), y=pvalue)) +
    geom_bar(aes(fill=factor((pvalue>0)+1)),stat="identity", width=0.7, position=position_dodge(0.7)) +
    coord_flip() +
    scale_fill_manual(values=c("#0072B2", "#B20072"), guide="none") +
    labs(x="", y="" ) +
    theme_pander()  +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #axis.ticks.x = element_blank(),
          axis.line.x = element_line(size = 0.3, colour = "black"),#x轴连线
          axis.ticks.length.x = unit(-0.20, "cm"),#修改x轴刻度的高度，负号表示向上
          axis.text.x = element_text(margin = margin(t = 0.3, unit = "cm")),##线与数字不要重叠
          axis.ticks.x = element_line(colour = "black",size = 0.3) ,#修改x轴刻度的线                         
          axis.ticks.y = element_blank(),
          axis.text.y  = element_text(hjust=0),
          panel.background = element_rect(fill=NULL, colour = 'white')
    )
}


