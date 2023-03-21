p1<-ggplot(data=data, aes(x=risk, y=TMB,fill=risk)) +
  geom_violin(trim=FALSE)+
    geom_boxplot(width=0.1, fill="white")+
   labs(y = "tumor mutational burden")+
  ylim(0,30)+
  theme_classic()+theme(axis.text = element_text(size = 12))+
  theme(legend.position="top")

wilcox.test(TMB~risk,data)

library(dplyr)
library(ggpubr)
tideFile="MSI.txt"          #MSI的打分文件
riskFile="risk.all.txt"      #风险文件
setwd("D:\\STAG2 博士后项目一\\148cuproptosis\\34.MSI")     #设置工作目录

#读取TIDE数据
MSI=read.table(tideFile, header=T, sep="\t", check.names=F, row.names=1)
MSI1<-filter(MSI,CancerType%in%c('UCEC'))
rownames(MSI1)<-substr(rownames(MSI1),1,12)
#读取风险数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(MSI1), row.names(risk))
MSI1=MSI1[sameSample, , drop=F]
risk=risk[sameSample, "risk", drop=F]
data=cbind(tide, risk)

p2<-ggplot(data=data, aes(x=risk, y=MSI,fill=risk)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(y = "microsatellite instability")+
 # ylim(0,30)+
  theme_classic()+theme(axis.text = element_text(size = 12))+
  theme(legend.position="top")

wilcox.test(MSI~risk,data)

library(cowplot)
plot_grid(p1,p2, ncol = 2)
