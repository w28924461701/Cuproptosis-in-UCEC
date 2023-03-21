library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)

setwd("D:\\panCancer\\25.GSEA")      #设置工作目录
gene="BRCA1"                               #基因名称
gmtFile="c5.all.v7.1.symbols.gmt"          #基因集文件

#读入gmt文件
gmt=read.gmt(gmtFile)

#获取目录下的所有肿瘤数据文件
files=dir()
files=grep("^symbol.",files,value=T)

for(i in files){
	#读取肿瘤数据文件
	rt=read.table(i,sep="\t",header=T,check.names=F)
	CancerType=gsub("^symbol\\.|\\.txt$","",i)
	
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	#删除正常样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	data=t(data)
	rownames(data)=gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
	data=avereps(data)
	data=t(data)
	data=data[rowMeans(data)>0.05,]
	
	#按风险将样品分成高低表达两组

	
	#读取risk文件
	risk=read.table("risk.all.txt" , header=T, sep="\t", check.names=F, row.names=1)
	sameSample=intersect(colnames(data), row.names(risk))
	data=data[,sameSample]
	risk=risk[sameSample,]
	
	riskLow=risk[risk$risk=="low",]
	riskHigh=risk[risk$risk=="high",]
	dataL=data[,row.names(riskLow)]
	dataH=data[,row.names(riskHigh)]
	
	meanL=rowMeans(dataL)
	meanH=rowMeans(dataH)
	meanL[meanL<0.00001]=0.00001
	meanH[meanH<0.00001]=0.00001
	logFC=log2(meanH/meanL)
	logFC=sort(logFC,decreasing=T)
	
	
	
    #富集分析
    kk=GSEA(logFC,TERM2GENE=gmt, nPerm=100,pvalueCutoff = 1)

library(GseaVis)

terms <-  kkTab$ID[c(1,2,11,25,26,27)]


# plot
lapply(terms, function(x){
  gseaNb(object = kk,
         geneSetID =x,
         addPval = T,
         pvalX = 0.75,pvalY = 0.75,
         pCol = 'black',
         pHjust = 0)
}) -> gseaList

# combine
cowplot::plot_grid(plotlist = gseaList,ncol = 3,align = 'hv')