#引用包
library(limma)
setwd("D:\\panCancer\\21.CIBERSORT")        #设置工作目录
pFilter=0.05

#读取目录下的文件
files=dir()
files=grep("^symbol.",files,value=T)

outTab=data.frame()
for(i in files){
	#读取文件
	CancerType=gsub("symbol\\.|\\.txt","",i)
	rt=read.table(files,sep="\t",header=T,check.names=F)
	
	#如果一个基因占了多行，取均值
	rt=as.matrix(rt)
	rownames(rt)=rt[,1]
	exp=rt[,2:ncol(rt)]
	dimnames=list(rownames(exp),colnames(exp))
	data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
	data=avereps(data)
	
	#删除正常，只保留肿瘤样品
	group=sapply(strsplit(colnames(data),"\\-"),"[",4)
	group=sapply(strsplit(group,""),"[",1)
	group=gsub("2","1",group)
	data=data[,group==0]
	data=data[rowMeans(data)>0,]

	#数据矫正
	v <-voom(data, plot = F, save.plot = F)
	out=v$E
	out=rbind(ID=colnames(out),out)
	write.table(out,file="uniq.symbol.txt",sep="\t",quote=F,col.names=F)        #输出文件
	
	#运行CIBERSORT，得到免疫细胞含量结果
	source("panCancer21.CIBERSORT.R")
	results=CIBERSORT("ref.txt", "uniq.symbol.txt", perm=100, QN=TRUE)

	#输出每个样品的打分
	immune=read.table("CIBERSORT-Results.txt",sep="\t",header=T,row.names=1,check.names=F)
	immune=immune[immune[,"P-value"]<pFilter,]
	immune=as.matrix(immune[,1:(ncol(immune)-3)])
	outTab=rbind(outTab,cbind(immune,CancerType))
	file.remove("CIBERSORT-Results.txt")
	file.remove("uniq.symbol.txt")
}
out=cbind(ID=row.names(outTab),outTab)
write.table(out,file="CIBERSORT.result.txt",sep="\t",quote=F,row.names=F)
rm(list = ls())
outpdf="barplot.pdf"

data <- read.table('CIBERSORT.result.txt',header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)

pdf(outpdf,height=10,width=25)
par(las=1,mar=c(8,4,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=60,xpd=T);text(a1,-0.02,colnames(data),adj=1,cex=0.6);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),col=col,pch=15,bty="n",cex=1.3)
dev.off()

library(corrplot)
rt=read.table("CIBERSORT.result.txt",sep="\t",header=T,row.names=1,check.names=F)
rt<-na.omit(rt)
rt<-rt[,colMeans(rt)>0]
pdf("corHeatmap.pdf",height=13,width=13)              #保存图片的文件名称
corrplot(corr=cor(rt),
         method = "color",
         order = "hclust",
         tl.col="black",
         addCoef.col = "black",
         number.cex = 0.8,
         col=colorRampPalette(c("blue", "white", "red"))(50),
)
dev.off()