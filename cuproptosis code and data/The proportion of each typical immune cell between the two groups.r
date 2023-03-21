rownames(rt)<-substr(rownames(rt),1,12)
risk=read.table('risk.all.txt', header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(rownames(rt), row.names(risk))
rt=rt[sameSample,]
risk=risk[sameSample,]
library(vioplot)
#提取低风险组和高风险组的样品
riskLow=risk[risk$risk=="low",]
riskHigh=risk[risk$risk=="high",]
dataLow=rt[row.names(riskLow),]
dataHigh=rt[row.names(riskHigh),]

pdf("vioplot.pdf",height=8,width=15) #保存图片的文件名称
par(las=1,mar=c(10,6,3,3))
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
xlim=c(0,63),ylim=c(min(rt),max(rt)+0.02),
main="",xlab="", ylab="Fraction",
pch=21,
col="white",
xaxt="n")

#对每个免疫细胞循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
for(i in 1:ncol(rt)){
normalData=dataLow[,i]
tumorData=dataHigh[,i]
vioplot(normalData,at=3*(i-1),lty=1,add = T,col = 'blue')
vioplot(tumorData,at=3*(i-1)+1,lty=1,add = T,col = 'red')
wilcoxTest=wilcox.test(normalData,tumorData)
p=round(wilcoxTest$p.value,3)
mx=max(c(normalData,tumorData))
lines(c(x=3*(i-1)+0.2,x=3*(i-1)+0.8),c(mx,mx))
text(x=3*(i-1)+0.5,y=mx+0.02,labels=ifelse(p<0.001,paste0("p<0.001"),paste0("p=",p)),cex = 0.8)
text(seq(1,64,3),-0.05,xpd = NA,labels=colnames(rt),cex = 1,srt = 45,pos=2)
legend("topleft", legend=c("low risk", "high risk"),
fill=c("blue", "red"), cex = 1)
}
dev.off()