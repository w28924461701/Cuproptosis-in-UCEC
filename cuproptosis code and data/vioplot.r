par(mfrow=c(1,4))

#循环，绘制vioplot，正常用蓝色表示，肿瘤用红色表示
library(vioplot)
tt$risk<-as.factor(tt$risk)
for(i in c(4:7)){
vioplot(tt[,i]~tt$risk,col = c('red','blue'),
xlab=c("high", "low"),
ylab =colnames(tt)[i],lwd =2,cex=2,cex.axis =1.5,font.axis =1,font.lab=2)

}

wilcox.test(StromalScore ~ risk, data = tt)