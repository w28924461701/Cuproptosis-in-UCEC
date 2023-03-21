#引用包
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
library(tidyverse)

coxPfilter=0.05        #单因素cox方法显著性的过滤标准
setwd("D:\\__easyHelper__\\148cuproptosis\\14.model")      #设置工作目录

#读取输入文件
rt=read.table("expTime.txt", header=T, sep="\t", check.names=F, row.names=1)
rt$days_to_last_followup[rt$days_to_last_followup<=0]=1
rt$days_to_last_followup=rt$days_to_last_followup/365
rt[,21:ncol(rt)]=log2(rt[,21:ncol(rt)]+1)
rt1<-rt
rt<-rt[,c(2,9,21:ncol(rt))]

fix(rt)
rt2<-rt[,c(3:ncol(rt))]

rt2=rt2[,colMeans(rt2)>0.1]
#data<-as.data.frame(data)
rt2<-cbind(rt[,c(1,2)],rt2)
rt<-rt2
############定义森林图函数############
bioForest=function(coxFile=null,forestFile=null,forestCol=null){
	#读取输入文件
	rt <- read.table(coxFile,header=T,sep="\t",row.names=1,check.names=F)
	gene <- rownames(rt)
	hr <- sprintf("%.3f",rt$"HR")
	hrLow  <- sprintf("%.3f",rt$"HR.95L")
	hrHigh <- sprintf("%.3f",rt$"HR.95H")
	Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
	pVal <- ifelse(rt$pvalue<0.001, "<0.001", sprintf("%.3f", rt$pvalue))
		
	#输出图形
	pdf(file=forestFile, width=7, height=6)
	n <- nrow(rt)
	nRow <- n+1
	ylim <- c(1,nRow)
	layout(matrix(c(1,2),nc=2),width=c(3,2.5))
		
	#绘制森林图左边的临床信息
	xlim = c(0,3)
	par(mar=c(4,2.5,2,1))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
	text.cex=0.8
	text(0,n:1,gene,adj=0,cex=text.cex)
	text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,adj=1)
	text(3,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,adj=1,)
		
	#绘制森林图
	par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
	LOGindex = 10 
	hrLow = log(as.numeric(hrLow),LOGindex)
	hrHigh = log(as.numeric(hrHigh),LOGindex)
	hr = log(as.numeric(hr),LOGindex)
	xlim = c(floor(min(hrLow,hrHigh)),ceiling(max(hrLow,hrHigh)))
	plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
	arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
	abline(v=log(1,LOGindex),col="black",lty=2,lwd=2)
	boxcolor = ifelse(as.numeric(hr) > log(1,LOGindex), forestCol[1],forestCol[2])
	points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.3)
	a1 = axis(1,labels=F,tick=F)
	axis(1,a1,10^a1)
	dev.off()
}
############绘制森林图函数############

#对数据进行分组，构建模型
n=1  #分组的数目
for(i in 1:n){
	#############对数据进行分组#############
	inTrain<-createDataPartition(y=rt[,2], p=0.7, list=F)
	train<-rt
	test<-rt[inTrain,]
	trainOut=cbind(id=row.names(train),train)
	testOut=cbind(id=row.names(test),test)
	
	#单因素cox分析
	outUniTab=data.frame()
	sigGenes=c("futime","fustat")
	for(i in colnames(train[,3:ncol(train)])){
		#cox分析
		cox <- coxph(Surv(futime, fustat) ~ train[,i], data = train)
	
		coxSummary = summary(cox)
		coxP=coxSummary$coefficients[,"Pr(>|z|)"] 

		#保留显著性基因
		if(coxP<0.01){
		    sigGenes=c(sigGenes,i)
			outUniTab=rbind(outUniTab,
				         cbind(id=i,
				         HR=coxSummary$conf.int[,"exp(coef)"],
				         HR.95L=coxSummary$conf.int[,"lower .95"],
				         HR.95H=coxSummary$conf.int[,"upper .95"],
				         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
				         )
		}
	}
	
	#outUniTab<-filter(outUniTab,pvalue<0.001)
	
	
	uniSigExp=train[,sigGenes]
	uniSigExpOut=cbind(id=row.names(uniSigExp),uniSigExp)
	if(ncol(uniSigExp)<10){next}
	
	#

	
	#lasso回归 LASSO—交叉验证选lambda
	x=as.matrix(uniSigExp[,c(3:ncol(uniSigExp))])
	
	y=data.matrix(Surv(uniSigExp$futime,uniSigExp$fustat))
	fit <- glmnet(x, y, family = "cox", maxit = 1000)
#	fit1 <- glmnet(x, y, family = "cox", maxit = 1000,alpha = 0) 岭回归
	pdf("lambda.pdf")
	plot(fit, xvar = "lambda", label = TRUE)
	dev.off()
	cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
	pdf("cvfit.pdf")
	plot(cvfit)
	abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
	dev.off()

	
	#LASSO—建模与验证
	
	# 根据最佳lambda构建岭回归模型

	coeff_lasso <- predict(fit, s = cvfit$lambda.min, type = 'coefficients')
	coeff_lasso
	
	# 模型评估
	newx1 = as.matrix(test[,sigGenes])
	
	pred_lasso <- predict(fit, s =cvfit$lambda.min, newx = as.matrix(newx1[,c(3:ncol(newx1))]))
	RMSE <- sqrt(mean((test$futime
	                     -pred_lasso)**2))
	RMSE
	
	
	
	
	coef <- coef(fit, s = cvfit$lambda.min)
	index <- which(coef != 0)
	actCoef <- coef[index]
	lassoGene=row.names(coef)[index]
	lassoSigExp=uniSigExp[,c("futime", "fustat", lassoGene)]
	lassoSigExpOut=cbind(id=row.names(lassoSigExp), lassoSigExp)
	geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
	if(nrow(geneCoef)<2){next}
	
	#############构建COX模型#############
	multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
	multiCox=step(multiCox, direction = "both")
	multiCoxSum=summary(multiCox)
	
	#输出模型的公式
	outMultiTab=data.frame()
	outMultiTab=cbind(
		               coef=multiCoxSum$coefficients[,"coef"],
		               HR=multiCoxSum$conf.int[,"exp(coef)"],
		               HR.95L=multiCoxSum$conf.int[,"lower .95"],
		               HR.95H=multiCoxSum$conf.int[,"upper .95"],
		               pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
	outMultiTab=cbind(id=row.names(outMultiTab),outMultiTab)
	outMultiTab=outMultiTab[,1:2]
	
	#输出train组风险文件
	riskScore=predict(multiCox,type="risk",newdata=train)      #利用train得到模型预测train样品风险
	coxGene=rownames(multiCoxSum$coefficients)
	coxGene=gsub("`","",coxGene)
	outCol=c("futime","fustat",coxGene)
	medianTrainRisk=median(riskScore)
	risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
	trainRiskOut=cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk))
		
	#输出test组风险文件
	riskScoreTest=predict(multiCox,type="risk",newdata=test)     #利用train得到模型预测test样品风险
	riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
	testRiskOut=cbind(id=rownames(cbind(test[,outCol],riskScoreTest,riskTest)),cbind(test[,outCol],riskScore=riskScoreTest,risk=riskTest))
	
	#比较高低风险组的生存差异，得到差异的pvalue	
	diff=survdiff(Surv(futime, fustat) ~risk,data = train)
	pValue=1-pchisq(diff$chisq, df=1)
	diffTest=survdiff(Surv(futime, fustat) ~riskTest,data = test)
	pValueTest=1-pchisq(diffTest$chisq, df=1)
	
	
	
	#ROC曲线下面积
	predictTime=1    #预测时间
	roc=timeROC(T=train$futime, delta=train$fustat,
	            marker=riskScore, cause=1,
	            times=c(predictTime), ROC=TRUE)
	
	rocTest=timeROC(T=test$futime, delta=test$fustat,
	                marker=riskScoreTest, cause=1,
	                times=c(predictTime), ROC=TRUE)	
	
	plot(roc,time=1)   
	plot(rocTest,time=1,add=TRUE,col="blue")  
	legend("bottomright",c("train group","test group"),col=c("red","blue"),lty=1,lwd=2)
	

	
	
	
#if((pValue<0.05) & (roc$AUC[2]>0.68) & (pValueTest<0.05) & (rocTest$AUC[2]>0.68)){
		#输出分组结果
		write.table(trainOut,file="data.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testOut,file="data.test.txt",sep="\t",quote=F,row.names=F)
		#输出单因素结果
		write.table(outUniTab,file="uni.trainCox.txt",sep="\t",row.names=F,quote=F)
		write.table(uniSigExpOut,file="uni.SigExp.txt",sep="\t",row.names=F,quote=F)
		bioForest(coxFile="uni.trainCox.txt",forestFile="uni.foreast.pdf",forestCol=c("red","green"))
	    #lasso结果
	    write.table(lassoSigExpOut,file="lasso.SigExp.txt",sep="\t",row.names=F,quote=F)
		pdf("lasso.lambda.pdf")
		plot(fit, xvar = "lambda", label = TRUE)
		dev.off()
		pdf("lasso.cvfit.pdf")
		plot(cvfit)
		abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)), lty="dashed")
		dev.off()
	    #输出多因素结果
		write.table(outMultiTab,file="multiCox1.txt",sep="\t",row.names=F,quote=F)
		write.table(trainRiskOut,file="risk.train.txt",sep="\t",quote=F,row.names=F)
		write.table(testRiskOut,file="risk.test.txt",sep="\t",quote=F,row.names=F)
	#	bioForest(coxFile="multiCox1.txt",forestFile="Multi.foreast.pdf",forestCol=c("red","green"))
		#所有样品的风险值
		allRiskOut=rbind(trainRiskOut, testRiskOut)
		write.table(allRiskOut,file="risk.all.txt",sep="\t",quote=F,row.names=F)
		break
	}
}