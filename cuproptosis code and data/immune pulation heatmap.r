library(RColorBrewer)
library(circlize)
library(gplots)
library(viridis)
library(oompaBase)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

# 加载热图数据
hmdat <- read.csv("easy_input.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
# 或许你的数据行列与此相反，就这样转置一下
#hmdat <- t(hmdat)

# 加载分组
risk <- read.csv("easy_input_group.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
rownames(risk) <- paste0(rownames(risk),"-01")
dim(hmdat)
hmdat[1:3,1:3]
dim(risk)
head(risk)
risk<-risk[,-c(3,4)]
# 从easy_input_type.csv读取分类信息
type <- read.csv("easy_input_type.csv", row.names = 1)
head(type)
# 取出共有样本更新数据
comsam <- intersect(rownames(risk),rownames(hmdat))
hmdat <- hmdat[comsam,]
risk <- risk[comsam,,drop = F]
dim(hmdat)

# 拆分不同算法结果，获得类的名字
#immMethod <- sapply(strsplit(colnames(hmdat),"_",fixed = T),"[",2) #用easy_input.csv列名里的算法信息
immMethod <- type$Methods # 用easy_input_type.csv的算法那一列
# 用pheatmap画图
library(pheatmap)

# 定义颜色
methods.col <- brewer.pal(n = length(unique(immMethod)),name = "Paired")

# 创建注释
# 列注释，位于热图顶端
annCol <- data.frame(RiskScore = scale(risk$riskScore),
                     RiskType = risk$risk,
                     # 以上是risk score和risk type两种注释，可以按照这样的格式继续添加更多种类的注释信息，记得在下面的annColors里设置颜色
                     row.names = rownames(risk),
                     stringsAsFactors = F)

# 行注释，位于热图左侧
annRow <- data.frame(Methods = factor(immMethod,levels = unique(immMethod)),
                     row.names = colnames(hmdat),
                     stringsAsFactors = F)

# 为各注释信息设置颜色
annColors <- list(Methods = c("TIMER" = methods.col[1], #行注释的颜色
                              "CIBERSORT" = methods.col[2],
                              "CIBERSORT-ABS" = methods.col[3],
                              "QUANTISEQ" = methods.col[4],
                              "MCPCOUNTER" = methods.col[5],
                              "XCELL" = methods.col[6],
                              "EPIC" = methods.col[7]),
                  # 下面是列注释的颜色，可依此设置更多注释的颜色
                  "RiskScore" = greenred(64), 
                  "RiskType" = c("high" = "red","low" = "blue"))

# 数据标准化
indata <- t(hmdat)
indata <- indata[,colSums(indata) > 0] # 确保没有富集全为0的细胞
plotdata <- standarize.fun(indata,halfwidth = 2)

# 样本按risk score排序
samorder <- rownames(risk[order(risk$riskScore),])

# pheatmap绘图
pheatmap::pheatmap(mat = as.matrix(plotdata[,samorder]), # 标准化后的数值矩阵
                   border_color = NA, # 无边框色
                   color = bluered(64), # 热图颜色为红蓝
                   cluster_rows = F, # 行不聚类
                   cluster_cols = F, # 列不聚类
                   show_rownames = T, # 显示行名
                   show_colnames = F, # 不显示列名
                   annotation_col = annCol[samorder,,drop = F], # 列注释
                   annotation_row = annRow, # 行注释
                   annotation_colors = annColors, # 注释颜色
                   gaps_col = table(annCol$RiskType)[2], # 列分割
                   gaps_row = cumsum(table(annRow$Methods)), # 行分割
                   cellwidth = 0.8, # 元素宽度
                   cellheight = 10, # 元素高度
                   filename = "immune heatmap by pheatmap.pdf")