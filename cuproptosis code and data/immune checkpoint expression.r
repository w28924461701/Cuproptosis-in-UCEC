library(tidyverse)
library(gapminder)
library(ggsci)
library(ggprism)
library(rstatix)
library(ggpubr)


data=cbind(dataLow,dataHigh)
data=data[rowMeans(data)>1,]
conNum=ncol(dataLow)
treatNum=ncol(dataHigh)
Type=c(rep('low-risk',conNum), rep('high-risk',treatNum))
View(data)
data<-t(data)
data<-as.data.frame(data)
data$Type<-Type

library(reshape2)
short2long = melt(data2, id=c("Type"),
                   variable.name= 'genes', value.name = 'value')

short2long$Type<-as.factor(short2long$Type)

ggplot(short2long,aes(x=genes,y=value,fill=Type))+
  geom_boxplot(width=0.9)+
  geom_jitter(position=position_jitter(0.05))+
  theme_classic()+
  theme(legend.position="top",axis.text.x = element_text(size=10))+
  ylim(0,15)+
  labs(title="immune checkpoint",y = "immune checkpoint expression value")


df_p_val1 <- short2long %>% group_by(genes) %>%
  wilcox_test(value  ~ Type) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>% 
  add_xy_position(x = "genes", dodge = 0.8) 


p1<-ggplot(short2long,aes(x=genes,y=value,fill=Type))+
  stat_boxplot(geom="errorbar",position=position_dodge(width=0.2),width=0.8)+
  geom_boxplot(position=position_dodge(width =1.0),width=0.6)+
  geom_line(aes(group=genes),position = position_dodge(0.2),color="red") +
  geom_point(aes(fill=genes,group=Type,alpha=value),pch=21,
             position = position_dodge(0.2))+
  theme_classic()+
  ylim(0,8)+theme(legend.position="top",axis.text.x = element_text(size=10))+
  labs(title="immune checkpoint",y = "immune checkpoint expression gene value")