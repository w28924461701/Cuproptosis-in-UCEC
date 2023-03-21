rm(list = ls())
library(tidyverse)
library(gghalves)
library(dplyr)

setwd('D:\\STAG2 博士后项目一\\148cuproptosis\\04.symbol\\gdc_download_20220621_041850.268719')

rt<-read.table('symbol.txt',sep = '\t',header = T,check.names = F)
rt<-filter(rt,id%in%c('AL512353.1','ACOXL-AS1'))

rownames(rt)<-rt$id
rt$id<-NULL
rt<-t(rt)  #转置操作
rt<-as.data.frame(rt)
rt%>%mutate(Species=c(rep('normal',19),rep('tumor',408)))->rt
 
# 添加配对分组
rt$id<- c(rep(1:19),rep(1:408))


# 统计均值标准差
summ_rt<- rt %>% 
  group_by(Species) %>% 
  summarise(
    mean = mean(`ACOXL-AS1`),
    sd = sd(`ACOXL-AS1`),
    n = n()
  ) %>% 
  mutate(se = sd/sqrt(n),
         Species = factor(Species, levels = c('normal', 'tumor')))
summ_rt$id <- c(1,1)

# 数据分割 
rt_normal <- rt %>% 
  filter(Species == "normal")
rt_tumor <- rt %>% 
  filter(Species == "tumor")

summ_rt_normal <- summ_rt %>% 
  filter(Species == "normal")
summ_rt_tumor <- summ_rt %>% 
  filter(Species == "tumor")
# 绘图
library(gghalves)
library(ggpubr)
library(ggsignif)
library(ggsci)
p1<-ggplot(rt , aes(x = Species, y = `ACOXL-AS1`, fill = Species))+
  #配对连线
  geom_line(aes(group=id), color="gray" ,position = position_dodge(0.2)) +
  # 散点图
  geom_point(aes(color=Species, group=id), 
             size = 2, 
             position = position_dodge(0.2))+
  # 两侧的小提琴图，position设置对称位置
  geom_half_violin(aes(fill = Species),data = rt_normal, side = 'r',
                   position = position_nudge(x = .20, y = 0))+
  geom_half_violin(aes(fill = Species),data = rt_tumor, side = 'l',
                   position = position_nudge(x = -0.20, y = 0))+
  # 两侧的箱线图，position设置对称位置
  geom_boxplot(data =  rt_normal,
               aes(x = Species,y = `ACOXL-AS1`, fill = Species),
               outlier.shape = NA,
               width = .05,
               color = "black",
               position = position_nudge(x = 0.15, y = 0))+
  geom_boxplot(data = rt_tumor,
               aes(x = Species,y = `ACOXL-AS1`, fill = Species),
               outlier.shape = NA,
               width = .05,
               color = "black",
               position = position_nudge(x = -0.15, y = 0))+
  #两侧的误差棒图图，position设置对称位置
  geom_errorbar(data = summ_rt_normal,
                aes(x = Species, y = mean, group = Species, colour = Species,
                    ymin = mean-sd, ymax = mean+sd),
                width=0.1,size=3,
                position=position_nudge(x = -0.2, y = 0)
  ) +
  # 均值的点 大圆套小圆即两个点大小不一样 图层叠加
  geom_point(data=summ_rt_normal,
             aes(x=Species,y = mean,group = Species, color = Species),
             size = 10,
             position = position_nudge(x = -0.2,y = 0)) +
  geom_point(data=summ_rt_normal,
             aes(x=Species,y = mean,group = Species),
             color = "black",
             size = 6,
             position = position_nudge(x = -0.2,y = 0))+
  geom_errorbar(data = summ_rt_tumor,
                aes(x = Species, y = mean, group = Species, colour = Species,
                    ymin = mean-sd, ymax = mean+sd),
                width=0.1,size=3,
                position=position_nudge(x = 0.2, y = 0)
  ) +
  geom_point(data=summ_rt_tumor,
             aes(x=Species,y = mean,group = Species, color = Species),
             size = 10,
             position = position_nudge(x = 0.2,y = 0)) +
  geom_point(data=summ_rt_tumor,
             aes(x=Species,y = mean,group = Species),
             color = "black",
             size = 6,
             position = position_nudge(x = 0.2,y = 0))+
  # 连线 暴力修改连线长度
  geom_line(data =summ_rt, aes(x=-as.numeric(Species)*0.6,
                                 y=mean,group=id), color="black",
            position = position_nudge(x = 2.4,y = 0),
            size = 2)+
  scale_color_jco() +
  scale_fill_jco() +
  # 差异检验
  geom_signif(comparisons = list(c("normal", "tumor"))) +
  theme_classic()
ylim1<-boxplot.stats(rt$`ACOXL-AS1`)$stats[c(1, 3)] 
p1<-p1+coord_cartesian(ylim = ylim1*4) 






# 统计均值标准差
summ_rt<- rt %>% 
  group_by(Species) %>% 
  summarise(
    mean = mean(AL512353.1),
    sd = sd(`AL512353.1`),
    n = n()
  ) %>% 
  mutate(se = sd/sqrt(n),
         Species = factor(Species, levels = c('normal', 'tumor')))
summ_rt$id <- c(1,1)

# 数据分割 
rt_normal <- rt %>% 
  filter(Species == "normal")
rt_tumor <- rt %>% 
  filter(Species == "tumor")

summ_rt_normal <- summ_rt %>% 
  filter(Species == "normal")
summ_rt_tumor <- summ_rt %>% 
  filter(Species == "tumor")
p2<-ggplot(rt , aes(x = Species, y = AL512353.1, fill = Species))+
  #配对连线
  geom_line(aes(group=id), color="gray" ,position = position_dodge(0.2)) +
  # 散点图
  geom_point(aes(color=Species, group=id), 
             size = 2, 
             position = position_dodge(0.2))+
  # 两侧的小提琴图，position设置对称位置
  geom_half_violin(aes(fill = Species),data = rt_normal, side = 'r',
                   position = position_nudge(x = .20, y = 0))+
  geom_half_violin(aes(fill = Species),data = rt_tumor, side = 'l',
                   position = position_nudge(x = -0.20, y = 0))+
  # 两侧的箱线图，position设置对称位置
  geom_boxplot(data =  rt_normal,
               aes(x = Species,y = AL512353.1, fill = Species),
               outlier.shape = NA,
               width = .05,
               color = "black",
               position = position_nudge(x = 0.15, y = 0))+
  geom_boxplot(data = rt_tumor,
               aes(x = Species,y = AL512353.1, fill = Species),
               outlier.shape = NA,
               width = .05,
               color = "black",
               position = position_nudge(x = -0.15, y = 0))+
  #两侧的误差棒图图，position设置对称位置
  geom_errorbar(data = summ_rt_normal,
                aes(x = Species, y = mean, group = Species, colour = Species,
                    ymin = mean-sd, ymax = mean+sd),
                width=0.1,size=3,
                position=position_nudge(x = -0.2, y = 0)
  ) +
  # 均值的点 大圆套小圆即两个点大小不一样 图层叠加
  geom_point(data=summ_rt_normal,
             aes(x=Species,y = mean,group = Species, color = Species),
             size = 10,
             position = position_nudge(x = -0.2,y = 0)) +
  geom_point(data=summ_rt_normal,
             aes(x=Species,y = mean,group = Species),
             color = "black",
             size = 6,
             position = position_nudge(x = -0.2,y = 0))+
  geom_errorbar(data = summ_rt_tumor,
                aes(x = Species, y = mean, group = Species, colour = Species,
                    ymin = mean-sd, ymax = mean+sd),
                width=0.1,size=3,
                position=position_nudge(x = 0.2, y = 0)
  ) +
  geom_point(data=summ_rt_tumor,
             aes(x=Species,y = mean,group = Species, color = Species),
             size = 10,
             position = position_nudge(x = 0.2,y = 0)) +
  geom_point(data=summ_rt_tumor,
             aes(x=Species,y = mean,group = Species),
             color = "black",
             size = 6,
             position = position_nudge(x = 0.2,y = 0))+
  # 连线 暴力修改连线长度
  geom_line(data =summ_rt, aes(x=-as.numeric(Species)*0.6,
                               y=mean,group=id), color="black",
            position = position_nudge(x = 2.4,y = 0),
            size = 2)+
  scale_color_jco() +
  scale_fill_jco() +
  # 差异检验
  geom_signif(comparisons = list(c("normal", "tumor"))) +
  theme_classic()
ylim1<-boxplot.stats(rt$AL512353.1)$stats[c(1, 3)] 
p2<-p2+coord_cartesian(ylim = ylim1*4) 
library(cowplot)
plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

ggsave(filename = "test.pdf", device="pdf", height=6, width=9, useDingbats=FALSE)