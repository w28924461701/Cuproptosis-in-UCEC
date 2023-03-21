p2<-ggplot(data=ann, aes(x=risk, y=TIDE,fill=risk)) +
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs( y = "TIDE scores")+
 theme_classic()+theme(axis.text = element_text(size = 12))+
  theme(legend.position="top")


ann2<-dcast(ann1,risk~Response)
ann2<-melt(ann2,id=c("risk"))
 

df_cumsum <- ddply(ann2, "variable",
                   transform, 
                   label_ypos=cumsum(value) - 0.5*value)

p1<-ggplot(data=df_cumsum, aes(x=risk, y=value,fill=variable)) +
  geom_bar(stat="identity")+
  theme_classic()+
 theme(axis.text = element_text(size = 12))+
  theme(legend.position="top")+

  labs(fill='variable')+
  
  labs(y = "response to immune checkpoint therapy")


plot_grid(p1,p2)