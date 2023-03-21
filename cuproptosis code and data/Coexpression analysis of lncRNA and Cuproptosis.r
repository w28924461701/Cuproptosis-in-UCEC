rt1 <- data.frame(rt)
head(rt)
gg1<-ggplot(data = rt1,aes(axis1 = lncRNA,axis2 = Cuproptosis,
           weight = cor)) +
  scale_x_discrete(limits = c("lncRNA", "Cuproptosis"), expand = c(.1, .05)) +
  geom_alluvium(aes(fill = Cuproptosis)) +
  geom_stratum(alpha = .5) +
  theme_minimal() +
  ggtitle("Coexpression analysis of lncRNA and Cuproptosis")


pdf(file="gg1.pdf", width=6, height=4.75)
print(gg1)
dev.off()