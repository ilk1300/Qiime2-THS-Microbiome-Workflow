ggplot(s, aes(x = category, y = s$count)) +
  geom_boxplot(width = 0.4, fill = "white") +
  geom_jitter(aes(color = category, shape = category), 
              width = 0.1, size = 1) +
  scale_color_manual(values = c("#00AFBB", "#E7B800")) + 
  labs(x = NULL)   # Remove x axis label

#anova#
anova_Result=aov(a$count ~category, data = a)
summary(anova_Result)

###option 2###
boxplot(NUMS ~ GRP, data = ddf, lwd = 2, ylab = 'NUMS')
stripchart(NUMS ~ GRP, vertical = TRUE, data = ddf, 
           method = "jitter", add = TRUE, pch = 20, col = 'blue')

s=read.csv("s.csv",row.names=1,check.names=FALSE)

#Taking the mean
aggregate(s$count, list(s$cat),mean) #6677.93 new mean difference.
#5968.38 (count)

taketwo=read.csv("merged_otu.csv",row.names=1, check.names=FALSE)
t=taketwo[sort(colnames(taketwo),)]
tNew=t[,1:74]
tNew2=t[,75:126]
meanS=colMeans(tNew)
meanNS=colMeans(tNew2)
summary(meanS)
summary(meanNS)
#difference in mean: 1.02594