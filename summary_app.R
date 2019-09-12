
library(tidyverse)
dat = read.csv("Dropbox/Tanja&XinW/Newdata/AggObese1990_2017.csv")
dat = dat[dat$AGE!=80,]


pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup/doc/figures/datapattern.pdf",width = 5,height = 4)
ggplot(data = dat, aes(x = IYEAR, y = PropObese, group = AGE, color = AGE)) +
  geom_line() + 
  labs(x = "Year",y = "Proportion of Obesity") +
  theme_bw()
dev.off()

groupest = read.table("Google Drive File Stream/My Drive/Research/FDA_subgroup/result/groupest.txt")




meanfunest = read.table("Google Drive File Stream/My Drive/Research/FDA_subgroup/result/meanfunest52.txt")

datpred = dat
datpred$pred = meanfunest$V1
datpred = arrange(datpred, AGE, IYEAR) %>%
  mutate(groupest = rep(groupest$V1,each = length(unique(datpred$IYEAR))))%>%
  mutate(groupest = factor(groupest, 1:7, labels = c("18-20","21-23","24-30","31-38","39-50","51-73","74-79")))


cbind(unique(dat$AGE), groupest$V1)


##### figure #####
pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup/doc/figures/groupfit.pdf",width = 7,height = 5)
ggplot(data = datpred) + 
  geom_line(aes(x = IYEAR, y = PropObese, group = AGE)) +
  geom_line(aes(x = IYEAR, y = pred, color = "red"), size = 1) +
  facet_wrap(~groupest, ncol = 4) +
  labs( x= "Year", y = "Proportion of Obesity") + 
  theme_bw() +
  theme(legend.position = "none")
dev.off()




