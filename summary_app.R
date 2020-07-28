
library(tidyverse)
dat = read.csv("Dropbox/Tanja&XinW/Newdata/AggObese1990_2017.csv")
dat = dat[dat$AGE!=80,]

tm1 = (dat$IYEAR - 1989)/(2018-1989)
  
y1 = scale(dat$PropObese)



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



##### curve for eigenfunctions ###

eigenmat = read.table("Google Drive File Stream/My Drive/Research/FDA_subgroup/result/eigenfunmat.txt")
colnames(eigenmat) = c("value","time")

pdf("Google Drive File Stream/My Drive/Research/FDA_subgroup/doc/figures/eigenfunctions.pdf",width = 7,height = 5)
ggplot(data = eigenmat, aes(x = time, y = value)) + geom_line() + theme_bw()
dev.off()



#### correlation ####

Bmtm = read.table("Google Drive File Stream/My Drive/Research/FDA_subgroup/result/Bmtm.txt")
theta = read.table("Google Drive File Stream/My Drive/Research/FDA_subgroup/result/theta.txt")

eigenfuntm = as.matrix(Bmtm) %*% theta[,1]
covmat = corrmat = matrix(0, 28,28)

for(i in 1:27)
{
  for(j in (i+1):28)
    covmat[i,j] = eigenfuntm[i] * eigenfuntm[j]
}

covmat = covmat + t(covmat)
diag(covmat) = eigenfuntm* eigenfuntm
lamj = 0.032756463159806126


