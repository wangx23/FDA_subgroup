##### output for resultnew_v2 #######
library(ggplot2)
library(gridExtra)
setwd("Research/FDA_subgroup/resultnew_v2/")

#### m = 10 ####
resultsim3m10 = read.csv("resultsim3m10.csv",header = FALSE)

round(colMeans(resultsim3m10[,1:4]),2)
round(apply(resultsim3m10[,1:4],2,sd),3)

round(colMeans(resultsim3m10[,5:8]),2)
round(apply(resultsim3m10[,5:8],2,sd),2)
colMeans(resultsim3m10[,5:8]==3)

mean(resultsim3m10[,17])
sd(resultsim3m10[,17])
mean(resultsim3m10[,17]==2)

rmse3m10 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m10[,c(13,15)])))

plot10 = qplot(x = model, y = RMSE, data = rmse3m10, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


#### m = 20 ####
resultsim3m20 = read.csv("resultsim3m20.csv",header = FALSE)

round(colMeans(resultsim3m20[,1:4]),2)
round(apply(resultsim3m20[,1:4],2,sd),3)

round(colMeans(resultsim3m20[,5:8]),2)
round(apply(resultsim3m20[,5:8],2,sd),2)
colMeans(resultsim3m20[,5:8]==3)


mean(resultsim3m20[,17])
sd(resultsim3m20[,17])
mean(resultsim3m20[,17]==2)

rmse3m20 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m20[,c(13,15)])))

plot20 = qplot(x = model, y = RMSE, data = rmse3m20, main = "m=20") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


#### m = 30 ####
resultsim3m30 = read.csv("resultsim3m30.csv",header = FALSE)
resultsim3m30 = resultsim3m30[-23,]

round(colMeans(resultsim3m30[,1:4]),2)
round(apply(resultsim3m30[,1:4],2,sd),3)

round(colMeans(resultsim3m30[,5:8]),2)
round(apply(resultsim3m30[,5:8],2,sd),2)
colMeans(resultsim3m30[,5:8]==3)


mean(resultsim3m30[,17])
sd(resultsim3m30[,17])
mean(resultsim3m30[,17]==2)

rmse3m30 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m30[,c(13,15)])))

plot30 = qplot(x = model, y = RMSE, data = rmse3m30, main = "m=30") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


plotdf = rbind(rmse3m10, rmse3m20, rmse3m30) %>%
  mutate(m = rep(c(10,20,30),each = nrow(rmse3m10))) %>%
  mutate(m = factor(m, c(10,20,30),labels = c("m=10","m=20","m=30")))

pdf("../doc/figures/rmsem10v2.pdf",width = 6,height = 3)
ggplot(data = plotdf, aes(x = model, y= RMSE)) + 
  geom_boxplot() +
  facet_wrap(~m, nrow = 1) +
  theme_bw()
dev.off()



##### spatial result ####

resultsim3s10 = read.csv("resultsim3s10.csv",header = FALSE)

round(colMeans(resultsim3s10[,1:2]),2)
round(apply(resultsim3s10[,1:2],2,sd),3)

round(colMeans(resultsim3s10[,3:4]),2)
round(apply(resultsim3s10[,3:4],2,sd),2)
colMeans(resultsim3s10[,3:4]==3)

colMeans(resultsim3s10[,9:10])
round(apply(resultsim3s10[,9:10],2,sd),3)

rmse3s10 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3s10[,c(7,8)])))

plots10 = qplot(x = model, y = RMSE, data = rmse3s10, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


#### m = 20 ####
resultsim3s20 = read.csv("resultsim3s20.csv",header = FALSE)

round(colMeans(resultsim3s20[,1:2]),3)
round(apply(resultsim3s20[,1:2],2,sd),3)

round(colMeans(resultsim3s20[,3:4]),2)
round(apply(resultsim3s20[,3:4],2,sd),2)
colMeans(resultsim3s20[,3:4]==3)

colMeans(resultsim3s20[,9:10])
round(apply(resultsim3s20[,9:10],2,sd),3)

rmse3s20 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3s20[,c(7,8)])))

plots20 = qplot(x = model, y = RMSE, data = rmse3s20, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


#### m = 30 ####
resultsim3s30 = read.csv("resultsim3s30.csv",header = FALSE)

round(colMeans(resultsim3s30[,1:2]),3)
round(apply(resultsim3s30[,1:2],2,sd),3)

round(colMeans(resultsim3s30[,3:4]),2)
round(apply(resultsim3s30[,3:4],2,sd),2)
colMeans(resultsim3s30[,3:4]==3)

colMeans(resultsim3s30[,9:10])
round(apply(resultsim3s30[,9:10],2,sd),3)

rmse3s30 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3s30[,c(7,8)])))

plots30 = qplot(x = model, y = RMSE, data = rmse3s30, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


plotdfs = rbind(rmse3s10, rmse3s20, rmse3s30) %>%
  mutate(m = rep(c(10,20,30),each = nrow(rmse3s10))) %>%
  mutate(m = factor(m, c(10,20,30),labels = c("m=10","m=20","m=30")))

pdf("../doc/figures/rmsespv2.pdf",width = 6,height = 3)
ggplot(data = plotdfs, aes(x = model, y= RMSE)) + 
  geom_boxplot() +
  facet_wrap(~m, nrow = 1) +
  theme_bw()
dev.off()
