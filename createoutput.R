#### summary for EM FDA #####
library(ggplot2)
library(gridExtra)
setwd("Research/FDA_subgroup/resultnew/")

resultsim3m10 = read.csv("resultsim3m10.csv",header = FALSE)

round(colMeans(resultsim3m10[,1:4]),2)
round(apply(resultsim3m10[,1:4],2,sd),3)

round(colMeans(resultsim3m10[,5:8]),2)
round(apply(resultsim3m10[,5:8],2,sd),2)
colMeans(resultsim3m10[,5:8]==3)

mean(resultsim3m10[,13])
sd(resultsim3m10[,13])
mean(resultsim3m10[,13]==2)

rmse3m10 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m10[,c(9,11)])/sqrt(150)))

plot10 = qplot(x = model, y = RMSE, data = rmse3m10, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


resultsim3m20 = read.csv("resultsim3m20.csv",header = FALSE)
round(colMeans(resultsim3m20[,1:4]),2)
round(apply(resultsim3m20[,1:4],2,sd),3)

round(colMeans(resultsim3m20[,5:8]),2)
round(apply(resultsim3m20[,5:8],2,sd),2)
colMeans(resultsim3m20[,5:8]==3)

mean(resultsim3m20[,13])
sd(resultsim3m20[,13])


rmse3m20 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m20[,c(9,11)])/sqrt(150)))

plot20 = qplot(x = model, y = RMSE, data = rmse3m20, main = "m = 20") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))



resultsim3m30 = read.csv("resultsim3m30.csv",header = FALSE)
round(colMeans(resultsim3m30[,1:4]),2)
round(apply(resultsim3m30[,1:4],2,sd),3)

round(colMeans(resultsim3m30[,5:8]),2)
round(apply(resultsim3m30[,5:8],2,sd),2)
colMeans(resultsim3m30[,5:8]==3)

mean(resultsim3m30[,13])
sd(resultsim3m30[,13])


rmse3m30 = data.frame(model = c(rep("IND",100), rep("FDA",100)), 
                      RMSE = c(as.matrix(resultsim3m30[,c(9,11)])/sqrt(150)))

plot30 = qplot(x = model, y = RMSE, data = rmse3m30, main = "m = 30") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))

pdf("../doc/figures/rmsem10.pdf",width = 6,height = 3)
grid.arrange(plot10, plot20, plot30, ncol = 3)
dev.off()


resultsim3m10ncl = read.csv("resultsim3m10ncl100.csv",header = FALSE)

round(colMeans(resultsim3m10ncl[,1:3]),2)
round(apply(resultsim3m10ncl[,1:3],2,sd),3)


round(colMeans(resultsim3m10ncl[,4:6]),2)
round(apply(resultsim3m10ncl[,4:6],2,sd),2)
colMeans(resultsim3m10ncl[,4:6]==3)

mean(resultsim3m10ncl[,10])
sd(resultsim3m10ncl[,10])


resultsim3m20ncl = read.csv("resultsim3m20ncl100.csv",header = FALSE)

round(colMeans(resultsim3m20ncl[,1:3]),2)
round(apply(resultsim3m20ncl[,1:3],2,sd),3)


round(colMeans(resultsim3m20ncl[,4:6]),2)
round(apply(resultsim3m20ncl[,4:6],2,sd),2)
colMeans(resultsim3m20ncl[,4:6]==3)

mean(resultsim3m20ncl[,10])
sd(resultsim3m20ncl[,10])

##### spatial result ####

resultsim3s10 = read.csv("resultsim3s10.csv",header = FALSE)

round(colMeans(resultsim3s10[,1:2]),3)
round(apply(resultsim3s10[,1:2],2,sd),3)

round(colMeans(resultsim3s10[,3:4]),2)
round(apply(resultsim3s10[,3:4],2,sd),2)

colMeans(resultsim3s10[,7:8])
round(apply(resultsim3s10[,7:8],2,sd),3)


rmse3s10 = data.frame(model = c(rep("equal",100), rep("sp",100)), 
                      RMSE = c(as.matrix(resultsim3s10[,c(5,6)])/sqrt(48*3)))

plots10 = qplot(x = model, y = RMSE, data = rmse3s10, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))


resultsim3s20 = read.csv("resultsim3s20.csv",header = FALSE)

round(colMeans(resultsim3s20[,1:2]),3)
round(apply(resultsim3s20[,1:2],2,sd),3)

round(colMeans(resultsim3s20[,3:4]),2)
round(apply(resultsim3s20[,3:4],2,sd),2)

colMeans(resultsim3s20[,7:8])
apply(resultsim3s20[,7:8],2,sd)


rmse3s20 = data.frame(model = c(rep("equal",100), rep("sp",100)), 
                      RMSE = c(as.matrix(resultsim3s20[,c(5,6)])/sqrt(48*3)))

plots20 = qplot(x = model, y = RMSE, data = rmse3s20, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))



resultsim3s30 = read.csv("resultsim3s30.csv",header = FALSE)

round(colMeans(resultsim3s30[,1:2]),3)
round(apply(resultsim3s30[,1:2],2,sd),3)

round(colMeans(resultsim3s30[,3:4]),2)
round(apply(resultsim3s30[,3:4],2,sd),2)

colMeans(resultsim3s30[,7:8])
apply(resultsim3s30[,7:8],2,sd)


rmse3s30 = data.frame(model = c(rep("equal",100), rep("sp",100)), 
                      RMSE = c(as.matrix(resultsim3s30[,c(5,6)])/sqrt(48*3)))

plots30 = qplot(x = model, y = RMSE, data = rmse3s30, main = "m=10") + 
  geom_boxplot() + theme_bw() +  theme(plot.title = element_text(hjust = 0.5))

pdf("../doc/figures/rmsesp.pdf",width = 6,height = 3)
grid.arrange(plots10, plots20, plots30, ncol = 3)
dev.off()





