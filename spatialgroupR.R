gridm <- expand.grid(12:1,1:12)[,2:1]

f1 <- function(x)
{
  0.7*x[1] + x[2] - 13
}

f2 <- function(x)
{
  0.75*x[1] - x[2]
}

value1 <- apply(gridm,1,f1)
value2 <- apply(gridm, 1, f2)

group <- rep(0, 12*12)

group[value1<0 & value2 <=0& gridm[,1]<7] <- 1
group[ gridm[,2] >=7 & group==0] <- 2
group[group==0] <- 3

table(group)
library(ggplot2)

df1 <- data.frame(x = gridm[,1], y = gridm[,2], group = group)

pdf("../doc/figures/groupsp.pdf",width = 6,height = 5)
ggplot(data = df1) + geom_point(aes(x = x, y= y, shape = as.factor(group),color = as.factor(group)), 
                           size = 3) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.title = element_blank())
dev.off()




### four groups ####


f1 <- function(x)
{
  sqrt(2) * sin(4*pi*x) + 1
}


f2 <- function(x)
{
  2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)
}


f3 <- function(x)
{
  sqrt(2) * sin(4*pi*x) + 0.5
}

f4 <- function(x)
{
  2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2) + 0.5
}


ngrid = 14
n = ngrid * ngrid
gridm = matrix(0,ngrid*ngrid, 2)
gridm[,1] = rep(1:ngrid, each =ngrid)
gridm[,2] = rep(ngrid:1,ngrid)

group = rep(0,n)

group[(gridm[,1] <= 7) & (gridm[,2]>7)] = 1
group[(gridm[,1] <= 7) & (gridm[,2]<=7)] = 3
group[(gridm[,1] > 7) & (gridm[,2]>7)] = 2
group[(gridm[,1] >7) & (gridm[,2]<=7)] = 4

library(ggplot2)

df2 <- data.frame(x = gridm[,1], y = gridm[,2], group = group)

g1 <- ggplot(data = df2) + 
  geom_point(aes(x = x, y= y, shape = as.factor(group),color = as.factor(group)), 
                                size = 1) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.title = element_blank())


data4s30_v2 <- readRDS("../simdata/data4s30_v2.rds")

dati <- data4s30_v2[,,1]
colnames(dati)  <- c("group","ind","time","obs","mean")
dati <- as.data.frame(dati)
dati$group0 <- dati$group


dati <- dati %>% mutate(group = case_when(group0 ==2 ~ 3,
                                          group0 ==3 ~2,
                                          group0 ==1 ~1,
                                          group0==4 ~4)) %>%
  mutate(group = as.factor(group))


g2 <- ggplot(dati, aes(x = time, y =obs, group = ind, color = group)) + 
  geom_line(alpha = 0.5) + theme_bw()

library(gridExtra)

pdf("../doc/figures/gridexample4.pdf",width = 7,height = 3)
grid.arrange(g1, g2, ncol = 2)
dev.off()

pdf("../doc/figures/grid4.pdf",width = 6,height = 5)
g1
dev.off()

pdf("../doc/figures/example4.pdf",width = 6,height = 5)
g2
dev.off()








### 4 groups more complicated example 


ngrid = 14
n = ngrid * ngrid
gridm = matrix(0,ngrid*ngrid, 2)
gridm[,1] = rep(1:ngrid, each =ngrid)
gridm[,2] = rep(ngrid:1,ngrid)

f1 <- function(x)
{
  -x[1] + x[2] - 6
}

f2 <- function(x)
{
  -x[1]+x[2]
}

f3 <- function(x)
{
  -x[1]+x[2] + 6
}



group = rep(0,n)

value1 <- apply(gridm,1,f1)
value2 <- apply(gridm, 1, f2)
value3 <- apply(gridm, 1, f3)
group <- rep(0, 12*12)

group[value1>=0 ] <- 1
group[value1 <0 & value2 >=0] <- 2
group[value2 <0 & value3 >0] <- 3
group[value3 <=0] <- 4



df3 <- data.frame(x = gridm[,1], y = gridm[,2], group = group)
df3$group1 <- df3$group

df3$group1[group==2] <- 3
df3$group1[group==3] <- 2

ggplot(data = df3) + 
  geom_point(aes(x = x, y= y, shape = as.factor(group1),color = as.factor(group1)), 
             size = 3) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.title = element_blank())


write.table(df3,"../simdata/group4_g2.csv", row.names = FALSE)
