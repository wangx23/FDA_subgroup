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

pdf("../doc/figures/groupsp.pdf",width = 6,height = 5)
qplot(x = gridm[,1],y = gridm[,2], shape = as.factor(group)) +
  theme_bw() +   theme(axis.title.x=element_blank(),
                       axis.text.x=element_blank(),
                       axis.ticks.x=element_blank())+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),legend.title = element_blank())
dev.off()
