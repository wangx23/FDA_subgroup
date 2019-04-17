#### test  #####

## refitFDA ##

include("simdat.jl")
include("refitFDA.jl")
include("initial.jl")

m = 50
ncl = 50
sig2 = 0.1
lamj = [0.1,0.2]

data = simdat(sig2, lamj, m = m, ncl = ncl)


indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 5, stop = 1))[2:4]
boundary = [0,1]
maxiter = 1000
P = 2

tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]
mean1 = m1.(tvec)
mean2 = m2.(tvec)

group = unique(data[:,1:2])[:,1]


using LinearAlgebra
using Clustering


resg1 = refitFDA(indexy, tm, y, knots, group, 2,maxiter = 10)



resg0 = refitFDA(indexy, tm, y, knots, group, 2, maxiter = 1000)


fixalpm = resg0.alpm
fixlamj = resg0.lamj
fixsig2 = resg0.sig2
fixtheta = resg0.theta
