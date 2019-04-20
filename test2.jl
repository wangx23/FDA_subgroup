### test GrFDA and GrFDA2 ####

include("initial.jl")
include("simdat.jl")
include("scad.jl")
include("GrFDA.jl")
include("GrFDA2.jl")



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
nu = 1
gam = 3
tolabs = 1e-4
tolrel = 1e-2
wt = ones(convert(Int,100*99/2))

tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]
mean1 = m1.(tvec)
mean2 = m2.(tvec)

group = unique(data[:,1:2])[:,1]

betam0 = initial2(indexy, tm, y, knots, lam = 10)


res2 = GrFDA2(indexy,tm,y,knots,2, wt, betam0, lam = 0.3)
