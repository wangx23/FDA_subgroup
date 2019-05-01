### test GrFDA and GrFDA2 ####

include("initial.jl")
include("simdat.jl")
include("scad.jl")
include("GrInd.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrFDAproxy.jl")
include("refitFDA.jl")
include("BIC.jl")



m = 50
ncl = 50
sig2 = 0.1
lamj = [0.1,0.2]

data = simdat(sig2, lamj, m = m, ncl = ncl,seed = 2)

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
group = unique(data[:,1:2])[:,1]

betam0 = initial2(indexy, tm, y, knots, lam = 10)



lamvec = collect(range(0.1,0.5,step = 0.02))
nlam = length(lamvec)


### Ind
BICvec0 = zeros(nlam)
for l = 1:nlam
    res0l = GrInd(indexy, tm, y, knots, wt, betam0, lam = lamvec[l])
    BICvec0[l] = BICind(res0l,1)
end

argmin(BICvec0)
res0 = GrInd(indexy, tm, y, knots, wt, betam0, lam = lamvec[1])
group0 = getgroup(res0.deltam,100)
randindex(group,group0)


## EM
BICvec1 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0,lam = lamvec[l],maxiter = 1000)
        BICvec1[l,P] = BICem(res1l,1)
    end
end

argmin(BICvec1)[1]

res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = lamvec[12],maxiter = 1000)
group1 = getgroup(res1.deltam,100)
randindex(group,group1)



### iterative steps
BICvec2 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res2l = GrFDA2(indexy,tm,y,knots,P,wt,betam0,lam = lamvec[l],
        maxiter2 = 100, maxiter = 50)
        BICvec2[l,P] = BIC2(res2l,1)
    end
end
argmin(BICvec2)
res2 = GrFDA2(indexy,tm,y,knots,2,wt,betam0,lam = lamvec[6], maxiter2 = 100, maxiter = 50)
group2 = getgroup(res2.deltam,100)
randindex(group,group2)


### proxy
BICvec3 = zeros(nlam)
for l = 1:nlam
    res3l = GrFDAproxy(indexy, tm, y, knots, wt, betam0, lam = lamvec[l], maxiter =1000)
    BICvec3[l] = BICproxy(res3l,1)
end

argmin(BICvec3)

res3 = GrFDAproxy(indexy, tm, y, knots, wt, betam0, lam = lamvec[10], maxiter =1000)

group3 = getgroup(res3.deltam,100)
randindex(group,group3)

BICvec32 = zeros(3)

for P = 1:3
    refit3p = refitFDA(indexy, tm, y, knots, group3, P)
    BICvec32[P] = BICrefit(refit3p)
end

refit3 = refitFDA(indexy, tm, y, knots, group3, 2)
