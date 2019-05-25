### test GrFDA and GrFDA2 ####

include("initial.jl")
include("simdat3.jl")
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

data3 = simdat3(sig2, lamj, m = m, ncl = ncl,seed = 10)

indexy3 = data3.ind
tm3 = data3.time
y3 = data3.obs
knots3 = collect(range(0,length = 7, stop = 1))[2:6]
boundary = [0,1]
maxiter = 1000
P = 2
nu = 1
gam = 3
tolabs = 1e-4
tolrel = 1e-2
wt = ones(convert(Int,150*149/2))
group = unique(data3[:,1:2])[:,1]

lamv = collect(range(0,10,step=0.5))

betam03v5 = initial5(indexy3,tm3, y3, knots3, lamv = lamv)


lamvec = collect(range(0.15,0.4,step = 0.01))
nlam = length(lamvec)


uniqtm3 = unique(tm)
Bmi3 = orthogonalBsplines(uniqtm3,knots)

tr(Bmi3* inv(transpose(Bmi3) * Bmi3) * transpose(Bmi3))



### Ind
BICvec0 = zeros(nlam)
for l = 1:nlam
    res0l = GrInd(indexy3, tm3, y3, knots3, wt, betam03v5, lam = lamvec[l])
    BICvec0[l] = BICind(res0l,1)
end

argmin(BICvec0)
res0 = GrInd(indexy3, tm3, y3, knots3, wt, betam03v5, lam = lamvec[1])
group0 = getgroup(res0.deltam,150)
randindex(group,group0)


## EM

BICvec1 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy3,tm3,y3,knots3,P,wt,betam03v5,lam = lamvec[l],
        K0=12,maxiter = 1000)
        BICvec1[l,P] = BICem(res1l)
    end
end

argmin(BICvec1)

res1 = GrFDA(indexy3,tm3,y3,knots3,2,wt,betam03v5,lam = lamvec[14],
K0=12,maxiter = 1000)

group1 = getgroup(res1.deltam,150)
randindex(group,group1)



### iterative steps
BICvec2 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res2l = GrFDA2(indexy3,tm3,y3,knots3,P,wt,betam03v5,lam = lamvec[l],
        K0 = 12, maxiter2 = 100, maxiter = 50)
        BICvec2[l,P] = BIC2(res2l,1)
    end
end

argmin(BICvec2)
res2 = GrFDA2(indexy3,tm3,y3,knots3,2,wt,betam03v5,lam = lamvec[16],
K0=12, maxiter2 = 100, maxiter = 50)
group2 = getgroup(res2.deltam,150)
randindex(group,group2)


### proxy
BICvec3 = zeros(nlam)
for l = 1:nlam
    res3l = GrFDAproxy(indexy3, tm3, y3, knots3, wt, betam03v5,
    lam = lamvec[l], maxiter =1000)
    BICvec3[l] = BICproxy(res3l,1)
end

argmin(BICvec3)

res3 = GrFDAproxy(indexy3, tm3, y3, knots3, wt, betam03v5,
lam = lamvec[22], maxiter =1000)

group3 = getgroup(res3.deltam,150)
randindex(group,group3)

BICvec32 = zeros(3)

for P = 1:3
    refit3p = refitFDA(indexy, tm, y, knots, group3, P)
    BICvec32[P] = BICrefit(refit3p)
end

refit3 = refitFDA(indexy, tm, y, knots, group3, 2)
