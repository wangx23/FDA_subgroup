### test GrFDA and GrFDA2 ####

include("initial.jl")
include("simdat4.jl")
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

data4 = simdat4(sig2, lamj, m = m, ncl = ncl,seed = 10)

indexy4 = data4.ind
tm4 = data4.time
y4 = data4.obs
knots4 = collect(range(0,length = 9, stop = 1))[2:7]
boundary = [0,1]
maxiter = 1000
P = 2
nu = 1
gam = 3
tolabs = 1e-4
tolrel = 1e-2
wt = ones(convert(Int,200*199/2))
group = unique(data4[:,1:2])[:,1]

lamv = collect(range(0,10,step=0.5))

betam04v5 = initial5(indexy4,tm4, y4, knots4, lamv = lamv)
betam04v3 = initial3(indexy4,tm4, y4, knots4, lamv = lamv)


lamvec = collect(range(0.15,0.4,step = 0.01))
nlam = length(lamvec)


uniqtm4 = unique(tm4)
Bmi4 = orthogonalBsplines(uniqtm4,knots4)

tr(Bmi4* inv(transpose(Bmi4) * Bmi4) * transpose(Bmi4))



### Ind
BICvec0 = zeros(nlam)
for l = 1:nlam
    res0l = GrInd(indexy4, tm4, y4, knots4, wt, betam04v5, lam = lamvec[l])
    BICvec0[l] = BICind(res0l,1)
end

argmin(BICvec0)
res0 = GrInd(indexy4, tm4, y4, knots4, wt, betam04v5, lam = lamvec[1])
group0 = getgroup(res0.deltam,200)
randindex(group,group0)


BICvec01 = zeros(nlam)
for l = 1:nlam
    res0l = GrInd(indexy4, tm4, y4, knots4, wt, betam04v5, lam = lamvec[l])
    BICvec01[l] = BICind0(res0l)
end

index01 = argmin(BICvec01)
res01 = GrInd(indexy4, tm4, y4, knots4, wt, betam04v5, lam = lamvec[index01])
group01 = getgroup(res01.deltam,200)
randindex(group,group01)


## EM

BICvec1 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy4,tm4,y4,knots4,P,wt,betam04v5,
        lam = lamvec[l],
        K0=10,maxiter = 1000)
        BICvec1[l,P] = BICem(res1l)
    end
end

argmin(BICvec1)

res1 = GrFDA(indexy4,tm4,y4,knots4,2,wt,betam04v5,lam = lamvec[11],
K0=10,maxiter = 1000)

group1 = getgroup(res1.deltam,200)
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
