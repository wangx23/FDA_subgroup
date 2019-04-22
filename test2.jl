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

res0 = GrInd(indexy, tm, y, knots, 2, wt, betam0, lam = 0.3)
group0 = getgroup(res0.deltam, 100)
BICind(res0,1)

res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.3,maxiter = 1000)
group1 = getgroup(res1.deltam, 100)
BICem(res1,1)

res2 = GrFDA2(indexy,tm,y,knots,2,wt,betam0,lam = 0.3,
maxiter2 = 10, maxiter = 100)
group2 = getgroup(res2.deltam, 100)
BIC2(res2,1)

res3 = GrFDAproxy(indexy, tm, y, knots, wt, betam0, lam = 0.28, maxiter =1000)
group3 = getgroup(res3.deltam, 100)
BICproxy(res3,1)

refit3 = refitFDA(indexy, tm, y, knots, group3,2)

randindex(group,group3) ### ARI, RI


lamvec = collect(range(0.1,0.5,step = 0.02))
nlam = length(lamvec)
BICvec0 = zeros(nlam)
BICvec1 = zeros(nlam,3)
BICvec2 = zeros(nlam,3)



for l = 1:nlam
    res0l = GrInd(indexy, tm, y, knots, wt, betam0, lam = lamvec[l])
    BICvec0[l] = BICind(res0l,1)
end

for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0,lam = lamvec[l],maxiter = 1000)
        BICvec1[l,P] = BICem(res1l,1)
    end
end

res11 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = lamvec[l],maxiter = 1000)
BICem(res11,1)

res12 = GrFDA(indexy,tm,y,knots,3,wt,betam0,lam = lamvec[l],maxiter = 1000)

res13 = GrFDA(indexy,tm,y,knots,4,wt,betam0,lam = lamvec[l],maxiter = 1000)

res14 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = lamvec[16],maxiter = 1000)

BICem(res13,1)

for l = 1:nlam
    for P = 1:3
        res2l = GrFDA2(indexy,tm,y,knots,P,wt,betam0,lam = lamvec[l],
        maxiter2 = 100, maxiter = 200)
        BICvec2[l,P] = BIC2(res2l,1)
    end
end



BICvec3 = zeros(nlam)
for l = 1:nlam
    res3l = GrFDAproxy(indexy, tm, y, knots, wt, betam0, lam = lamvec[l], maxiter =1000)
    BICvec3[l] = BICproxy(res3l,1)
end



BICvec32 = zeros(3)

for P = 1:3
    refit3 = refitFDA(indexy, tm, y, knots, group3,2)
