##### a function to simulate and return results ###
include("initial.jl")
include("scad.jl")
include("GrInd.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrFDAproxy.jl")
include("refitFDA.jl")
include("BIC.jl")
include("simdat3.jl")
include("complam.jl")
using Dates
using Distributed
using RCall


m = 30
ncl = 50
sig2 = 0.04
lamj = [0.1,0.2]

data = simdat3(sig2, lamj, m = m, ncl = ncl, seed = 1)

indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 7, stop = 1))[2:6]

nobstotal = length(unique(indexy))
wt = ones(convert(Int,nobstotal*(nobstotal)/2))
group = unique(data[:,1:2])[:,1]
meanfun = data.meanfun


#betam0 = initial2(indexy, tm, y, knots, lam = 5)
lamv = collect(range(0,20,step=0.5))
betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)

resor = refitFDA(indexy,tm,y,knots,group,2, betam0v5)
betaor = transpose(resor.alpm[:,group])

lamvec = collect(range(0.2,0.5,step = 0.01))
nlam = length(lamvec)


## EM
BICvec1 = zeros(nlam,3)
ngmat1 = zeros(nlam,3)
arimat1 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0v5,lam = lamvec[l],
        K0 = 12,maxiter = 1000)
        BICvec1[l,P] = BICem4(res1l,1)
        group1 = getgroup(res1l.deltam,nobstotal)
        ngmat1[l,P] = size(unique(group1))[1]
        arimat1[l,P] = randindex(group,group1)[1]
    end
end

index1 = argmin(BICvec1)

res1 = GrFDA(indexy,tm,y,knots,index1[2],wt,betam0v5,
lam = lamvec[index1[1]],
K0=12,maxiter = 1000)
group1 = getgroup(res1.deltam,nobstotal)
ng1 = size(unique(group1))[1]
ari1 = randindex(group,group1)[1]
norm1 = norm(res1.betaest - betaor)/sqrt(150)
estpc1 = index1[2]
rmse1 = norm(res1.meanfunest - meanfun)/sqrt(150*30)

### accuray of lamj
rmselam = complam(res1.lamj,lamj,index1[2])

@rlibrary fclust2

### something needs to redone #####
R" datlist = list(x = $y,
                 curve = $indexy,
                 timeindex = match(dat$time,grids))"

R"fitfclust(data=datlist,q=7,h=Kjs[j]-1,p=8,K=Kjs[j],
    maxit=30,grid=grids,plot=F,trace=F)"
