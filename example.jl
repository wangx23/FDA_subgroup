#### a example based on the given data set ####

using CSV
using DataFrames
data = CSV.File("data.csv")|> DataFrame

include("initial.jl")
include("GrFDA.jl")
include("refitFDA.jl")
include("BIC.jl")
include("neigh.jl")


indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 7, stop = 1))[2:6]

nobstotal = length(unique(indexy))
wt = ones(convert(Int,nobstotal*(nobstotal)/2))
group = unique(data[:,1:2])[:,1]
meanfun = data.meanfun

### initial 
lamv = collect(range(0,20,step=0.5))
betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)
K0 = 12
res1l = GrFDA(indexy,tm,y,knots,5,wt,betam0v5,lam = 0,
K0 = K0,maxiter = 1000)
estpc1 = minimum((1:5)[cumsum(reverse(res1l.lamj))/sum(res1l.lamj) .>=0.95])


lamvec = collect(range(0.1,1,step = 0.01))
nlam = length(lamvec)

# equal weights
BICvec1= zeros(nlam)
for l = 1:nlam
    res1l = GrFDA(indexy,tm,y,knots,estpc1,wt,betam0v5,lam = lamvec[l],
    K0 = K0,maxiter = 500)
    BICvec1[l] = BICemCn(res1l)
end

index1 = argmin(BICvec1) ### not bounday of the grid of tuning parameters
res1 = GrFDA(indexy,tm,y,knots,estpc1,wt,betam0v5,
lam = lamvec[index1],K0=K0,maxiter = 500)
group1 = getgroup(res1.deltam,nobstotal)
ng1 = size(unique(group1))[1]
ari1 = randindex(group,group1)[1]
    

### spatial weights
Cmat = cmatfun(12,12)
ordermat = zeros(Int64, 12*12 , 12*12)
for i = 1:nobstotal
    ordermat[:,i] = calorder(i,Cmat)
end

alp2 = collect(range(0.05,1,step = 0.05))
nalp = length(alp2)


BICvec2 = zeros(nlam,nalp)
for l1 = 1:nalp
    wt2 = exp.(alp2[l1] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
    for l = 1:nlam
        res1l = GrFDA(indexy,tm,y,knots,estpc1,wt2,betam0v5,lam = lamvec[l],
        K0=K0,maxiter = 500)
        BICvec2[l,l1] = BICemCn(res1l)
    end
end

index2 = argmin(BICvec2) ### not bounday of the grid of tuning parameters
wt2 = exp.(alp2[index2[2]] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
res2 = GrFDA(indexy,tm,y,knots,estpc1,wt2,betam0v5,
lam = lamvec[index2[1]],
K0=K0,maxiter = 500)
group2 = getgroup(res2.deltam,nobstotal)
ng2 = size(unique(group2))[1]
ari2 = randindex(group,group2)[1]
 

[ng1,ari1, ng2, ari2]