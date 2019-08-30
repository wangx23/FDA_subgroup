#### application ###
using CSV
using Statistics
using Random
include("initial.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrInd.jl")
include("GrFDAproxy.jl")
include("refitFDA.jl")
include("BIC.jl")




function standfun(x::Number, lower::Number, upper::Number)
xnew = (x - lower)/(upper - lower)
end

datp = CSV.read("../doc/dat656.csv")
indexy = datp.Station
tm = standfun.(datp.Year, 1960, 2015)
y = (datp.max88 .- mean(datp.max88))./std(datp.max88)
knots = collect(range(0,length = 5, stop = 1))[2:4]
betam0 = initial2(indexy, tm, y, knots, lam = 10)

wt = ones(convert(Int,656*655/2))

res0 = GrInd(indexy,tm,y,knots,wt,betam0,lam = 0.3,maxiter = 1000)
res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.2,maxiter = 1000)

groupest = getgroup(res1.deltam, 656)


using FreqTables

datp1 = CSV.read("../doc/AggObese1990_2017.csv")
datp1 = datp1[datp1.AGE .!=80,:]

lamv = collect(range(0,20,step=0.5))
indexy1 = datp1.AGE
tm1 = standfun.(datp1.IYEAR, 1989, 2018)
y1 = (datp1.PropObese .- mean(datp1.PropObese))./std(datp1.PropObese)
nage = size(unique(datp1.AGE))[1]



knots1 = collect(range(0,length = 6, stop = 1))[2:5]
betam01 = initial5(indexy1, tm1, y1, knots1, lamv = lamv)

wt1 = ones(convert(Int,62*61/2))

res10 = GrInd(indexy1,tm1,y1,knots1,wt1,betam01, lam  = 0.2, maxiter = 1000)
group10 = getgroup(res10.deltam,62)
freqtable(group10)

res11 = GrFDA(indexy1,tm1,y1,knots1,2,wt1,betam01,
lam = 0.2,maxiter = 1000)
group11 = getgroup(res11.deltam,62)
freqtable(group11)

freqtable(group10, group11)

lamvec = collect(range(0.05,0.45,step = 0.01))
nlam = length(lamvec)

BICvec0 = zeros(nlam)
BICvec01 = zeros(nlam)
for l = 1:nlam
    res0l = GrInd(indexy1,tm1,y1,knots1,wt1,betam01, lam  = lamvec[l])
    BICvec0[l] = BICind2(res0l,1)
    BICvec01[l] = BICind0(res0l)
end

argmin(BICvec0)
res10 = GrInd(indexy1,tm1,y1,knots1,wt1,betam01, lam  = lamvec[4], maxiter = 1000)
group10 = getgroup(res10.deltam,62)
freqtable(group10)


BICvec11 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res11l = GrFDA(indexy1,tm1,y1,knots1,P,wt1,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec11[l,P] = BICem2(res11l)
    end
end

argmin(BICvec11)


res11 = GrFDA(indexy1,tm1,y1,knots1,2,wt1,betam01,
lam = lamvec[32],maxiter = 1000)

res12 = GrFDA(indexy1,tm1,y1,knots1,2,wt1,betam01,
lam = lamvec[10],maxiter = 1000)

res13 = GrFDA(indexy1,tm1,y1,knots1,1,wt1,betam01,
lam = lamvec[10],maxiter = 1000)

res14 = GrFDA(indexy1,tm1,y1,knots1,3,wt1,betam01,
lam = lamvec[30],maxiter = 1000)


group11 = getgroup(res11.deltam,62)


#### considering weights? #######

wmat =  zeros(nage, nage)
for i in 1:(nage-1)
  for j in (i+1):(nage)
      wmat[i,j] = abs(i - j)
    end
end


ind1 = 1
ordervec = zeros(convert(Int, nage*(nage-1)/2))
for i in 1:(nage - 1)
    ind2 = convert(Int,(2*nage - i - 1)*i/2)
    ordervec[ind1:ind2] = wmat[i,(i+1):nage]
    global ind1 = ind2+1
end


wt2 =  exp.(0.5.*(1 .- ordervec))

BICvec12 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res12l = GrFDA(indexy1,tm1,y1,knots1,P,wt2,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec12[l,P] = BICem2(res12l)
    end
end

wt2 =  exp.(0.5*(1 .- ordervec))
res13 = GrFDA(indexy1,tm1,y1,knots1,3,wt2,betam01,lam = lamvec[33],
K0 = 10)

freqtable(res13.groupest)


### select part of the 656 stations #####
