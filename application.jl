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

lamv = collect(range(0,20,step=0.5))

datp = CSV.read("../doc/dat656.csv")
datp = datp[datp.Station .<= 55683,:]
indexy = datp.Station
tm = standfun.(datp.Year, 1960, 2015)
y = (datp.max88 .- mean(datp.max88))./std(datp.max88)
knots = collect(range(0,length = 5, stop = 1))[2:4]
betam0 = initial5(indexy, tm, y, knots, lamv = lamv)

wt = ones(convert(Int,294*293/2))

res0 = GrInd(indexy,tm,y,knots,wt,betam0,lam = 0.3,maxiter = 1000)
res1 = GrFDA(indexy,tm,y,knots,1,wt,betam0,lam = 0.2,maxiter = 1000)




subdat = CSV.read("../doc/subdat.csv")
indexy = subdat.Station[:,1]
tm = standfun.(subdat.Year, 1960, 2015)
y = (subdat.max88 .- mean(subdat.max88))./std(subdat.max88)
knots = collect(range(0,length = 5, stop = 1))[2:4]
betam0 = initial5(indexy, tm, y, knots, lamv = lamv)


wt = ones(convert(Int,85*84/2))

res0 = GrInd(indexy,tm,y,knots,wt,betam0,lam = 0.22,maxiter = 1000)
freqtable(getgroup(res0.deltam,85))

res1 = GrFDA(indexy,tm,y,knots,1,wt,betam0,lam = 0.22,maxiter = 1000)
freqtable(res1.groupest)



subdat1 = CSV.read("../doc/subdat1.csv")
indexy = subdat1.Station[:,1]
tm = standfun.(subdat1.Year, 1960, 2015)
y = (subdat1.max88 .- mean(subdat1.max88))./std(subdat1.max88)
knots = collect(range(0,length = 6, stop = 1))[2:5]
betam0 = initial5(indexy, tm, y, knots, lamv = lamv)


wt = ones(convert(Int,121*120/2))

res0 = GrInd(indexy,tm,y,knots,wt,betam0,lam = 0.22,maxiter = 1000)
freqtable(getgroup(res0.deltam,121))

res1 = GrFDA(indexy,tm,y,knots,1,wt,betam0,lam = 0.2,maxiter = 1000)
freqtable(res1.groupest)


urban = CSV.read("../doc/urban.csv")

indexy = urban.STATE[:,1]
tm = standfun.(urban.year, 0, 13)
y = (urban.ratio .- mean(urban.ratio))./std(urban.ratio)
knots = collect(range(0,length = 5, stop = 1))[2:4]
betam0 = initial5(indexy, tm, y, knots, lamv = lamv)


wt = ones(convert(Int,48*47/2))

res0 = GrInd(indexy,tm,y,knots,wt,betam0,lam = 0.22,maxiter = 1000)
freqtable(getgroup(res0.deltam,48))

res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.21,maxiter = 1000)
freqtable(res1.groupest)


using FreqTables

datp1 = CSV.read("../doc/AggObese1990_2017.csv")
datp1 = datp1[datp1.AGE .!=80,:]

lamv = collect(range(0,20,step=0.5))
indexy1 = datp1.AGE[:,1]
tm1 = standfun.(datp1.IYEAR, 1989, 2018)
y1 = (datp1.PropObese .- mean(datp1.PropObese))./std(datp1.PropObese)
nage = size(unique(datp1.AGE))[1]

datp2  = datp1
datp2.IYEAR = tm1
datp2.PropObese = y1
CSV.write("../doc/datp2.csv", datp2)



knots1 = collect(range(0,length = 5, stop = 1))[2:4]
betam01 = initial5(indexy1, tm1, y1, knots1, lamv = lamv)

wt1 = ones(convert(Int,62*61/2))

res10 = GrInd(indexy1,tm1,y1,knots1,wt1,betam01, lam  = 0.2, maxiter = 1000)
group10 = getgroup(res10.deltam,62)
freqtable(group10)

res11 = GrFDA(indexy1,tm1,y1,knots1,1,wt1,betam01,
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
        BICvec11[l,P] = BICemlog3(res11l)
    end
end

argmin(BICvec11)
res11 = GrFDA(indexy1,tm1,y1,knots1,1,wt1,betam01,lam = lamvec[30],
K0 = 10)
freqtable(res11.groupest)

BICvec2 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res2l = GrFDA(indexy1,tm1,y1,knots1,P,wt1,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec2[l,P] = BICem4(res2l)
    end
end
argmin(BICvec2)

res2 = GrFDA(indexy1,tm1,y1,knots1,1,wt1,betam01,lam = lamvec[21],
K0 = 10)
freqtable(res2.groupest)


argmin(BICvec11)



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


wt21 =  exp.(0.25.*(1 .- ordervec))

BICvec21 = zeros(nlam,3)
for l = 1:nlam
    for P = 3
        res12l = GrFDA(indexy1,tm1,y1,knots1,P,wt21,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec21[l,P] = BICem4(res12l)
    end
end

argmin(BICvec21)

res21 = GrFDA(indexy1,tm1,y1,knots1,3,wt21,betam01,lam = lamvec[25],
K0 = 10)
freqtable(res21.groupest)



wt22 =  exp.(0.5.*(1 .- ordervec))

BICvec22 = zeros(nlam,3)
for l = 1:nlam
    for P = 3
        res12l = GrFDA(indexy1,tm1,y1,knots1,P,wt22,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec22[l,P] = BICem4(res12l)
    end
end

argmin(BICvec22)

res22 = GrFDA(indexy1,tm1,y1,knots1,3,wt2,betam01,lam = lamvec[32],
K0 = 10)
freqtable(res22.groupest)



wt23 =  exp.(0.75.*(1 .- ordervec))

BICvec23 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res23l = GrFDA(indexy1,tm1,y1,knots1,P,wt23,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec23[l,P] = BICem4(res23l)
    end
end

argmin(BICvec23)

res23 = GrFDA(indexy1,tm1,y1,knots1,1,wt3,betam01,lam = lamvec[12],
K0 = 10)
freqtable(res23.groupest)


wt24 =  exp.((1 .- ordervec))

BICvec24 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res24l = GrFDA(indexy1,tm1,y1,knots1,P,wt24,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec24[l,P] = BICem4(res24l)
    end
end

argmin(BICvec24)
res24 = GrFDA(indexy1,tm1,y1,knots1,1,wt3,betam01,lam = lamvec[17],
K0 = 10)

freqtable(res24.groupest)
