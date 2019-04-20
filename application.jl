#### application ###
using CSV
using Statistics
using Random
include("initial.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrInd.jl")




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

res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.3,maxiter = 100)
res2 = GrFDA2(indexy,tm,y,knots,2,wt,betam0,lam = 0.3,
maxiter2 = 2, maxiter = 100)
groupest = getgroup(res1.deltam, 656)


datp1 = CSV.read("../doc/AggObese1990_2017.csv")
datp1 = datp1[datp1.AGE .!=80,:]

indexy1 = datp1.AGE
tm1 = standfun.(datp1.IYEAR, 1989, 2018)
y1 = (datp1.PropObese .- mean(datp1.PropObese))./std(datp1.PropObese)
knots = collect(range(0,length = 5, stop = 1))[2:4]
betam01 = initial2(indexy1, tm1, y1, knots, lam = 10)

wt1 = ones(convert(Int,62*61/2))

res10 = GrInd(indexy1,tm1,y1,knots,2,wt1,betam01, lam  = 0.3, maxiter = 100)

res11 = GrFDA(indexy1,tm1,y1,knots,2,wt1,betam01,lam = 0.3,maxiter = 100)
group11 = getgroup(res11.deltam,62)

res12 = GrFDA2(indexy1,tm1,y1,knots,2,wt1,betam01,lam = 0.3,
maxiter2 = 10, maxiter = 100)
group12 = getgroup(res12.deltam,62)
res12.lamj
