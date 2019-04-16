#### application ###
using CSV
using Statistics
using Random

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

res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.5,maxiter = 100)
groupest = getgroup(res1.deltam, 100)
