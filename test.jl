#### test file ###

include("simdat.jl")
include("Bsplinestd.jl")

sig2 = 0.2
lamj = [0.6,0.3,0.1]

m = 20
ncl = 50
knots = collect(range(0,length = 9,stop = 1))[2:8]

data = simdat(sig2, lamj, m = m, ncl = ncl)
Bm = Bsplinestd(data.time, knots, g = 20,boundary = [0,1])
CSV.write("data.csv",data)
