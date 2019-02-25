#### test file ###
using CSV
using DelimitedFiles
include("simdat.jl")
include("Bsplinestd.jl")

sig2 = 0.2
lamj = [0.6,0.3,0.1]

m = 20
ncl = 50
data = simdat(sig2, lamj, m = m, ncl = ncl)


indexy = data.ind
tm = data.time
y = data.obs

Bm = Bsplinestd(data.time, knots, g = 20,boundary = [0,1])
CSV.write("data.csv",data)
writedlm("Bm.txt",Bm)

include("initial.jl")

betam0 = initial(indexy,tm,y,knots)
writedlm("betam0.txt",betam0)

betam1 = initial(indexy,tm,y,knots,K0 = 20)
writedlm("betam1.txt",betam1)
