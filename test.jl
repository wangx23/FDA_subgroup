#### test file ###
using CSV
using DelimitedFiles
include("simdat.jl")
include("Bsplinestd.jl")
include("orthogonalBsplines.jl")

sig2 = 0.2
lamj = [0.3,0.2,0.1]

m = 30
ncl = 50
data = simdat(sig2, lamj, m = m, ncl = ncl)


indexy = data.ind
tm = data.time
y = data.obs

gridx = collect(range(0,length = 100 + 2, stop = 1))[2:101]

knots = collect(range(0,length = 7, stop = 1))[2:6]
Bm = Bsplinestd(data.time, knots, g = 100,boundary = [0,1])
Bmo2 = orthogonalBsplines(gridx,knots)
writedlm("Bmo2.txt",Bmo2)

CSV.write("data.csv",data)
writedlm("Bm.txt",Bm)

include("initial.jl")

betam0 = initial(indexy,tm,y,knots)
writedlm("betam0.txt",betam0)

betam1 = initial(indexy,tm,y,knots,K0 = 20)
writedlm("betam1.txt",betam1)

writedlm("detlam.txt",res1.deltam)
writedlm("groupest.txt", groupest)
