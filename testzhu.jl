##### export to csv and use zhu's code to test the results ####
using CSV
include("simdat3.jl")

m = 20
ncl = 50
sig2 = 0.1
lamj = [0.1,0.2]

data = simdat3(sig2, lamj, m = m, ncl = ncl, seed = 1)
CSV.write("datz1.csv",data)
