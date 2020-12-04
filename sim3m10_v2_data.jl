##### sim3m10_v2 data array ###


include("simdat3.jl")
using RCall


m = 10
ncl = 50
# 0.04 [0.1,0.2]
sig2 = 0.1
lamj = [0.4,0.3]

data3m10_v2 = zeros(1500,5,100)


for i = 1:100
    datai = simdat3(sig2, lamj, m = m, ncl = ncl, seed = i)
    data3m10_v2[:,1,i] = datai.group
    data3m10_v2[:,2,i] = datai.ind
    data3m10_v2[:,3,i] = datai.time
    data3m10_v2[:,4,i] = datai.obs
    data3m10_v2[:,5,i] = datai.meanfun
end

using JLD
save("../simdata/data3m10_v2.jld", "data3m10_v2", data3m10_v2)

using RCall
@rput data3m10_v2
R"saveRDS(data3m10_v2, '../simdata/data3m10_v2.rds')"
