##### sim3m30_v2 data array ###


include("simdat3.jl")
using RCall


m = 30
ncl = 50
sig2 = 0.04
lamj = [0.1,0.2]



data3m30_v2 = zeros(4500,5,100)


for i = 1:100
    datai = simdat3(sig2, lamj, m = m, ncl = ncl, seed = i)
    data3m30_v2[:,1,i] = datai.group
    data3m30_v2[:,2,i] = datai.ind
    data3m30_v2[:,3,i] = datai.time
    data3m30_v2[:,4,i] = datai.obs
    data3m30_v2[:,5,i] = datai.meanfun
end


using RCall
@rput data3m30_v2
R"saveRDS(data3m30_v2, '../simdata/data3m30_v2.rds')"
