#### simulate for spatial 4 groups with 30 replicates 

include("simdat4s_v2.jl")
using RCall

m = 30
sig2 = 0.04
lamj = [0.1,0.2]

ngrid = 14
n = ngrid * ngrid
gridm = zeros(ngrid*ngrid, 2)
gridm[:,1] = repeat(1:ngrid, inner= ngrid)
gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)


group = zeros(Int64,n)

group[(gridm[:,1] .<= 7) .& (gridm[:,2].>7)] .= 1
group[(gridm[:,1] .<= 7) .& (gridm[:,2].<=7)] .= 2
group[(gridm[:,1] .> 7) .& (gridm[:,2].>7)] .= 3
group[(gridm[:,1] .>7) .& (gridm[:,2].<=7)] .= 4


data4s30_v2 = zeros(14*14*30,5,100)

for i = 1:100
    datai = simdat4s_v2(sig2, lamj, group, m = m,seed = i)
    data4s30_v2[:,1,i] = datai.group
    data4s30_v2[:,2,i] = datai.ind
    data4s30_v2[:,3,i] = datai.time
    data4s30_v2[:,4,i] = datai.obs
    data4s30_v2[:,5,i] = datai.meanfun
end


using RCall
@rput data4s30_v2
R"saveRDS(data4s30_v2, '../simdata/data4s30_v2.rds')"



####### m = 20 #####
m = 20
sig2 = 0.04
lamj = [0.1,0.2]

ngrid = 14
n = ngrid * ngrid
gridm = zeros(ngrid*ngrid, 2)
gridm[:,1] = repeat(1:ngrid, inner= ngrid)
gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)


group = zeros(Int64,n)

group[(gridm[:,1] .<= 7) .& (gridm[:,2].>7)] .= 1
group[(gridm[:,1] .<= 7) .& (gridm[:,2].<=7)] .= 2
group[(gridm[:,1] .> 7) .& (gridm[:,2].>7)] .= 3
group[(gridm[:,1] .>7) .& (gridm[:,2].<=7)] .= 4


data4s20 = zeros(14*14*20,5,100)

for i = 1:100
    datai = simdat4s(sig2, lamj, group, m = m,seed = i)
    data4s20[:,1,i] = datai.group
    data4s20[:,2,i] = datai.ind
    data4s20[:,3,i] = datai.time
    data4s20[:,4,i] = datai.obs
    data4s20[:,5,i] = datai.meanfun
end


@rput data4s20
R"saveRDS(data4s20, '../simdata/data4s20.rds')"





###### m = 10 #######
m = 10
sig2 = 0.04
lamj = [0.1,0.2]

ngrid = 14
n = ngrid * ngrid
gridm = zeros(ngrid*ngrid, 2)
gridm[:,1] = repeat(1:ngrid, inner= ngrid)
gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)


group = zeros(Int64,n)

group[(gridm[:,1] .<= 7) .& (gridm[:,2].>7)] .= 1
group[(gridm[:,1] .<= 7) .& (gridm[:,2].<=7)] .= 2
group[(gridm[:,1] .> 7) .& (gridm[:,2].>7)] .= 3
group[(gridm[:,1] .>7) .& (gridm[:,2].<=7)] .= 4


data4s10 = zeros(14*14*10,5,100)

for i = 1:100
    datai = simdat4s(sig2, lamj, group, m = m,seed = i)
    data4s10[:,1,i] = datai.group
    data4s10[:,2,i] = datai.ind
    data4s10[:,3,i] = datai.time
    data4s10[:,4,i] = datai.obs
    data4s10[:,5,i] = datai.meanfun
end


@rput data4s10
R"saveRDS(data4s10, '../simdata/data4s10.rds')"



