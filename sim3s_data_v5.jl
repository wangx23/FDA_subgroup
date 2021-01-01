#### simulate for spatial version 5 with 30 replicates 

include("simdat3s_v5.jl")
using RCall

m = 30
sig2 = 0.04
lamj = [0.1,0.2]

ngrid = 12
n = ngrid * ngrid
gridm = zeros(ngrid*ngrid, 2)
gridm[:,1] = repeat(1:ngrid, inner= ngrid)
gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)

function f1(x::Vector)
    0.7*x[1] + x[2] - 13
end

function f2(x::Vector)
    0.75*x[1] - x[2]
end

value1 = mapslices(f1,gridm,dims = 2)
value2 = mapslices(f2,gridm,dims = 2)

group = zeros(Int64,n)

group[((value1.<0) .& (value2 .<=0) .& (gridm[:,1] .<7))[:,1]] .= 1
group[((gridm[:,2] .>=7) .& (group .==0))[:,1]] .= 2
group[group.==0] .= 3



data3s30_v5 = zeros(4320,5,100)

for i = 1:100
    datai = simdat3s_v5(sig2, lamj, group, m = m,seed = i)
    data3s30_v5[:,1,i] = datai.group
    data3s30_v5[:,2,i] = datai.ind
    data3s30_v5[:,3,i] = datai.time
    data3s30_v5[:,4,i] = datai.obs
    data3s30_v5[:,5,i] = datai.meanfun
end


using RCall
@rput data3s30_v5
R"saveRDS(data3s30_v5, '../simdata/data3s30_v5.rds')"



####### m = 20 #####
m = 20
sig2 = 0.04
lamj = [0.1,0.2]

ngrid = 12
n = ngrid * ngrid
gridm = zeros(ngrid*ngrid, 2)
gridm[:,1] = repeat(1:ngrid, inner= ngrid)
gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)

function f1(x::Vector)
    0.7*x[1] + x[2] - 13
end

function f2(x::Vector)
    0.75*x[1] - x[2]
end

value1 = mapslices(f1,gridm,dims = 2)
value2 = mapslices(f2,gridm,dims = 2)

group = zeros(Int64,n)

group[((value1.<0) .& (value2 .<=0) .& (gridm[:,1] .<7))[:,1]] .= 1
group[((gridm[:,2] .>=7) .& (group .==0))[:,1]] .= 2
group[group.==0] .= 3


data3s20_v5 = zeros(2880,5,100)

for i = 1:100
    datai = simdat3s_v5(sig2, lamj, group, m = m,seed = i)
    data3s20_v5[:,1,i] = datai.group
    data3s20_v5[:,2,i] = datai.ind
    data3s20_v5[:,3,i] = datai.time
    data3s20_v5[:,4,i] = datai.obs
    data3s20_v5[:,5,i] = datai.meanfun
end


@rput data3s20_v5
R"saveRDS(data3s20_v5, '../simdata/data3s20_v5.rds')"






