#### simulate for spatial version 5 with 30 replicates 

include("sim_onegroup.jl")
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



data_onegroup_m30 = zeros(4320,5,100)

for i = 1:100
    datai = sim_onegroup(sig2, lamj, group, m = m,seed = i)
    data_onegroup_m30[:,1,i] = datai.group
    data_onegroup_m30[:,2,i] = datai.ind
    data_onegroup_m30[:,3,i] = datai.time
    data_onegroup_m30[:,4,i] = datai.obs
    data_onegroup_m30[:,5,i] = datai.meanfun
end


using RCall
@rput data_onegroup_m30
R"saveRDS(data_onegroup_m30, '../../simdata/onegroup/data_onegroup_m30.rds')"



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


data_onegroup_m20 = zeros(2880,5,100)

for i = 1:100
    datai = sim_onegroup(sig2, lamj, group, m = m,seed = i)
    data_onegroup_m20[:,1,i] = datai.group
    data_onegroup_m20[:,2,i] = datai.ind
    data_onegroup_m20[:,3,i] = datai.time
    data_onegroup_m20[:,4,i] = datai.obs
    data_onegroup_m20[:,5,i] = datai.meanfun
end


@rput data_onegroup_m20
R"saveRDS(data_onegroup_m20, '../../simdata/onegroup/data_onegroup_m20.rds')"



###### m = 10 #######
m = 10
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


data_onegroup_m10 = zeros(1440,5,100)

for i = 1:100
    datai = sim_onegroup(sig2, lamj, group, m = m,seed = i)
    data_onegroup_m10[:,1,i] = datai.group
    data_onegroup_m10[:,2,i] = datai.ind
    data_onegroup_m10[:,3,i] = datai.time
    data_onegroup_m10[:,4,i] = datai.obs
    data_onegroup_m10[:,5,i] = datai.meanfun
end


@rput data_onegroup_m10
R"saveRDS(data_onegroup_m10, '../../simdata/onegroup/data_onegroup_m10.rds')"

