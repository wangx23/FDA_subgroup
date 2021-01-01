##### simulated data ###
#### two eigenfunctions####
## three  groups with spatial distribution
using DataFrames
using Random
using StatsBase
using Statistics


function m31(x::Number)
    sqrt(2) * sin(4*pi*x) + 1
end

function m32(x::Number)
    2.5*exp(-25*(x-0.25)^2) + 2*exp(-50*(x-0.75)^2)
end

function m33(x::Number)
    sqrt(2) * sin(4*pi*x)
end



#x0 = collect(range(0,1,length = 100))
#plot(1:m, mean1)
#plot!(1:m, data.obs[data.ind.==2])
#plot!(x0, m32.(x0))
#plot!(x0, m33.(x0))


## eigenfunctions
function psi2(x::Number)
    sqrt(2)*sin(2*pi*x)
end

function psi3(x::Number)
    sqrt(2)*cos(2*pi*x)
end



## equal spaced grid
## the output is a tuple, first element is about data column information,
## the second element is about the data, the third element is the group information
## number of time series
function simdat3s_v5(sig2::Number, lamj::Vector, group::Vector;
    m::Int = 20, seed::Int = 2228)

    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    n = length(group)
    ntotal = n*m

    ng1 = sum(group.==1)
    ng2 = sum(group.==2)
    ng3 = sum(group.==3)

    data = DataFrame(group = repeat(group,inner = m),
             ind = repeat(1:n,inner = m),
             time = repeat(tvec, n),
             obs = zeros(ntotal),
             meanfun = zeros(ntotal))

    id1 = unique(data.ind[data.group.==1])
    id2 = unique(data.ind[data.group.==2])
    id3 = unique(data.ind[data.group.==3])

    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m31.(tvec)
    mean2 = m32.(tvec)
    mean3 = m33.(tvec)

    ### group 1
    for i= 1:ng1
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((data.group.==1) .& (data.ind .== id1[i]))[:,1]] = mean1 + vi
        data.meanfun[((data.group.==1) .& (data.ind .== id1[i]))[:,1]] = mean1
    end

    ### group 2
    for i = 1:ng2
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((data.group.==2) .& (data.ind .== id2[i]))[:,1]] = mean2 + vi
        data.meanfun[((data.group.==2) .& (data.ind .== id2[i]))[:,1]] = mean2
    end

    ### group 3

    for i = 1:ng3
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((data.group.==3) .& (data.ind .== id3[i]))[:,1]] = mean3 + vi
        data.meanfun[((data.group.==3) .& (data.ind .== id3[i]))[:,1]] = mean3
    end


    data.obs = data.obs .+ rande
    return data
end
