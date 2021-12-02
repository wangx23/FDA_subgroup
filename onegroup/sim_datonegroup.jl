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



####v6 0.4 v7 0.3

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
function sim_datonegroup(sig2::Number, n::Int, lamj::Vector, m::Int, seed::Int)

    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]
    ntotal = n*m

    data = DataFrame(group = ones(Int64, ntotal),
             ind = repeat(1:n,inner = m),
             time = repeat(tvec, n),
             obs = zeros(ntotal),
             meanfun = zeros(ntotal))

    id1 = unique(data.ind)
    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m31.(tvec)


    for i= 1:n
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[(data.ind .== id1[i])[:,1]] = mean1 + vi
        data.meanfun[(data.ind .== id1[i])[:,1]] = mean1
    end

    data.obs = data.obs .+ rande
    return data
end
