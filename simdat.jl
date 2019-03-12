##### simulated data ###
## two groups
using DataFrames
using Random
using StatsBase
using Statistics

function m1(x::Number)
4*(x - 0.5)^2 + 1
end

function m2(x::Number)
    2.5*exp(-25*(x - 0.25)^2) + 2*exp(-50*(x - 0.75)^2)
end

## eigenfunctions
function psi2(x::Number)
    sqrt(2)*sin(2*pi*x)
end

function psi3(x::Number)
    sqrt(2)*cos(2*pi*x)
end



## equal spaced grid
## ncl: number of observations for each cluster, assume same now.
## the output is a tuple, first element is about data column information,
## the second element is about the data, the third element is the group information
function simdat(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)
    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal = 2 * ncl * m

    data = DataFrame(group = repeat(1:2,inner = m*ncl),
             ind = repeat(1:(2*ncl),inner = m),
             time = repeat(tvec, ncl*2),
             obs = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    for i = 1:(ntotal)
        ti = data.time[i]
        groupi = data.group[i]
        vi = sqrt(lamj[1])*randn(1)[1] + sqrt(lamj[2])*randn(1)[1]*psi2(ti) +
                sqrt(lamj[3])*randn(1)[1]*psi3(ti)

        if groupi== 1
            data.obs[i] = m1(ti) + vi
        else
            data.obs[i] = m2(ti) + vi
        end
    end

    data.obs = data.obs .+ rande
    return data
end
