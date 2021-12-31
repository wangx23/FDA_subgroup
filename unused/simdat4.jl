##### simulated data ###
#### two eigenfunctions####
## four groups, the same as zhu and qu
### not finished yet
using DataFrames
using Random
using StatsBase
using Statistics


function m41(x::Number)
    cos(2*pi*x)
end

function m42(x::Number)
    1 - 2*exp(-6*x)
end

function m43(x::Number)
    - 1.5*x
end

function m44(x::Number)
    -1.5*x + 1.5
end



#x0 = collect(range(0,1,length = 100))
#plot(x0, m41.(x0))
#plot!(x0, m42.(x0))
#plot!(x0, m43.(x0))
#plot!(x0, m44.(x0))

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
function simdat4(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)

    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal = 4 * ncl * m

    data = DataFrame(group = repeat(1:4,inner = m*ncl),
             ind = repeat(1:(4*ncl),inner = m),
             time = repeat(tvec, ncl*4),
             obs = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m41.(tvec)
    mean2 = m42.(tvec)
    mean3 = m43.(tvec)
    mean4 = m44.(tvec)

    ### group 1
    for i= 1:ncl
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean1 + vi
    end

    ### group 2
    for i = (ncl+1):(2*ncl)
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean2 + vi
    end

    ### group 3

    for i = (2*ncl+1):(3*ncl)
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean3 + vi
    end

    for i = (3*ncl+1):(4*ncl)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean4 + vi
    end


    data.obs = data.obs .+ rande
    return data
end
