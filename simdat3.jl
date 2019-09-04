##### simulated data ###
#### two eigenfunctions####
## three groups
using DataFrames
using Random
using StatsBase
using Statistics


function m31(x::Number)
    sqrt(2) * sin(4*pi*x)
end

function m32(x::Number)
    exp(-10*(x-0.25)^2)
end

function m33(x::Number)
    1.5x - 1
end



#x0 = collect(range(0,1,length = 100))
#plot(x0, m31.(x0))
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
## ncl: number of observations for each cluster, assume same now.
## the output is a tuple, first element is about data column information,
## the second element is about the data, the third element is the group information
function simdat3(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)

    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal = 3 * ncl * m

    data = DataFrame(group = repeat(1:3,inner = m*ncl),
             ind = repeat(1:(3*ncl),inner = m),
             time = repeat(tvec, ncl*3),
             obs = zeros(ntotal),
             meanfun = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m31.(tvec)
    mean2 = m32.(tvec)
    mean3 = m33.(tvec)

    ### group 1
    for i= 1:ncl
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean1 + vi
        data.meanfun[((i-1)*m + 1):(i*m)] = mean1
    end

    ### group 2
    for i = (ncl+1):(2*ncl)
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean2 + vi
        data.meanfun[((i-1)*m + 1):(i*m)] = mean2
    end

    ### group 3

    for i = (2*ncl+1):(3*ncl)
        #vi = sqrt(lamj[1])*randn(1)[1] .+ sqrt(lamj[2])*randn(1)[1]*psi2.(tvec) +
        #        sqrt(lamj[3])*randn(1)[1]*psi3.(tvec)
        vi = sqrt(lamj[1])*randn(1)[1]*psi2.(tvec) +
                sqrt(lamj[2])*randn(1)[1]*psi3.(tvec)
        data.obs[((i-1)*m + 1):(i*m)] = mean3 + vi
        data.meanfun[((i-1)*m + 1):(i*m)] = mean3
    end


    data.obs = data.obs .+ rande
    return data
end
