##### simulated data ###
#### two eigenfunctions####
## two groups
using DataFrames
using Random
using StatsBase
using Statistics

function m11(x::Number)
4*(x - 0.5)^2 + 1
end

function m12(x::Number)
    2.5*exp(-25*(x - 0.25)^2) + 2*exp(-50*(x - 0.75)^2)
end

function m21(x::Number)
    sqrt(2) * sin(4*pi*x)
end

function m22(x::Number)
    exp(-10*(x-0.25)^2)
end

function m1(x::Number)
    cos(2*pi*x)
end

function m2(x::Number)
    1-2*exp(-6x)
end



#x0 = collect(range(0,1,length = 100))
#plot(x0, m1.(x0))
#plot!(x0, m2.(x0))

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

    mean1 = m1.(tvec)
    mean2 = m2.(tvec)

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


    data.obs = data.obs .+ rande
    return data
end


#### use m11 and m12

function simdat1(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)
    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal = 2 * ncl * m

    data = DataFrame(group = repeat(1:2,inner = m*ncl),
             ind = repeat(1:(2*ncl),inner = m),
             time = repeat(tvec, ncl*2),
             obs = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m11.(tvec)
    mean2 = m12.(tvec)

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


    data.obs = data.obs .+ rande
    return data
end



### use m21 and m22

function simdat2(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)
    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal = 2 * ncl * m

    data = DataFrame(group = repeat(1:2,inner = m*ncl),
             ind = repeat(1:(2*ncl),inner = m),
             time = repeat(tvec, ncl*2),
             obs = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    mean1 = m21.(tvec)
    mean2 = m22.(tvec)

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


    data.obs = data.obs .+ rande
    return data
end
