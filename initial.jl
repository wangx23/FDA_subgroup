##### calculate initial value ####
#include("Bsplinestd.jl")
include("inverseb.jl")
include("orthogonalBsplines.jl")

function initial(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; g::Int = 1000, boundary::Vector = [0,1],
    lam::Number = 0.001, K0::Int = 0)

    ntotal = length(y)
    #Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary)
    Bmt = orthogonalBsplines(tm,knots)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    np = n * p

    XtX = inverseb(indexy,Bmt,lam, p, n, np)

    Xty = zeros(np)

    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bmt[indexi,:]) * y[indexi]
    end

    betam0 = XtX * Xty

    betam = convert(Array, transpose(reshape(betam0,p,n)))

    if K0 ==0 || K0==n
        return betam
    else
        betamed = zeros(n)
        for i = 1:n
            betamed[i] = median(betam[i,:])
        end
        rankm = ordinalrank(betamed)
        group = ceil.(Int,rankm * K0/n)
        X0 = zeros(ntotal, K0*p)
        for i = 1:K0
            indexi = [x in uindex[group.==i] for x in indexy]
            X0[indexi,(i-1)*p + 1:i*p] = Bmt[indexi,:]
        end
        est0 = inv(transpose(X0) * X0)*transpose(X0)*y
        alphaest0 = transpose(reshape(est0,p,K0))
        betam = alphaest0[group,:]
        return betam
    end
end
