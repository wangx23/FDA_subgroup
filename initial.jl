##### calculate initial value ####
include("Bsplinestd.jl")
include("inverseb.jl")

function initial(;indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; g::Int = 20, boundary::Vector = [0,1],
    lam = 0.001)

    Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary)
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
