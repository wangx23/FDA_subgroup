### Find sub groups in FDA ####
# tm: observation time
# y: observation
# indexy: observation id
include("scad.jl")
include("Bsplinestd.jl")

indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 7, stop = 1))[2:6]
boundary = [0,1]
nu = 1
gam = 3
g = 20
maxiter = 1000
tolabs = 1e-4
tolrel = 1e-2
betam0 


function GrFDA(;indexy::Vector, tm::Vector, y::Vector,
    betam0::Array, lam = 0.5, nu = 1, gam = 3,
    knots::Vector; g::Int = 20, boundary::Vector = [0,1],
    maxiter::UInt = 1000, tolabs::Number = 1e-4, tolrel::Number = 1e-2)


    Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary) # Bspline matrix

    ntotal = length(y)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    np = n * p

    Ip = diagm(0 => ones(p))
    npair = n * (n - 1)/2

    Dmat = zeros(convert(Int,npair), n) # difference matrix
    i0 = 0
    for i = 1:(n - 1)
        for j = (i+1):n
            global i0 += 1
            Dmat[i0,i] = 1
            Dmat[i0,j] = -1
        end
    end

    ## initial values
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)

    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam0)
