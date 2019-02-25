### Find sub groups in FDA ####
# tm: observation time
# y: observation
# indexy: observation id
# P is the number of dimension of eigenfunctions
include("scad.jl")
include("Bsplinestd.jl")
include("initial.jl")


indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 7, stop = 1))[2:6]
boundary = [0,1]
nu = 1
gam = 3
g = 1000
maxiter = 1000
tolabs = 1e-4
tolrel = 1e-2
betam0  = initial(indexy,tm,y,knots)
P = 3


function GrFDA(;indexy::Vector, tm::Vector, y::Vector, P::Int,
    betam0::Array, lam = 0.5, nu = 1, gam = 3,
    knots::Vector; g::Int = 1000, boundary::Vector = [0,1],
    maxiter::UInt = 1000, tolabs::Number = 1e-4, tolrel::Number = 1e-2)


    Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary) # Bspline matrix

    ntotal = length(y)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    np = n * p

    Ip = diagm(0 => ones(p))
    npair = convert(Int,n * (n - 1)/2)

    Dmat = zeros(npair, n) # difference matrix
    i0 = 0
    for i = 1:(n - 1)
        for j = (i+1):n
            global i0 += 1
            Dmat[i0,i] = 1
            Dmat[i0,j] = -1
        end
    end


    # initial values of parameters
    # calculate covariance matrix, since this is grid data
    uniqtm = unique(tm)
    lent = length(uniqtm)
    residual = zeros(ntotal)
    for i = 1:n
        indexi = indexy.== uindex[i]
        residual[indexi] = y[indexi] - Bmt[indexi,:] * betam0[i,:]
    end

    e2 = mean(residual.^2)
    residualm = reshape(residual, lent, n)

    Rm = zeros(lent, lent)
    for i = 1:lent
        for j = (i+1):lent
            Rm[i,j] = mean(residualm[i,:].*residualm[j,:])
        end
    end

    Rm = Rm + transpose(Rm)
    for i = 1:lent
        Rm[i,i] = e2
    end

    decomp = eigen(Rm)
    lamj = decomp.values[end - P + 1:end]
    psim = decomp.vectors[:,end- P+ 1:end]* sqrt(20)
    Bmi = Bsplinestd(uniqtm,knots,g = g, boundary = boundary)

    theta = zeros(p, P)
    for i = 1:P
        thetai = inv(transpose(Bmi)*Bmi)* transpose(Bmi) * psim[:,i]
        #theta[:,i]= thetai/sqrt(sum(thetai.^2))
        theta[:,i]  = thetai
    end

    sig2 = mean(residual.^2)

    ## initial values
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)

    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam0)


    ## define some variables
    mhat = zeros(n,P)# a (n,P) matrix
    Vhat = zeros(P,P,n)  # a (P,P,n) array

    betam = betam0

    for m = 1:maxiter
        ## expectation
        Laminv = diagm(0=> 1 ./lamj)
        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            Vi = inv(transpose(theta) * transpose(Bmi) * Bmi * theta ./sig2
             + Laminv)
            Vhat[:,:,i] = Vi
            mi = 1/sig2 .* Vi * transpose(theta) * transpose(Bmi) *
            (y[indexi] - Bmi * betam[i,:])
            mhat[i,:] = mi
        end



    end
