####### an algorithm without EM algorithm #####

using Calculus
using LinearAlgebra
using Clustering

include("initial.jl")
include("simdat.jl")
include("scad.jl")


m = 50
ncl = 50
sig2 = 0.1
lamj = [0.1,0.2]

data = simdat(sig2, lamj, m = m, ncl = ncl)

indexy = data.ind
tm = data.time
y = data.obs
knots = collect(range(0,length = 5, stop = 1))[2:4]
boundary = [0,1]
maxiter = 1000
P = 2
nu = 1
gam = 3
tolabs = 1e-4
tolrel = 1e-2
wt = ones(convert(Int,100*99/2))

tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]
mean1 = m1.(tvec)
mean2 = m2.(tvec)

group = unique(data[:,1:2])[:,1]

betam0 = initial2(indexy, tm, y, knots, lam = 10)




ynew = zeros(ntotal)
Bminew = transP * Bmi
Bmtnew = zeros(ntotal, p)
for i = 1:n
    indexi = indexy .== uindex[i]
    ynew[indexi] = transP * y[indexi]
    Bmtnew[indexi,:] = Bminew
end




Bmt = orthogonalBsplines(tm, knots)
uniqtm = unique(tm)
lent = length(uniqtm)
Bmi = orthogonalBsplines(uniqtm,knots)

function GrFDA2(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    boundary::Vector = [0,1], maxiter2::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)

    Bmt = orthogonalBsplines(tm, knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)

    ntotal = length(y)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    np = n * p
    lent = length(uniqtm)

    # define the difference matrix
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

    ### initial assuming independent
    betam = GrFDA0(indexy,uindex,Bmt,Bmi,y,Dmat,Ip,
    ntotal,p,n,np,lent,npair,P,wt,betam0,lam = lam)

    ### given betam update sig2, theta and lamj








end


####  a function to find groups without consideration of any covariance #
function GrFDA0(indexy::Vector, uindex::Vector, Bmt::Array, Bmi::Array,
    y::Vector, Dmat::Array, Ip::Array,
    ntotal::Int, p::Int, n::Int, np::Int, lent::Int, npair::Int,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    maxiter::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)


    ##  initial values
    betam = betam0
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)
    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam)

    XtXinv = inverseb(indexy,Bmt,nu*lent, p, n, np)
    Xty = zeros(np)
    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bmi) * y[indexi]
    end

    b1 = XtXinv * Xty

    normbd = zeros(2)

    niteration = 0
    rvalue = 0
    svalue = 0
    tolpri = n
    toldual = n

    for m = 1:maxiter

        temp = reshape((deltamold - 1/nu * vm)*Dmat, np,1)
        betanew = b1 + nu*lent* XtXinv * temp
        betam = transpose(reshape(betanew, p, n))
        betadiff = transpose(Dmat * betam)

        deltam = betadiff + (1/nu) * vm

        # update deltam
        for i = 1:npair
            deltam[:,i] = scad(deltam[:,i], wt[i]*lam, nu,gam)
        end
        vm =  vm + nu * (betadiff - deltam)

        normbd[1] = norm(betadiff)
        normbd[2] = norm(deltam)

        tolpri = tolabs*sqrt(npair*p) + tolrel* maximum(normbd)
        toldual = tolabs*sqrt(n * p) + tolrel * norm(vm * Dmat)

        rvalue = norm(betadiff - deltam)
        svalue = nu * norm((deltam - deltamold)*Dmat)

        deltamold = deltam
        niteration += 1

       if (rvalue <= tolpri) && (svalue <= toldual)
          break
       end

    end

    flag = 0
    if niteration == maxiter
        flag = 1
    end

    return betam
end

function estEM(Bmi::Array,betam::Array)


end

res2 = GrFDA2(indexy,tm,y,knots, 2, wt, betam0, maxiter = 10,lam =0.4)

groupest = getgroup(res2.deltam, 100)
plot(groupest)
