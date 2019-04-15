### Find sub groups in FDA ####
# tm: observation time
# y: observation
# indexy: observation id
# P is the number of dimension of eigenfunctions
include("scad.jl")
#include("Bsplinestd.jl")
include("initial.jl")
include("orthogonalBsplines.jl")
include("getgroup.jl")

using Clustering


include("simdat.jl")
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
nu = 1
gam = 3
g = 1000
maxiter = 1000
tolabs = 1e-4
tolrel = 1e-2
betam0  = initial(indexy,tm,y,knots)
betam0 = initial2(indexy, tm, y, knots, lam = 10)
P = 2
wt = ones(convert(Int,100*99/2))

tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]
mean1 = m1.(tvec)
mean2 = m2.(tvec)



function GrFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    boundary::Vector = [0,1], maxiter::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)

    #Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary) # Bspline matrix
    Bmt = orthogonalBsplines(tm, knots)
    #Bmi = orthogonalBsplines(uniqtm,knots)

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
            i0 += 1
            Dmat[i0,i] = 1
            Dmat[i0,j] = -1
        end
    end

    # initial values of parameters
    # calculate covariance matrix, since this is grid data
    betam0bar = mapslices(mean, betam0,dims = 1)
    Cm = zeros(p,p)
    uniqtm = unique(tm)
    lent = length(uniqtm)
    residual = zeros(ntotal)
    Random.seed!(1256)
    reskm = kmeans(transpose(betam0), 20; maxiter = 200)
    groupkm = reskm.assignments
    centerskm = reskm.centers
    for i = 1:n
        indexi = indexy.== uindex[i]
        residual[indexi] = y[indexi] - Bmt[indexi,:] * betam0[i,:]
        #cv = betam0[i,:] - betam0bar[1,:]
        #cv = betam0[i,:]  - fixbetam[i,:]
        cv = betam0[i,:] - centerskm[:,groupkm[i]]
        Cm = Cm + cv * transpose(cv)/n
    end

    decomp = eigen(Cm)
    lamj = decomp.values[end - P + 1:end]
    theta = decomp.vectors[:,end- P+ 1:end]

    sig2 = mean(residual.^2)

    ## initial values
    betam = transpose(centerskm[:,groupkm])
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)
    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam)


    ## define some variables
    mhat = zeros(n,P)# a (n,P) matrix
    Sigma = zeros(P,P)
    InvE = zeros(p,p,P)
    BtyE = zeros(p,P)
    residv = zeros(n)
    #Vii = zeros(P,P)

    XtXinv = inverseb(indexy,Bmt,nu*lent, p, n, np)
    Xty = zeros(np)
    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bmt[indexi,:]) * y[indexi]
    end

    b1 = XtXinv * Xty
    Btheta = zeros(np)
    normbd = zeros(2)

    niteration = 0
    rvalue = 0
    svalue = 0
    tolpri = n
    toldual = n

    lamjold = lamj

    #if i<=n/2
    #    Bty = transpose(Bmi) * (y[indexi] - mean1)
    #else
    #    Bty = transpose(Bmi) * (y[indexi] - mean2)
    #end

    for m = 1:maxiter
        ## expectation
        lamjold = lamj
        Laminv = diagm(0=> 1 ./lamj)

        #theta = fixtheta
        #lamj = fixlamj
        #sig2 = fixsig2

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            BtB = transpose(Bmi) * Bmi
            Bty = transpose(Bmi) * (y[indexi] - Bmi * betam[i,:])


            Vi = inv(transpose(theta) * BtB * theta ./sig2
             + Laminv)
            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] = mi
            Sigma = Sigma +  mi * transpose(mi) + Vi

            # for updating sig2
            residv[i] = sum((y[indexi] - Bmi *betam[i,:] - Bmi *theta * mi).^2) +
            tr(Bmi * theta * Vi * transpose(theta) * transpose(Bmi))


            # for updating theta
            for j = 1:P
                InvE[:,:,j] = InvE[:,:,j] + BtB .* (mi[j]^2 + Vi[j,j])
                temp = zeros(p,1)
                for j1 = (1:P)[1:P .!=j]
                    temp = temp + theta[:,j1] .*(mi[j1]*mi[j] + Vi[j1,j])
                end
                BtyE[:,j] = BtyE[:,j] + Bty .* mi[j] - BtB * temp
            end
        end

        # update sig2
        sig2 = sum(residv)/ntotal


        # update theta and lamj
        Sigma = Sigma ./n
        for j = 1:P
            theta[:,j] = inv(InvE[:,:,j]) * BtyE[:,j]
        end

        M0 = Symmetric(theta * Sigma * transpose(theta))
        decompm = eigen(M0)
        lamj = decompm.values[end - P + 1:end]
        theta = decompm.vectors[:,end- P+ 1:end]

        # update betam, vm, deltam
        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            Btheta[(i-1)*p + 1:(i*p)] = transpose(Bmi) * Bmi *theta * mhat[i,:]
        end

        temp = reshape((deltamold - 1/nu * vm)*Dmat, np,1)
        betanew = b1 - XtXinv * Btheta + nu*lent* XtXinv * temp
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

    res = (index = uindex, beta = betam, sig2 = sig2, theta = theta, lamj = lamj, lamjold = lamjold,
    deltam = deltam, rvalue = rvalue, svalue = svalue,
    tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag)

    return res
end

betam01 = convert(Array,transpose(resg0.alpm[:,group]))
res1 = GrFDA(indexy,tm,y,knots,2,wt,betam0,lam = 0.3,maxiter = 1000)

res2 = GrFDA1(indexy,tm,y,knots,2,wt,betam0,fixtheta, fixlamj, fixsig2,lam = 0.3,maxiter = 1000)

groupest = getgroup(res1.deltam, 100)
plot(groupest)


using Plots
