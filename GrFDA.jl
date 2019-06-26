### Find sub groups in FDA ####
# tm: observation time
# y: observation
# indexy: observation id
# P is the number of dimension of eigenfunctions
include("header.jl")


function GrFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    boundary::Vector = [0,1], K0 = 20, maxiter::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)

    uniqtm = unique(tm)
    Bmt = orthogonalBsplines(tm, knots)
    Bmi = orthogonalBsplines(uniqtm,knots)


    lent = length(uniqtm)

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
    Cm = zeros(p,p)
    residual = zeros(ntotal)
    Random.seed!(1256)
    reskm = kmeans(transpose(betam0), K0; maxiter = 200)
    groupkm = reskm.assignments
    centerskm = reskm.centers

    resf = refitFDA(indexy,tm,y,knots,groupkm, P, betam0)

    lamj = resf.lamj
    theta = resf.theta
    sig2 = resf.sig2


    ## initial values
    betam = transpose(resf.alpm[:,groupkm])
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)
    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam)


    ## define some variables
    mhat = zeros(n,P)# a (n,P) matrix
    #postE = zeros(n,P)
    Sigma = zeros(P,P)
    residv = zeros(n)
    #Vii = zeros(P,P)

    XtXinv = inverseb(indexy,Bmt,nu*lent, p, n, np)
    Xty = zeros(np)
    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bmi) * y[indexi]
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

    BtB = transpose(Bmi) * Bmi


    temp1 = zeros(p,1)


    for m = 1:maxiter
        ## expectation
        lamjold = lamj
        Laminv = diagm(0=> 1 ./lamj)
        Vi = inv(transpose(theta) * BtB * theta ./sig2
         + Laminv)
        ysubm = reshape(y, lent, n) - Bmi * transpose(betam)

        fill!(Sigma,0.0)

        for i = 1:n
            Bty = transpose(Bmi) * ysubm[:,i]
            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] .= mi
            #postE[i,:] = mi.^2 + diag(Vi)
            Sigma .= Sigma .+  mi * transpose(mi) .+ Vi

            # for updating sig2
            residv[i] = sum((ysubm[:,i] - Bmi *theta * mi).^2) +
            tr(Bmi * theta * Vi * transpose(theta) * transpose(Bmi))
        end

        for j = 1:P
            InvE = BtB .* (sum(mhat[:,j].^2) + n*Vi[j,j])

            fill!(temp1,0.0)
            for j1 = (1:P)[1:P .!=j]
                temp1 .= temp1 .+ theta[:,j1] .* sum((mhat[:,j1].*mhat[:,j] .+ Vi[j1,j]))
            end

            BtyE = transpose(Bmi) * (ysubm * mhat[:,j])  - BtB * temp1

            theta[:,j] = inv(InvE) * BtyE
        end


        # update sig2
        sig2 = sum(residv)/ntotal

        # update theta and lamj
        Sigma = Sigma ./n
        M0 = Symmetric(theta * Sigma * transpose(theta))
        decompm = eigen(M0)
        lamj = decompm.values[end - P + 1:end]
        theta = decompm.vectors[:,end- P+ 1:end]
        theta = mapslices(max2pos,theta,dims = 1)

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
        vm .+= nu * (betadiff - deltam)

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


    residf = zeros(n)
    for i = 1:n
        indexi = indexy .== uindex[i]
        residf[i] = sum((y[indexi] - Bmi * betam[i,:] -
        Bmi * theta * mhat[i,:]).^2)
    end

    residsum = sum(residf)

    flag = 0
    if niteration == maxiter
        flag = 1
    end


    ### recalculate beta as betaest ####

    groupest = getgroup(deltam, n)
    ugroupest = unique(groupest)
    ng = length(ugroupest)

    betaavg = zeros(ng,p)

    for j = 1:ng
        indexj = groupest.==ugroupest[j]
        nj = sum(indexj)
        betaavg[j,:] = mapslices(mean,betam[indexj,:],dims =1)
    end

    betaest = betaavg[groupest,:]

    res = (index = uindex, beta = betam, betaest = betaest, betaavg = betaavg,
    groupest = groupest, sig2 = sig2, theta = theta, lamj = lamj, lamjold = lamjold,
    residsum = residsum, lent = lent, deltam = deltam, rvalue = rvalue, svalue = svalue,
    tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag)

    return res
end
