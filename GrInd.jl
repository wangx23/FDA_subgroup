### default subgroups ####

function GrInd(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
     wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    maxiter::Int = 1000,
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

    residv = zeros(ntotal)

    for i=1:n
        indexi = indexy.== uindex[i]
        residv[i] = sum((y[indexi] - Bmi * betam[i,:]).^2)
    end

    sig2 = sum(residv)/ntotal

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

    ### mean function estimation
    meanfunest = zeros(ntotal)
    for i=1:n
        indexi = indexy.== uindex[i]
        meanfunest[indexi] = Bmi * betaest[i,:]
    end


    res = (index = uindex, beta = betam, betaest = betaest, betaavg = betaavg,
    sig2 = sig2, lent = lent, meanfunest = meanfunest,
    deltam = deltam, rvalue = rvalue, svalue = svalue,
    tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag)

    return res
end
