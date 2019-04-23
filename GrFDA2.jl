#######  two steps iterative algorithms #####
### GrFDA2 uses GrFDA0, estEM and GrFDA1

include("header.jl")


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

    return convert(Array,betam)
end

function estEM(y::Vector, Bmi::Array,betam::Array, theta::Array,
    lamj::Vector, sig2::Number, ntotal::Int, p::Int, P::Int;
    maxiter = 1000, tolem = 1e-3)

    BtB = transpose(Bmi) * Bmi
    ysubm = reshape(y, size(Bmi)[1],size(betam)[1]) - Bmi * transpose(betam)
    n = size(betam)[1]



    mhat = zeros(n,P)# a (n,P) matrix
    Sigma = zeros(P,P)
    residv = zeros(n)

    maxv = 1
    niteration = 0


    for m in 1:maxiter

        sig2old = sig2
        thetaold = theta
        lamjold = lamj

        Laminv = diagm(0=> 1 ./lamj)

        Vi = inv(transpose(theta) * BtB * theta ./sig2
         + Laminv)


        for i = 1:n
            Bty = transpose(Bmi) * ysubm[:,i]

            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] = mi
            Sigma = Sigma +  mi * transpose(mi) + Vi

            # for updating sig2
            residv[i] = sum((ysubm[:,i] - Bmi *theta * mi).^2) +
            tr(Bmi * theta * Vi * transpose(theta) * transpose(Bmi))

        end

        # for updating theta

        for j = 1:P
            InvE = BtB .* (sum(mhat[:,j].^2) + n*Vi[j,j])

            temp = zeros(p,1)
            for j1 = (1:P)[1:P .!=j]
                temp = temp + theta[:,j1] .* sum((mhat[:,j1].*mhat[:,j] .+ Vi[j1,j]))
            end

            BtyE = transpose(Bmi) * (ysubm * mhat[:,j])  - BtB * temp

            theta[:,j] = inv(InvE) * BtyE
        end

        # update sig2
        sig2 = sum(residv)/ntotal

        Sigma = Sigma ./n
        M0 = Symmetric(theta * Sigma * transpose(theta))
        decompm = eigen(M0)
        lamj = decompm.values[end - P + 1:end]
        theta = decompm.vectors[:,end- P+ 1:end]

        #maxv = maximum([maximum(abs.(theta - thetaold)./(abs.(thetaold) .+ delta1)),
        #maximum(abs.(lamj - lamjold)./(abs.(lamjold) .+ delta1)),
        #abs(sig2 - sig2old)/(abs(sig2old) + delta1)])

        #maxv = norm(thetaold - theta) + norm(lamj - lamjold) + abs(sig2 - sig2old)

        maxv = maximum([maximum(abs.(theta - thetaold)),
        maximum(abs.(lamj - lamjold)),abs(sig2 - sig2old)])
        niteration += 1

        if maxv <= tolem
            break
        end
    end

    res = (theta = theta, lamj = lamj, sig2 = sig2, niteration = niteration)
    return res

end

#### a function to find groups based on given theta, lamj and sig2
function GrFDA1(indexy::Vector, uindex::Vector, Bmt::Array, Bmi::Array,
    y::Vector, theta::Array, lamj::Vector, sig2::Number, Dmat::Array, Ip::Array,
    ntotal::Int, p::Int, n::Int, np::Int, lent::Int, npair::Int,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    maxiter::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)


    ### transform matrix
    Lamm = diagm(0 => lamj)
    In = diagm(0 => ones(lent))
    Covm = inv(In + 1/sig2 * (Bmi * theta * Lamm * transpose(theta) * transpose(Bmi)))
    transP = sqrt(Symmetric(Covm))


    ynew = zeros(ntotal)
    Bminew = transP * Bmi
    Bmtnew = zeros(ntotal, p)
    for i = 1:n
        indexi = indexy .== uindex[i]
        ynew[indexi] = transP * y[indexi]
        Bmtnew[indexi,:] = Bminew
    end


    ##  initial values
    betam = betam0
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)
    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam)

    XtXinv = inverseb(indexy,Bmtnew,nu*lent, p, n, np)
    Xty = zeros(np)
    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bminew) * ynew[indexi]
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

    resid = zeros(n)
    for i = 1:n
        indexi = indexy .== uindex[i]
        resid[i] = sum((ynew[indexi] - Bminew * betam[i,:]).^2)
    end

    residsum = sum(resid)

    flag = 0
    if niteration == maxiter
        flag = 1
    end

    res = (beta = convert(Array,betam), deltam = deltam, residsum = residsum,
    niteration = niteration, flag = flag)


    return res
end

##
function GrFDA2(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    boundary::Vector = [0,1], K0 = 20,maxiter2::Int = 100, maxiter= 100,
    tol2 = 1e-3)

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
            i0 += 1
            Dmat[i0,i] = 1
            Dmat[i0,j] = -1
        end
    end

    ### initial assuming independent
    betam = GrFDA0(indexy,uindex,Bmt,Bmi,y,Dmat,Ip,
    ntotal,p,n,np,lent,npair,P,wt,betam0,lam = lam,maxiter = maxiter)

    ## initial value of sig2, theta and lamj
    Cm = zeros(p,p)
    residual = zeros(ntotal)
    Random.seed!(1256)
    reskm = kmeans(transpose(betam0), K0; maxiter = 200)
    groupkm = reskm.assignments
    centerskm = reskm.centers
    for i = 1:n
        indexi = indexy.== uindex[i]
        residual[indexi] = y[indexi] - Bmi * betam0[i,:]
        cv = betam0[i,:] - centerskm[:,groupkm[i]]
        Cm = Cm + cv * transpose(cv)/n
    end

    decomp = eigen(Cm)
    lamj = decomp.values[end - P + 1:end]
    theta = decomp.vectors[:,end- P+ 1:end]
    sig2 = mean(residual.^2)

    ### given betam update sig2, theta and lamj

    maxtol = 1
    niteration = 0
    deltam = zeros(p, npair)
    residsum = 0


    for m = 1:maxiter2

        lamjold = lamj
        thetaold = theta
        sig2old = sig2
        betamold = betam

        resm = estEM(y,Bmi,betam,theta,lamj,sig2,ntotal,p,P, maxiter = maxiter)

        theta = resm.theta
        lamj = resm.lamj
        sig2 = resm.sig2

        resm1 = GrFDA1(indexy,uindex,Bmt,Bmi,y,theta, lamj, sig2,Dmat,Ip,
        ntotal,p,n,np,lent,npair,P,wt,betam,lam = lam,maxiter = maxiter)
        betam = resm1.beta
        deltam = resm1.deltam
        residsum = resm1.residsum

        maxtol = maximum([maximum(abs.(betam - betamold)),
        maximum(abs.(theta - thetaold)),
        maximum(abs.(lamj - lamjold)),abs(sig2 - sig2old)])

         niteration += 1

        if maxtol <= tol2
            break
        end
    end

    flag = 0
    if niteration == maxiter2
        flag = 1
    end


    res2 = (beta = betam, deltam = deltam, theta = theta, lamj = lamj, sig2 = sig2,
    residsum = residsum, maxtol = maxtol, lent = lent,niteration = niteration, flag = flag)
    return res2
end
