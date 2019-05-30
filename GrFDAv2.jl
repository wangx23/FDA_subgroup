### Find sub groups in FDA ####
### this algorithm is almost the same as GrFDA
## the only difference is Pspline is considered
# tm: observation time
# y: observation
# indexy: observation id
# P is the number of dimension of eigenfunctions
include("header.jl")
include("refitFDAv2.jl")

function GrFDAv2(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam1::Number = 0.5, lam2::Number = 0.5, nu::Number = 1, gam::Number = 3,
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

## psline penalty
    Dm1 = zeros(p-1, p)
    for j=1:(p-1)
        Dm1[j,j] = 1
        Dm1[j,j+1] = -1
    end

    Dm2 = Dm1[1:p-2,1:p-1] * Dm1
    Dm = transpose(Dm2) * Dm2

    inv(transpose(Bmi)* Bmi + 1.5 * Dm)

    temp = vcat(Bmi, sqrt(1.5)*Dm2)

    inv(transpose(temp) * temp) - inv(transpose(Bmi)* Bmi + 1.5 * Dm)


    # initial values of parameters
    # calculate covariance matrix, since this is grid data
    Cm = zeros(p,p)
    residual = zeros(ntotal)
    Random.seed!(1256)
    reskm = kmeans(transpose(betam0), K0; maxiter = 200)
    groupkm = reskm.assignments
    centerskm = reskm.centers

    resf = refitFDAv2(indexy,tm,y,knots,groupkm, P, lam1 = lam1)
    alpm = resf.alpm

    for i = 1:n
        indexi = indexy.== uindex[i]
        residual[indexi] = y[indexi] - Bmi * alpm[:,groupkm[i]]
        cv = betam0[i,:] - alpm[:,groupkm[i]]
        Cm = Cm + cv * transpose(cv)/n
    end

    decomp = eigen(Cm)
    lamj = decomp.values[end - P + 1:end]
    theta = decomp.vectors[:,end- P+ 1:end]

    sig2 = mean(residual.^2)

    ## initial values
    #betam = transpose(centerskm[:,groupkm])
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

    XtXinv = inversebv2(indexy,Bmt,nu*lent, p, n, np, Dm2, lam1)


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


    for m = 1:maxiter
        ## expectation
        lamjold = lamj
        Laminv = diagm(0=> 1 ./lamj)
        Vi = inv(transpose(theta) * BtB * theta ./sig2
         + Laminv)
        ysubm = reshape(y, lent, n) - Bmi * transpose(betam)

        for i = 1:n
            Bty = transpose(Bmi) * ysubm[:,i]
            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] = mi
            #postE[i,:] = mi.^2 + diag(Vi)
            Sigma = Sigma +  mi * transpose(mi) + Vi

            # for updating sig2
            residv[i] = sum((ysubm[:,i] - Bmi *theta * mi).^2) +
            tr(Bmi * theta * Vi * transpose(theta) * transpose(Bmi))
        end

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


        # update theta and lamj
        Sigma = Sigma ./n
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
            deltam[:,i] = scad(deltam[:,i], wt[i]*lam2, nu,gam)
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

    res = (index = uindex, beta = betam, sig2 = sig2, theta = theta, lamj = lamj, lamjold = lamjold,
    residsum = residsum, lent = lent, deltam = deltam, rvalue = rvalue, svalue = svalue,
    tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag)

    return res
end
