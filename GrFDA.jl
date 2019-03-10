### Find sub groups in FDA ####
# tm: observation time
# y: observation
# indexy: observation id
# P is the number of dimension of eigenfunctions
include("scad.jl")
#include("Bsplinestd.jl")
include("initial.jl")
include("orthogonalBsplines.jl")


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
wt = ones(convert(Int,100*99/2))



function GrFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int, wt::Vector, betam0::Array;
    lam::Number = 0.5, nu::Number = 1, gam::Number = 3,
    boundary::Vector = [0,1], maxiter::Int = 1000,
    tolabs::Number = 1e-4, tolrel::Number = 1e-2)

    #Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary) # Bspline matrix
    Bmt = orthogonalBsplines(tm, knots)

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
    psim = decomp.vectors[:,end- P+ 1:end]
    #Bmi = Bsplinestd(uniqtm,knots,g = g, boundary = boundary)
    Bmi = orthogonalBsplines(uniqtm,knots)

    theta = zeros(p, P)
    for i = 1:P
        thetai = inv(transpose(Bmi)*Bmi)* transpose(Bmi) * psim[:,i]
        theta[:,i]= thetai/sqrt(sum(thetai.^2))
        #theta[:,i]  = thetai
    end

    sig2 = mean(residual.^2)

    ## initial values
    deltam = zeros(p, npair)
    deltamold = zeros(p, npair)
    betadiff = zeros(p, npair)
    vm = zeros(p, npair)
    deltamold = transpose(Dmat * betam0)
    betam = betam0

    ## define some variables
    mhat = zeros(n,P)# a (n,P) matrix
    #Vhat = zeros(P,P,n)  # a (P,P,n) array
    #Mhat = zeros(P,P,n)
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

    for m = 1:maxiter
        ## expectation
        Laminv = diagm(0=> 1 ./lamj)
        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            BtB = transpose(Bmi) * Bmi
            Bty = transpose(Bmi) * (y[indexi] - Bmi * betam[i,:])
            Vi = inv(transpose(theta) * BtB * theta ./sig2
             + Laminv)
            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] = mi
            #Vhat[:,:,i] = Vi
            #Mhat[:,:,i] = mi * transpose(mi) + Vi
            Sigma = Sigma +  mi * transpose(mi) + Vi

            #global Vii = Vii + Vi

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

    res = (index = uindex, beta = betam, sig2 = sig2,
    deltam = deltam, rvalue = rvalue, svalue = svalue,
    tolpri = tolpri, toldual = toldual, niteration = niteration, flag = flag)

    return res
end


res1 = GrFDA(indexy,tm,y,knots,3,wt,betam0,lam = 0.3,maxiter = 1000)
sort(res1.beta[:,1])

temp = 0
for i = 1:10
    global temp += 1
end
