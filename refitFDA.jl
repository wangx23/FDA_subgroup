##### If group information is given refit the model ####
include("header.jl")

function refitFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector, group::Vector,
    P::Int; boundary::Vector = [0,1], K0 = 20, maxiter::Int = 1000, tol = 1e-4)

    Bmt = orthogonalBsplines(tm, knots)
    ntotal = length(y)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    Ip = diagm(0 => ones(p))

    ng = size(unique(group))[1]
    Ux = zeros(ntotal, ng*p)

# a new matrix for two groups U
    for i = 1:n
        groupi = group[i]
        indexi = indexy.== uindex[i]
        Ux[indexi, (groupi - 1)*p + 1 : groupi*p] = Bmt[indexi,:]
    end

    # initial value of beta

    est = inv(transpose(Ux)*Ux)*transpose(Ux)*y
    alpm = reshape(est,p,ng)

    betam0 = initial2(indexy, tm, y, knots)
    betam0bar = mapslices(mean, betam0,dims = 1)
    Random.seed!(1256)
    reskm = kmeans(transpose(betam0), K0; maxiter = 200)
    groupkm = reskm.assignments
    centerskm = reskm.centers
    Cm = zeros(p,p)
    for i=1:n
        #cv = betam0[i,:] - alpm[:,group[i]]
        #cv = betam0[i,:] - betam0bar[1,:]
        cv = betam0[i,:] - centerskm[:,groupkm[i]]
        Cm = Cm + cv * transpose(cv)/n
    end

    decomp = eigen(Cm)
    lamj = decomp.values[end - P + 1:end]
    theta = decomp.vectors[:,end- P+ 1:end]

    sig2 = mean((y - Ux * est).^2)

    mhat = zeros(n,P)
    Sigma = zeros(P,P)
    InvE = zeros(p,p,P)
    BtyE = zeros(p,P)
    residv = zeros(n)
    Btheta = zeros(ntotal)
    postE = zeros(n,P)

    niteration = 0
    lamjold = lamj

    for m = 1:maxiter
        lamjold = lamj
        Laminv = diagm(0=> 1 ./lamj)

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            BtB = transpose(Bmi) * Bmi
            Bty = transpose(Bmi) * (y[indexi] - Bmi * alpm[:,group[i]])
            Vi = inv(transpose(theta) * BtB * theta ./sig2
             + Laminv)
             mi = 1/sig2 .* Vi * transpose(theta) * Bty
             mhat[i,:] = mi
             postE[i,:] = mi.^ + diag(Vi)
             Sigma = Sigma +  mi * transpose(mi) + Vi

        # for updating sig2
            residv[i] = sum((y[indexi] - Bmi * alpm[:,group[i]] - Bmi *theta * mi).^2) +
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

    ## update est

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            Btheta[indexi] =  Bmi *theta * mhat[i,:]
        end


        est = inv(transpose(Ux)*Ux)*transpose(Ux)*( y - Btheta)
        alpm = reshape(est, p, ng)

        niteration += 1

        if norm(lamjold - lamj) <= tol
            break
        end
    end

    flag = 0
    if niteration == maxiter
        flag = 1
    end

    res = (index = uindex, lent = length(unique(tm)), sig2 = sig2, theta = theta, alpm = alpm,
    lamj = lamj, lamjold = lamjold, postE = postE,
     niteration = niteration, flag = flag)

    return res
end
