##### If group information is given refit the model ####
include("header.jl")

function refitFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector, group::Vector,
    P::Int; boundary::Vector = [0,1], K0 = 20, maxiter::Int = 1000, tol = 1e-4)

    Bmt = orthogonalBsplines(tm, knots)
    uniqtm = unique(tm)
    lent = length(uniqtm)
    Bmi = orthogonalBsplines(uniqtm, knots)
    BtB = transpose(Bmi) * Bmi

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
        cv = betam0[i,:] - alpm[:,group[i]]
        #cv = betam0[i,:] - betam0bar[1,:]
        #cv = betam0[i,:] - centerskm[:,groupkm[i]]
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

    niteration = 0
    lamjold = lamj

    for m = 1:maxiter
        lamjold = lamj
        Laminv = diagm(0=> 1 ./lamj)

        Vi = inv(transpose(theta) * BtB * theta ./sig2
         + Laminv)

         ysubm = reshape(y, lent, n) - Bmi * alpm[:,group]

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bty = transpose(Bmi) * (y[indexi] - Bmi * alpm[:,group[i]])

             mi = 1/sig2 .* Vi * transpose(theta) * Bty
             mhat[i,:] = mi
             Sigma = Sigma +  mi * transpose(mi) + Vi
        # for updating sig2
            residv[i] = sum((y[indexi] - Bmi * alpm[:,group[i]] - Bmi *theta * mi).^2) +
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



        sig2 = sum(residv)/ntotal

    # update theta and lamj
        Sigma = Sigma ./n
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

    residf = zeros(n)
    for i = 1:n
        indexi = indexy .== uindex[i]
        residf[i] = sum((y[indexi] - Bmi * alpm[:,group[i]] -
        Bmi * theta * mhat[i,:]).^2)
    end

    residsum = sum(residf)


    flag = 0
    if niteration == maxiter
        flag = 1
    end

    res = (index = uindex, lent = length(unique(tm)), sig2 = sig2, theta = theta, alpm = alpm,
    lamj = lamj, lamjold = lamjold, residsum = residsum,
     niteration = niteration, flag = flag)

    return res
end
