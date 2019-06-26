##### If group information is given refit the model ####
include("header.jl")

function refitFDA(indexy::Vector, tm::Vector, y::Vector, knots::Vector, group::Vector,
    P::Int, betam0::Array; boundary::Vector = [0,1], maxiter::Int = 50, tol = 1e-4)

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

    #betam0 = initial2(indexy, tm, y, knots)
    Cm = zeros(p,p)
    for i=1:n
        cv = betam0[i,:] - alpm[:,group[i]]
        Cm = Cm + cv * transpose(cv)/n
    end

    decomp = eigen(Cm)
    lamj = decomp.values[end - P + 1:end]
    theta = decomp.vectors[:,end- P+ 1:end]
    theta = mapslices(max2pos,theta,dims = 1)

    sig2 = mean((y - Ux * est).^2)

    mhat = zeros(n,P)
    InvE = zeros(p,p,P)
    BtyE = zeros(p,P)
    residv = zeros(n)
    Btheta = zeros(ntotal)
    Sigma = zeros(P,P)

    niteration = 0
    lamjold = lamj
    thetaold = theta
    sig2old = sig2

    normvalue = 1

    for m = 1:maxiter

        lamjold = lamj
        thetaold = theta
        sig2old = sig2

        Laminv = diagm(0=> 1 ./lamj)

        Vi = inv(transpose(theta) * BtB * theta ./sig2
         + Laminv)

         ysubm = reshape(y, lent, n) - Bmi * alpm[:,group]

         fill!(Sigma,0.0)

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bty = transpose(Bmi) * (y[indexi] - Bmi * alpm[:,group[i]])

             mi = 1/sig2 .* Vi * transpose(theta) * Bty
             mhat[i,:] = mi
             Sigma .= Sigma .+  mi * transpose(mi) .+ Vi
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
        theta = mapslices(max2pos,theta,dims = 1)

    ## update est

        for i = 1:n
            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            Btheta[indexi] =  Bmi *theta * mhat[i,:]
        end


        est = inv(transpose(Ux)*Ux)*transpose(Ux)*( y - Btheta)
        alpm = reshape(est, p, ng)

        niteration += 1

        normvalue = sqrt(norm(lamjold -lamj)^2 + norm(thetaold - theta)^2 +
        norm(sig2old - sig2)^2)


        if normvalue <= tol
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
    lamj = lamj, lamjold = lamjold, thetaold = thetaold, sig2old = sig2old,  normvalue = normvalue,
    residsum = residsum, niteration = niteration, flag = flag)

    return res
end
