#### test without mean function part and only one group ###

include("simdat.jl")
include("orthogonalBsplines.jl")
using LinearAlgebra
function simdat1(sig2::Number, lamj::Vector;
    m::Int = 20, ncl::Int = 50, seed::Int = 2228)
    Random.seed!(seed)
    tvec = collect(range(0, length = m + 2, stop = 1))[2:end-1]

    ntotal =  ncl * m

    data = DataFrame(ind = repeat(1:ncl,inner = m),
             time = repeat(tvec, ncl),
             obs = zeros(ntotal))

    rande = sqrt(sig2) * randn(ntotal)

    lamjsqrt = sqrt.(lamj)

    for i = 1:ncl
        #vi = lamjsqrt[1]*randn(1)[1] .+ lamjsqrt[2]*randn(1)[1] * psi2.(tvec) +
        #        lamjsqrt[3]*randn(1)[1] *psi3.(tvec)
       vi = lamjsqrt[1]*randn(1)[1] * psi2.(tvec) +
                        lamjsqrt[2]*randn(1)[1] *psi3.(tvec)
        data.obs[data.ind .== i] = vi
    end

    data.obs = data.obs .+ rande
    return data
end

sig2 = 0.2
lamj = [0.3,0.2]

dat1 = simdat1(sig2,lamj,m=30,ncl=100)

indexy = dat1.ind
tm = dat1.time
y = dat1.obs
knots = collect(range(0,length = 7, stop = 1))[2:6]
res1 = est(indexy, tm, y, knots,2)


function est(indexy::Vector, tm::Vector, y::Vector, knots::Vector,
    P::Int; boundary::Vector = [0,1], maxiter::Int = 1000, tol = 1e-4)

    Bmt = orthogonalBsplines(tm, knots)

    ntotal = length(y)
    p = size(Bmt, 2)

    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids

    uniqtm = unique(tm)
    lent = length(uniqtm)
    e2 = mean(y.^2)

    residualm = reshape(y, lent, n)

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

    thetaold = theta
    thetaold1 = thetaold

    mhat = zeros(n,P)# a (n,P) matrix
    Vhat = zeros(P,P,n)  # a (P,P,n) array
    #Mhat = zeros(P,P,n)
    Sigma = zeros(P,P)
    InvE = zeros(p,p,P)
    BtyE = zeros(p,P)
    residv = zeros(n)

    sig2 = mean(y.^2)
    niteration = 0
    diffth = 1

    for m = 1:maxiter
        #global lamj
        #global theta
        #global sig2
        #global Sigma
        #global niteration
        #global thetaold
        #global diffth
        #global thetaold1

        Laminv = diagm(0=> 1 ./lamj)

        for i = 1:n
            thetaold = theta

            indexi = indexy .== uindex[i]
            Bmi = Bmt[indexi,:]
            BtB = transpose(Bmi) * Bmi
            Bty = transpose(Bmi) * (y[indexi])
            Vi = inv(transpose(theta) * BtB * theta ./sig2
             + Laminv)
            mi = 1/sig2 .* Vi * transpose(theta) * Bty
            mhat[i,:] = mi
            Vhat[:,:,i] = Vi
            #Mhat[:,:,i] = mi * transpose(mi) + Vi
            Sigma = Sigma +  mi * transpose(mi) + Vi

            #global Vii = Vii + Vi

            # for updating sig2
            residv[i] = sum((y[indexi] - Bmi *theta * mi).^2) +
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
        niteration += 1

        diffth = norm(thetaold - theta)

        if diffth <= 1e-4
            break
        end
    end

    res = (sig2 = sig2, theta = theta, thetaold = thetaold, diffth = diffth,
    lamj = lamj, niteration = niteration)


return res
end
