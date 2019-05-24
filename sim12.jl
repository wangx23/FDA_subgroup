##### a function to simulate and return results ###
@everywhere include("initial.jl")
@everywhere include("scad.jl")
@everywhere include("GrInd.jl")
@everywhere include("GrFDA.jl")
@everywhere include("GrFDA2.jl")
@everywhere include("GrFDAproxy.jl")
@everywhere include("refitFDA.jl")
@everywhere include("BIC.jl")
@everywhere include("simdat.jl")
@everywhere using Dates
using Distributed


@everywhere function sim12(seed)

    m = 50
    ncl = 50
    sig2 = 0.1
    lamj = [0.1,0.2]

    data = simdat2(sig2, lamj, m = m, ncl = ncl, seed = seed)

    indexy = data.ind
    tm = data.time
    y = data.obs
    knots = collect(range(0,length = 7, stop = 1))[2:6]

    wt = ones(convert(Int,100*99/2))
    group = unique(data[:,1:2])[:,1]

    resor = refitFDA(indexy,tm,y,knots,group,2)
    betaor = transpose(resor.alpm[:,group])

    #betam0 = initial2(indexy, tm, y, knots, lam = 5)
    lamv = collect(range(0,10,step=1))
    betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)

    lamvec = collect(range(0.15,0.45,step = 0.01))
    nlam = length(lamvec)

    t0 = Dates.now()
    BICvec0 = zeros(nlam)
    for l = 1:nlam
        res0l = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[l])
        BICvec0[l] = BICind(res0l,1)
    end

    index0 = argmin(BICvec0)
    res0 = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index0])
    group0 = getgroup(res0.deltam,100)
    ng0 = size(unique(group0))[1]
    ari0 = randindex(group,group0)[1]
    norm0 = norm(res0.beta - betaor)

    ## EM
    t1 = Dates.now()
    BICvec1 = zeros(nlam,3)
    for l = 1:nlam
        for P = 1:3
            res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0v5,lam = lamvec[l],
            K0 = 10,maxiter = 1000)
            BICvec1[l,P] = BICem(res1l,1)
        end
    end

    index1 = argmin(BICvec1)

    res1 = GrFDA(indexy,tm,y,knots,index1[2],wt,betam0v5,lam = lamvec[index1[1]],
    K0=10,maxiter = 1000)
    group1 = getgroup(res1.deltam,100)
    ng1 = size(unique(group1))[1]
    ari1 = randindex(group,group1)[1]
    norm1 = norm(res1.beta - betaor)
    estpc1 = index1[2]


    ### iterative steps
    t2 = Dates.now()
    BICvec2 = zeros(nlam,3)
    for l = 1:nlam
        for P = 1:3
            res2l = GrFDA2(indexy,tm,y,knots,P,wt,betam0v5,lam = lamvec[l],
            K0 = 10, maxiter2 = 100, maxiter = 50)
            BICvec2[l,P] = BIC2(res2l,1)
        end
    end

    index2 = argmin(BICvec2)
    res2 = GrFDA2(indexy,tm,y,knots,index2[2],wt,betam0v5,lam = lamvec[index2[1]],
    K0=10, maxiter2 = 100, maxiter = 50)
    group2 = getgroup(res2.deltam,100)
    ng2 = size(unique(group2))[1]
    ari2 = randindex(group,group2)[1]
    norm2 = norm(res2.beta - betaor)
    estpc2 = index2[2]


    ### proxy
    t3 = Dates.now()
    BICvec3 = zeros(nlam)
    for l = 1:nlam
        res3l = GrFDAproxy(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[l], maxiter =1000)
        BICvec3[l] = BICproxy(res3l,1)
    end

    index3 = argmin(BICvec3)

    res3 = GrFDAproxy(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index3], maxiter =1000)
    group3 = getgroup(res3.deltam,100)
    ng3 = size(unique(group3))[1]
    ari3 = randindex(group,group3)[1]

    BICvec32 = zeros(3)
    for P = 1:3
        refit3p = refitFDA(indexy, tm, y, knots, group3, P)
        BICvec32[P] = BICrefit(refit3p)
    end

    estpc3 = argmin(BICvec32)
    refit3 = refitFDA(indexy, tm, y, knots, group3, estpc3)
    norm3 = norm(transpose(refit3.alpm[:,group3]) - betaor)

    t4 = Dates.now()

    ts0 = round(t1 - t0, Dates.Second)
    ts1 = round(t2 - t1, Dates.Second)
    ts2 = round(t3 - t2, Dates.Second)
    ts3 = round(t4 - t3, Dates.Second)



    resvec = [ari0, ari1, ari2, ari3, ng0, ng1, ng2, ng3,
    norm0, norm1, norm2, norm3,
    estpc1, estpc2, estpc3, ts0, ts1, ts2, ts3]
    return resvec
end

#res1 = sim1(1)

using DelimitedFiles
resultsim12init5 = pmap(sim12, 1:100)
writedlm("resultsim12init5.csv", resultsim12init5, ',')
