##### a function to simulate and return results ###
#### use BICem4 ####
@everywhere include("initial.jl")
@everywhere include("scad.jl")
@everywhere include("GrInd.jl")
@everywhere include("GrFDA.jl")
@everywhere include("GrFDA2.jl")
@everywhere include("GrFDAproxy.jl")
@everywhere include("refitFDA.jl")
@everywhere include("BIC.jl")
@everywhere include("simdat3.jl")
@everywhere include("complam.jl")
@everywhere using Dates
using Distributed


@everywhere function sim3m10(seed::Int64)

    m = 10
    ncl = 50
    sig2 = 0.1
    lamj = [0.4,0.3]

    data = simdat3(sig2, lamj, m = m, ncl = ncl, seed = seed)

    indexy = data.ind
    tm = data.time
    y = data.obs
    knots = collect(range(0,length = 5, stop = 1))[2:4]

    nobstotal = length(unique(indexy))
    wt = ones(convert(Int,nobstotal*(nobstotal)/2))
    group = unique(data[:,1:2])[:,1]
    meanfun = data.meanfun

    #betam0 = initial2(indexy, tm, y, knots, lam = 5)
    lamv = collect(range(0,20,step=0.5))

    betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)

    resor = refitFDA(indexy,tm,y,knots,group,2, betam0v5)
    betaor = transpose(resor.alpm[:,group])

    lamvec = collect(range(0.2,0.5,step = 0.01))
    nlam = length(lamvec)

    t0 = Dates.now()
    BICvec0 = zeros(nlam)
    BICvec01 = zeros(nlam)
    for l = 1:nlam
        res0l = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[l])
        BICvec0[l] = BICind2(res0l,1)
        BICvec01[l] = BICind0(res0l)
    end

    index0 = argmin(BICvec0)
    res0 = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index0])
    group0 = getgroup(res0.deltam,nobstotal)
    meanfunest0 = res0.meanfunest
    ng0 = size(unique(group0))[1]
    ari0 = randindex(group,group0)[1]
    norm0 = norm(res0.betaest - betaor)/sqrt(150)
    rmse0 =  norm(meanfunest0 - meanfun)/sqrt(150*10)



    t1 = Dates.now()
    index01 = argmin(BICvec01)
    res01 = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index01])
    meanfunest01 = res01.meanfunest
    group01 = getgroup(res01.deltam,nobstotal)
    ng01 = size(unique(group01))[1]
    ari01 = randindex(group,group01)[1]
    norm01 = norm(res01.betaest - betaor)/sqrt(150)
    rmse01 =  norm(meanfunest01 - meanfun)/sqrt(150*10)

    ## EM
    t2 = Dates.now()
    BICvec1 = zeros(nlam,3)
    for l = 1:nlam
        for P = 1:3
            res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0v5,lam = lamvec[l],
            K0 = 12,maxiter = 1000)
            BICvec1[l,P] = BICem4(res1l,1)
        end
    end

    index1 = argmin(BICvec1)

    res1 = GrFDA(indexy,tm,y,knots,index1[2],wt,betam0v5,
    lam = lamvec[index1[1]],
    K0=12,maxiter = 1000)
    meanfunest1 = res1.meanfunest
    group1 = getgroup(res1.deltam,nobstotal)
    ng1 = size(unique(group1))[1]
    ari1 = randindex(group,group1)[1]
    norm1 = norm(res1.betaest - betaor)/sqrt(150)
    estpc1 = index1[2]
    rmse1 =  norm(meanfunest1 - meanfun)/sqrt(150*10)

    ### accuray of lamj
    rmselam = complam(res1.lamj,lamj,index1[2])


    ### proxy
    t3 = Dates.now()
    BICvec3 = zeros(nlam)
    for l = 1:nlam
        res3l = GrFDAproxy(indexy, tm, y, knots, wt, betam0v5,
        lam = lamvec[l], maxiter =1000)
        BICvec3[l] = BICproxy2(res3l,1)
    end

    index3 = argmin(BICvec3)

    res3 = GrFDAproxy(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index3], maxiter =1000)
    group3 = getgroup(res3.deltam,nobstotal)
    ng3 = size(unique(group3))[1]
    ari3 = randindex(group,group3)[1]

    BICvec32 = zeros(3)
    for P = 1:3
        refit3p = refitFDA(indexy, tm, y, knots, group3, P, betam0v5)
        BICvec32[P] = BICrefit(refit3p)
    end

    estpc3 = argmin(BICvec32)
    refit3 = refitFDA(indexy, tm, y, knots, group3, estpc3, betam0v5)
    norm3 = norm(transpose(refit3.alpm[:,group3]) - betaor)/sqrt(150)
    rmse3 =  norm(refit3.meanfunest- meanfun)/sqrt(150*10)


    t4 = Dates.now()

    ts0 = round(t1 - t0, Dates.Second)
    ts2 = round(t3 - t2, Dates.Second)
    ts3 = round(t4 - t3, Dates.Second)



    resvec = vcat([ari0, ari01, ari1, ari3,
    ng0, ng01, ng1, ng3,
    norm0, norm01, norm1, norm3,
    rmse0, rmse01, rmse1, rmse3,
    estpc1, estpc3, ts0, ts2, ts3,
    rmselam], res1.lamj)
    return resvec
end

#res1 = sim1(1)

using DelimitedFiles
resultsim3m10 = pmap(sim3m10, 1:100)
writedlm("resultnew_v2/resultsim3m10.csv", resultsim3m10, ',')
