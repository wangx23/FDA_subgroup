##### a function to simulate and return results ###
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


@everywhere function sim3m10ncl100(seed)

    m = 10
    ncl = 100
    sig2 = 0.04
    lamj = [0.1,0.2]

    data = simdat3(sig2, lamj, m = m, ncl = ncl, seed = seed)

    indexy = data.ind
    tm = data.time
    y = data.obs
    knots = collect(range(0,length = 5, stop = 1))[2:4]

    nobstotal = length(unique(indexy))
    wt = ones(convert(Int,nobstotal*(nobstotal)/2))
    group = unique(data[:,1:2])[:,1]
    meanfun = data.meanfun

    lamv = collect(range(0,20,step=0.5))

    betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)

    resor = refitFDA(indexy,tm,y,knots,group,2, betam0v5)
    betaor = transpose(resor.alpm[:,group])

    #betam0 = initial2(indexy, tm, y, knots, lam = 5)
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
    ng0 = size(unique(group0))[1]
    ari0 = randindex(group,group0)[1]
    norm0 = norm(res0.beta - betaor)/sqrt(100*3)
    rmse0 = norm(res0.meanfunest - meanfun)/sqrt(100*3*10)

    t1 = Dates.now()


    index01 = argmin(BICvec01)
    res01 = GrInd(indexy, tm, y, knots, wt, betam0v5, lam = lamvec[index01])
    group01 = getgroup(res01.deltam,nobstotal)
    ng01 = size(unique(group01))[1]
    ari01 = randindex(group,group01)[1]
    norm01 = norm(res01.beta - betaor)/sqrt(100*3)
    rmse01 = norm(res01.meanfunest - meanfun)/sqrt(100*3*10)

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
    group1 = getgroup(res1.deltam,nobstotal)
    ng1 = size(unique(group1))[1]
    ari1 = randindex(group,group1)[1]
    norm1 = norm(res1.beta - betaor)/sqrt(100*3)
    rmse1 = norm(res1.meanfunest - meanfun)/sqrt(100*3*10)
    estpc1 = index1[2]

    t3 = Dates.now()


    ts0 = round(t1 - t0, Dates.Second)
    ts2 = round(t3 - t2, Dates.Second)

    resvec = [ari0, ari01, ari1,
    ng0, ng01, ng1,
    norm0, norm01, norm1,
    rmse0, rmse01, rmse1,
    estpc1, ts0,ts2,]
    return resvec
end

#res1 = sim1(1)

using DelimitedFiles
resultsim3m10ncl100 = pmap(sim3m10ncl100, 1:100)
writedlm("resultnew_v2/resultsim3m10ncl100.csv", resultsim3m10ncl100, ',')
