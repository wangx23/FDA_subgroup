##### a function to simulate and return results ###
using Distributed
addprocs(28)
@everywhere include("../initial.jl")
@everywhere include("../scad.jl")
@everywhere include("../GrFDA.jl")
@everywhere include("../refitFDA.jl")
@everywhere include("../BIC.jl")
@everywhere include("sim_datonegroup.jl")
@everywhere include("../neigh.jl")
@everywhere include("../complam.jl")


@everywhere function sim_onegroup_m30(seed::Int64)

    m = 30
    sig2 = 0.04
    lamj = [0.1,0.2]
    
    n = 150
    
    data = sim_datonegroup(sig2, n, lamj, m , seed )
    
    indexy = data.ind
    tm = data.time
    y = data.obs
    knots = collect(range(0,length = 7, stop = 1))[2:6]
    #knots = collect(range(0,length = 5, stop = 1))[2:4]
    
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
    
    
    ## EM
    BICvec1 = zeros(nlam,3)
    for l = 1:nlam
        for P = 1:3
            res1l = GrFDA(indexy,tm,y,knots,P,wt,betam0v5,lam = lamvec[l],
            K0 = 12,maxiter = 500)
            BICvec1[l,P] = BICem4(res1l,1)
        end
    end
    
    index1 = argmin(BICvec1)
    
    res1 = GrFDA(indexy,tm,y,knots,index1[2],wt,betam0v5,
    lam = lamvec[index1[1]],
    K0=12,maxiter = 500)
    group1 = getgroup(res1.deltam,nobstotal)
    ng1 = size(unique(group1))[1]
    ari1 = randindex(group,group1)[1]
    ri1 = randindex(group,group1)[2]
    estpc1 = index1[2]
    
    resvec = [ng1, ari1, ri1, estpc1]
    return resvec
end

#res1 = sim1(1)

using DelimitedFiles
result_onegroup_m30 = pmap(sim_onegroup_m30, 1:100)
writedlm("../../resultnew_v2/onegroup/result_onegroup_m30.csv", result_onegroup_m30, ',')
