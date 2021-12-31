
##### a function to simulate and return results ###
using Distributed
addprocs(48)
@everywhere include("initial.jl")
@everywhere include("scad.jl")
@everywhere include("GrFDA.jl")
@everywhere include("refitFDA.jl")
@everywhere include("BIC.jl")
@everywhere include("simdat4s.jl")
@everywhere include("neigh.jl")
@everywhere include("complam.jl")

@everywhere function sim4s10(seed::Int64)

    m = 10
    sig2 = 0.04
    lamj = [0.1,0.2]

    ngrid = 14
    n = ngrid * ngrid
    gridm = zeros(ngrid*ngrid, 2)
    gridm[:,1] = repeat(1:ngrid, inner= ngrid)
    gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)


    group = zeros(Int64,n)

    group[(gridm[:,1] .<= 7) .& (gridm[:,2].>7)] .= 1
    group[(gridm[:,1] .<= 7) .& (gridm[:,2].<=7)] .= 2
    group[(gridm[:,1] .> 7) .& (gridm[:,2].>7)] .= 3
    group[(gridm[:,1] .>7) .& (gridm[:,2].<=7)] .= 4


    data = simdat4s(sig2, lamj, group, m = m,seed = seed)

    indexy = data.ind
    tm = data.time
    y = data.obs
    #knots = collect(range(0,length = 7, stop = 1))[2:6]
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
    norm1 = norm(res1.betaest - betaor)/sqrt(nobstotal)
    rmse1 = norm(res1.meanfunest- meanfun)/sqrt(nobstotal*10)
    estpc1 = index1[2]


    ### spatial weights
    Cmat = cmatfun(14,14)
    ordermat = zeros(Int64, 14*14 , 14*14)
    for i = 1:(14*14)
        ordermat[:,i] = calorder(i,Cmat)
    end


    alp2 = [0.1, 0.5, 1]
    nalp = length(alp2)

    lamvec = collect(range(0.2,0.8,length= 41))
    nlam = length(lamvec)

    BICvec2 = zeros(nlam,3,nalp)

    for l1 = 1:nalp
    wt2 = exp.(alp2[l1] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
    for l = 1:nlam
        for P = 1:3
            res1l = GrFDA(indexy,tm,y,knots,P,wt2,betam0v5,lam = lamvec[l],
            K0=12,maxiter = 500)
            BICvec2[l,P,l1] = BICem4(res1l)
        end
    end
    end

    index2 = argmin(BICvec2)
    wt2 = exp.(alp2[index2[3]] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
    res2 = GrFDA(indexy,tm,y,knots,index2[2],wt2,betam0v5,
    lam = lamvec[index2[1]],
    K0=12,maxiter = 500)
    group2 = getgroup(res2.deltam,nobstotal)
    ng2 = size(unique(group2))[1]
    ari2 = randindex(group,group2)[1]
    norm2 = norm(res2.betaest - betaor)/sqrt(nobstotal)
    estpc2 = index2[2]
    rmse2 = norm(res2.meanfunest - meanfun)/sqrt(nobstotal*10)

    resvec = [ari1, ari2, ng1, ng2,
    norm1, norm2,rmse1, rmse2,
    estpc1, estpc2]
    return resvec
end

#res1 = sim1(1)

using DelimitedFiles
resultsim4s10 = pmap(sim4s10, 1:100)
writedlm("../resultnew_v4/resultsim4s10.csv", resultsim4s10, ',')

