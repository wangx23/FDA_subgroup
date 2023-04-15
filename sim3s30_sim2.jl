
##### a function to simulate and return results ###
#using Distributed
#addprocs(20)
#@everywhere include("initial.jl")
#@everywhere include("scad.jl")
#@everywhere include("GrFDA.jl")
#@everywhere include("refitFDA.jl")
#@everywhere include("BIC.jl")
#@everywhere include("simdat3s_sim2.jl")
#@everywhere include("neigh.jl")

# can use every to repeat multiple times



include("initial.jl")
include("scad.jl")
include("GrFDA.jl")
include("refitFDA.jl")
include("BIC.jl")
include("neigh.jl")

include("simdat3s_sim2.jl")




function sim3s30_sim2(seed::Int64)

    m = 30
    sig2 = 0.04
    lamj = [0.1,0.2]

    ngrid = 12
    n = ngrid * ngrid
    gridm = zeros(ngrid*ngrid, 2)
    gridm[:,1] = repeat(1:ngrid, inner= ngrid)
    gridm[:,2] = repeat(collect(range(ngrid, 1, step = -1)),outer = ngrid)

    function f1(x::Vector)
        0.7*x[1] + x[2] - 13
    end

    function f2(x::Vector)
        0.75*x[1] - x[2]
    end

    value1 = mapslices(f1,gridm,dims = 2)
    value2 = mapslices(f2,gridm,dims = 2)

    group = zeros(Int64,n)

    group[((value1.<0) .& (value2 .<=0) .& (gridm[:,1] .<7))[:,1]] .= 1
    group[((gridm[:,2] .>=7) .& (group .==0))[:,1]] .= 2
    group[group.==0] .= 3


    data = simdat3s_sim2(sig2, lamj, group, m = m,seed = seed)

    indexy = data.ind
    tm = data.time
    y = data.obs
    knots = collect(range(0,length = 7, stop = 1))[2:6]

    nobstotal = length(unique(indexy))
    wt = ones(convert(Int,nobstotal*(nobstotal)/2))
    group = unique(data[:,1:2])[:,1]
    meanfun = data.meanfun


    lamv = collect(range(0,20,step=0.5))
    betam0v5 = initial5(indexy, tm, y, knots, lamv = lamv)

    resor = refitFDA(indexy,tm,y,knots,group,2, betam0v5)
    betaor = transpose(resor.alpm[:,group])

    lamvec = collect(range(0.1,1,step = 0.01))
    nlam = length(lamvec)

    K0 = 12
    res1l = GrFDA(indexy,tm,y,knots,5,wt,betam0v5,lam = 0,
    K0 = K0,maxiter = 1000)

    estpc1 = minimum((1:5)[cumsum(reverse(res1l.lamj))/sum(res1l.lamj) .>=0.95])


    # equal weights
    BICvec1= zeros(nlam)
    for l = 1:nlam
        res1l = GrFDA(indexy,tm,y,knots,estpc1,wt,betam0v5,lam = lamvec[l],
        K0 = K0,maxiter = 500)
        BICvec1[l] = BICemCn(res1l)
    end

    index1 = argmin(BICvec1)

    res1 = GrFDA(indexy,tm,y,knots,estpc1,wt,betam0v5,
    lam = lamvec[index1],K0=K0,maxiter = 500)
    group1 = getgroup(res1.deltam,nobstotal)
    ng1 = size(unique(group1))[1]
    ari1 = randindex(group,group1)[1]
    norm1 = norm(res1.betaest - betaor)/sqrt(150)
    rmse1 = norm(res1.meanfunest- meanfun)/sqrt(150*30)


    ### spatial weights
    Cmat = cmatfun(12,12)
    ordermat = zeros(Int64, 12*12 , 12*12)
    for i = 1:144
        ordermat[:,i] = calorder(i,Cmat)
    end

    alp2 = collect(range(0.05,1,step = 0.05))
    nalp = length(alp2)


    BICvec2 = zeros(nlam,nalp)
    for l1 = 1:nalp
        wt2 = exp.(alp2[l1] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
        for l = 1:nlam
            res1l = GrFDA(indexy,tm,y,knots,estpc1,wt2,betam0v5,lam = lamvec[l],
            K0=K0,maxiter = 500)
            BICvec2[l,l1] = BICemCn(res1l)
        end
    end

    index2 = argmin(BICvec2)
    wt2 = exp.(alp2[index2[2]] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
    res2 = GrFDA(indexy,tm,y,knots,estpc1,wt2,betam0v5,
    lam = lamvec[index2[1]],
    K0=K0,maxiter = 500)
    group2 = getgroup(res2.deltam,nobstotal)
    ng2 = size(unique(group2))[1]
    ari2 = randindex(group,group2)[1]
    norm2 = norm(res2.betaest - betaor)/sqrt(150)
    rmse2 = norm(res2.meanfunest - meanfun)/sqrt(150*30)


    resvec = [ari1, ari2, ng1, ng2,norm1, norm2, rmse1, rmse2, estpc1]

return resvec
end

res1 = sim3s30_sim2(1)


#using DelimitedFiles
#resultsim3s30_sim2 = pmap(sim3s30_sim2, 1:100)
#writedlm("resultsim3s30_sim2.csv", resultsim3s30_sim2, ',')


#using CSV
#CSV.write("data.csv", data)


