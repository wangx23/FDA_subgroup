#### application based on obesity data ###
using CSV
using Statistics
using Random
include("initial.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrInd.jl")
include("GrFDAproxy.jl")
include("refitFDA.jl")
include("BIC.jl")

using FreqTables

function standfun(x::Number, lower::Number, upper::Number)
xnew = (x - lower)/(upper - lower)
end


datp1 = CSV.read("../doc/AggObese1990_2017.csv")
datp1 = datp1[datp1.AGE .!=80,:]

lamv = collect(range(0,20,step=0.5))
indexy1 = datp1.AGE[:,1]
tm1 = standfun.(datp1.IYEAR, 1989, 2018)
y1 = (datp1.PropObese .- mean(datp1.PropObese))./std(datp1.PropObese)



######## create a function with the number of knots ####
### input is the number of knots, alp value and lam value

lamvec = collect(range(0.05,0.45,step = 0.01))
alpvec = [0,0.25,0.5,0.75,1]

function estfun(nknots, alpvec, lamvec, lamv, indexy, tm, y)

    knots1 = collect(range(0,length = nknots, stop = 1))[2:(nknots-1)]
    betam01 = initial5(indexy, tm, y, knots1, lamv = lamv)
    nlam = length(lamvec)

    nage = size(unique(indexy))[1]

    wt1 = ones(convert(Int,nage*(nage-1)/2))

    # without covariance
    BICvec0 = zeros(nlam)
    BICvec01 = zeros(nlam)
    for l = 1:nlam
        res0l = GrInd(indexy1,tm1,y1,knots1,wt1,betam01, lam  = lamvec[l])
        BICvec0[l] = BICind2(res0l,1)
        BICvec01[l] = BICind0(res0l)
    end

    res0 = GrInd(indexy,tm,y,knots1,wt1,betam01,
    lam  = lamvec[argmin(BICvec0)], maxiter = 1000)
    group0 = getgroup(res0.deltam,nage)

    #### considering weights? #######

    wmat =  zeros(nage, nage)
    for i in 1:(nage-1)
      for j in (i+1):(nage)
          wmat[i,j] = abs(i - j)
        end
    end

    ind1 = 1
    ordervec = zeros(convert(Int, nage*(nage-1)/2))
    for i in 1:(nage - 1)
        ind2 = convert(Int,(2*nage - i - 1)*i/2)
        ordervec[ind1:ind2] = wmat[i,(i+1):nage]
        ind1 = ind2+1
    end

    nalp = length(alpvec)

    groupw = zeros(Int,nage,nalp)
    BICarray = zeros(nlam, 3, nalp)

    for k in 1:nalp
        wtk =  exp.(alpvec[k].*(1 .- ordervec))
        for l in 1:nlam
            for P = 1:3
                resl = GrFDA(indexy,tm,y,knots1,P,wtk,betam01,lam = lamvec[l],
                K0 = 10)
                BICarray[l,P,k] = BICem4(resl)
            end
        end
        indk = argmin(BICarray[:,:,k])
        resk = GrFDA(indexy,tm,y,knots1,indk[2],wtk,betam01,lam = lamvec[indk[1]],
        K0 = 10)
        groupk = getgroup(resk.deltam,nage)
        groupw[:,k] = groupk
    end

    res = (BICarray = BICarray, groupw = groupw,
    BICvec0 = BICvec0, group0 = group0)
    return res

end

estknots5 = estfun(5, alpvec[4], lamvec, lamv, indexy1, tm1, y1)






argmin(BICvec1)

res1 = GrFDA(indexy1,tm1,y1,knots1,1,wt1,betam01,lam = lamvec[20],
K0 = 10)
freqtable(res1.groupest)




wt21 =  exp.(0.25.*(1 .- ordervec))

BICvec21 = zeros(nlam,3)
for l = 1:nlam
    for P = 3
        res12l = GrFDA(indexy1,tm1,y1,knots1,P,wt21,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec21[l,P] = BICem4(res12l)
    end
end

argmin(BICvec21)

res21 = GrFDA(indexy1,tm1,y1,knots1,3,wt21,betam01,lam = lamvec[23],
K0 = 10)
freqtable(res21.groupest)



wt22 =  exp.(0.5.*(1 .- ordervec))

BICvec22 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res12l = GrFDA(indexy1,tm1,y1,knots1,P,wt22,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec22[l,P] = BICem4(res12l)
    end
end

argmin(BICvec22)

res22 = GrFDA(indexy1,tm1,y1,knots1,3,wt22,betam01,lam = lamvec[33],
K0 = 10)
freqtable(res22.groupest)



wt23 =  exp.(0.75.*(1 .- ordervec))

BICvec23 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res23l = GrFDA(indexy1,tm1,y1,knots1,P,wt23,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec23[l,P] = BICem4(res23l)
    end
end

argmin(BICvec23)

res23 = GrFDA(indexy1,tm1,y1,knots1,1,wt23,betam01,lam = lamvec[20],
K0 = 10)
freqtable(res23.groupest)


wt24 =  exp.((1 .- ordervec))

BICvec24 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res24l = GrFDA(indexy1,tm1,y1,knots1,P,wt24,betam01,lam = lamvec[l],
        K0 = 10)
        BICvec24[l,P] = BICem4(res24l)
    end
end

argmin(BICvec24)
res24 = GrFDA(indexy1,tm1,y1,knots1,1,wt3,betam01,lam = lamvec[17],
K0 = 10)

freqtable(res24.groupest)
