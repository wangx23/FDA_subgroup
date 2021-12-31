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

using DataFrames;
datp1 = CSV.File("../doc/AggObese1990_2017.csv")|> DataFrame!
datp1 = datp1[datp1.AGE .!=80,:]

lamv = collect(range(0,20,step=0.5))
indexy1 = datp1.AGE[:,1]
tm1 = standfun.(datp1.IYEAR, 1989, 2018)
y1 = (datp1.PropObese .- mean(datp1.PropObese))./std(datp1.PropObese)
nage = length(unique(indexy1))

y1 = datp1.PropObese./std(datp1.PropObese)

y = y1
tm = tm1
indexy = indexy1


######## create a function with the number of knots ####
### input is the number of knots, alp value and lam value

lamvec = collect(range(0.15,0.6,step = 0.01))
alpvec = [0,0.25,0.5,0.75,1]

knots1 = collect(range(0,length = 6, stop = 1))[2:5]
betam01 = initial5(indexy1, tm1, y1, knots1, lamv = lamv)
nlam = length(lamvec)

## kmeans to estimate the sig2 initially 

#reskm = kmeans(transpose(betam01), 10; maxiter = 200)
#groupkm = reskm.assignments
#centerskm = reskm.centers

#resf = refitFDA(indexy1,tm1,y1,knots1,groupkm, 1, betam01)

#y1 = y1/sqrt(resf.sig2)



function estfun(nknots, alpvec, lamvec, lamv, indexy, tm, y, ordervec)

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

    #### considering weights#######
    nalp = length(alpvec)

    lamvec2 = collect(range(0.1,1,length=51))

    nlam2 = length(lamvec2)
    groupw = zeros(Int,nage,nalp)
    BICarray = zeros(nlam2, 3, nalp)

    lam2 = zeros(nlam2,nalp)
    lam2[:,1] = lamvec2
    #lam2[:,2] = collect(range(2,20,length=60))
    #lam2[:,3] = collect(range(2,20,length=60))

    lam2[:,2] = lamvec2
    lam2[:,3] = lamvec2
    lam2[:,4] = lamvec2
    lam2[:,5] = lamvec2

    grouparray = zeros(nage, nlam2, 3, nalp)

    for k in 1:nalp
        wtk =  exp.(alpvec[k].*(1 .- ordervec))
        for P = 1:3
            for l in 1:nlam2
                resl = GrFDA(indexy,tm,y,knots1,P,wtk,betam01,lam = lam2[l,k],
                K0 = 8)
                BICarray[l,P,k] = BICem4(resl)
                groupl = getgroup(resl.deltam,nage)
                grouparray[:,l,P,k] = groupl
            end
        end
        indk = argmin(BICarray[:,:,k])
        resk = GrFDA(indexy,tm,y,knots1,indk[2],wtk,betam01,lam = lam2[indk[1],k],
        K0 = 8)
        groupk = getgroup(resk.deltam,nage)
        groupw[:,k] = groupk
    end

    res = (BICarray = BICarray, groupw = groupw, grouparray = grouparray,
    BICvec0 = BICvec0, group0 = group0)
    return res

end





##### calculate weights ####

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
    global ind1 = ind2+1
end

estknots5 = estfun(5, alpvec, lamvec, lamv, indexy1, tm1, y1, ordervec)

argmin(estknots5.BICarray)
minimum(estknots5.BICarray)
freqtable(estknots5.groupw[:,3])

freqtable(estknots5.group0)



##### refit ####
knots1 = collect(range(0,length = 5, stop = 1))[2:4]
betam01 = initial5(indexy1, tm1, y1, knots1, lamv = lamv)
nlam = length(lamvec)


refit52 = refitFDA(indexy1,tm1, y1, knots1, estknots5.groupw[:,3], 1,betam01)

meanfunest52 = refit52.meanfunest.*std(datp1.PropObese) .+ mean(datp1.PropObese)

refit52.theta
refit52.eigenfunmat
refit52.sig2
refit52.lamj

using JLD
using JLD2, FileIO
#save("../resultnew_v42/refit52.jld2", refit52)

save("../app_v3/sig2lam2.jld", "lamj", refit52.lamj, "sig2",refit52.sig2)


using DelimitedFiles

writedlm("../app_v3/meanfunest5.txt",meanfunest52)
writedlm("../app_v3/groupest5.txt",estknots5.groupw[:,5])
writedlm("../app_v3/eigenfunmat5.txt", refit52.eigenfunmat)

Bmtm = orthogonalBsplines(unique(tm1), knots1)
writedlm("../app_v3/Bmtm5.txt", Bmtm)
writedlm("../app_v3/theta5.txt", refit52.theta)


#### 
estknots6 = estfun(6, alpvec, lamvec, lamv, indexy1, tm1, y1, ordervec)

argmin(estknots6.BICarray)
minimum(estknots6.BICarray)
freqtable(estknots6.groupw[:,2])


knots1 = collect(range(0,length = 6, stop = 1))[2:5]
betam01 = initial5(indexy1, tm1, y1, knots1, lamv = lamv)
nlam = length(lamvec)


refit6 = refitFDA(indexy1,tm1, y1, knots1, estknots6.groupw[:,2], 1,betam01)
meanfunest6 = refit6.meanfunest .* std(datp1.PropObese) .+ mean(datp1.PropObese)


refit6.theta
refit6.eigenfunmat
refit6.sig2
refit6.lamj


writedlm("../app_v3/meanfunest6.txt",meanfunest6)
writedlm("../app_v3/groupest6.txt",estknots6.groupw[:,2])
writedlm("../app_v3/eigenfunmat6.txt", refit6.eigenfunmat)
writedlm("../app_v3/theta6.txt", refit6.theta)


###### another BIC #####
function estfunv2(nknots, alpvec, lamvec, lamv, indexy, tm, y, ordervec)

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

    #### considering weights#######
    nalp = length(alpvec)

    groupw = zeros(Int,nage,nalp)
    BICarray = zeros(nlam, 3, nalp)

    for k in 1:nalp
        wtk =  exp.(alpvec[k].*(1 .- ordervec))
        for l in 1:nlam
            for P = 1:3
                resl = GrFDA(indexy,tm,y,knots1,P,wtk,betam01,lam = lamvec[l],
                K0 = 10)
                BICarray[l,P,k] = BICem2(resl)
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


estknots5v2 = estfunv2(5, alpvec, lamvec, lamv, indexy1, tm1, y1, ordervec)
argmin(estknots5v2.BICarray)
minimum(estknots5v2.BICarray)
freqtable(estknots5v2.groupw[:,1])
