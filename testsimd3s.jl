### test GrFDA and GrFDA2 ####

include("initial.jl")
include("simdat3s.jl")
include("scad.jl")
include("GrInd.jl")
include("GrFDA.jl")
include("GrFDA2.jl")
include("GrFDAproxy.jl")
include("refitFDA.jl")
include("BIC.jl")
include("neigh.jl")


#### construct a spatial grid map ####
# from left top to right bottom
# 12*12 grid, each 48
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




m = 50
ncl = 50
sig2 = 0.1
lamj = [0.1,0.2]

data3 = simdat3s(sig2, lamj, group, m = m,seed = 10)

indexy3 = data3.ind
tm3 = data3.time
y3 = data3.obs
knots3 = collect(range(0,length = 8, stop = 1))[2:7]
boundary = [0,1]
maxiter = 1000
P = 2
nu = 1
gam = 3
tolabs = 1e-4
tolrel = 1e-2
wt = ones(convert(Int,144*143/2))

lamv = collect(range(0,20,step=0.5))

betam03v5 = initial5(indexy3,tm3, y3, knots3, lamv = lamv)


lamvec = collect(range(0.15,0.4,step = 0.01))
nlam = length(lamvec)




## EM equal weights

BICvec1 = zeros(nlam,3)
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy3,tm3,y3,knots3,P,wt,betam03v5,lam = lamvec[l],
        K0=12,maxiter = 1000)
        BICvec1[l,P] = BICem(res1l)
    end
end

argmin(BICvec1)

res1 = GrFDA(indexy3,tm3,y3,knots3,2,wt,betam03v5,lam = lamvec[13],
K0=12,maxiter = 1000)

group1 = getgroup(res1.deltam,144)
randindex(group,group1)



# unequal weights

Cmat = cmatfun(12,12)
ordermat = zeros(Int64, 12*12 , 12*12)
for i = 1:144
    ordermat[:,i] = calorder(i,Cmat)
end


alp2 = [0.1, 0.5, 1, 2]
nalp = length(alp2)
wt2 = exp.(1 .- ordermat[findall(tril(ordermat).!=0)])


BICvec2 = zeros(nlam,3,nalp)

for l1 = 1:nalp
wt2 = exp.(alp2[l1] .* (1 .- ordermat[findall(tril(ordermat).!=0)]))
for l = 1:nlam
    for P = 1:3
        res1l = GrFDA(indexy3,tm3,y3,knots3,P,wt2,betam03v5,lam = lamvec[l],
        K0=12,maxiter = 1000)
        BICvec2[l,P,l1] = BICem(res1l)
    end
end
end
