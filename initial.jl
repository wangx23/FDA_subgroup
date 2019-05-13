##### calculate initial value ####
#include("Bsplinestd.jl")
include("inverseb.jl")
include("orthogonalBsplines.jl")
include("header.jl")


#### original initial value method ####
function initial(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; g::Int = 1000, boundary::Vector = [0,1],
    lam::Number = 0.001, K0::Int = 0)

    ntotal = length(y)
    #Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary)
    Bmt = orthogonalBsplines(tm,knots)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    np = n * p

    XtX = inverseb(indexy,Bmt,lam, p, n, np)

    Xty = zeros(np)

    for i = 1:n
        indexi = indexy.==uindex[i]
        Xty[(i-1)*p + 1:(i*p)] = transpose(Bmt[indexi,:]) * y[indexi]
    end

    betam0 = XtX * Xty

    betam = convert(Array, transpose(reshape(betam0,p,n)))

    if K0 ==0 || K0==n
        return betam
    else
        betamed = zeros(n)
        for i = 1:n
            betamed[i] = median(betam[i,:])
        end
        rankm = ordinalrank(betamed)
        group = ceil.(Int,rankm * K0/n)
        X0 = zeros(ntotal, K0*p)
        for i = 1:K0
            indexi = [x in uindex[group.==i] for x in indexy]
            X0[indexi,(i-1)*p + 1:i*p] = Bmt[indexi,:]
        end
        est0 = inv(transpose(X0) * X0)*transpose(X0)*y
        alphaest0 = transpose(reshape(est0,p,K0))
        betam = alphaest0[group,:]
        return betam
    end
end


# for each i, estimate its own coefficients use pspline
function initial2(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],lam::Number = 0.001)

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids

    Dm1 = zeros(p-1, p)
    for j=1:(p-1)
        Dm1[j,j] = 1
        Dm1[j,j+1] = -1
    end

    Dm2 = Dm1[1:p-2,1:p-1] * Dm1

    Dmlam = lam*transpose(Dm2) * Dm2


    betam = zeros(n,p)

    for i = 1:n
        indexi = indexy.==uindex[i]
        Bmi = Bmt[indexi,:]
        BtB = transpose(Bmi) * Bmi
        betam[i,:] = inv(BtB + Dmlam) * transpose(Bmi) * y[indexi]
    end

    return betam
end


function gcv(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],lam::Number = 0.001)

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)

    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids

    Dm1 = zeros(p-1, p)
    for j=1:(p-1)
        Dm1[j,j] = 1
        Dm1[j,j+1] = -1
    end

    Dm2 = Dm1[1:p-2,1:p-1] * Dm1

    Dmlam = lam*transpose(Dm2) * Dm2

    Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
    Ini = diagm(0 => ones(convert(Int, ntotal/n)))
    Ih = Ini - Hmat
    trh = tr(Ih)

    yhv = zeros(n)
    for i = 1:n
        indexi = indexy.==uindex[i]
        yhi = Ih * y[indexi]
        yhv[i] = transpose(yhi) * yhi
    end

    gcvvalue = sum(yhv) * ntotal/(n^2*trh^2)

    return gcvvalue

end



#### initial value each i with different lam value ###

function initial3(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],
    lamv::Vector = collect(range(0,10,step =1)))

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)
    BtB = transpose(Bmi) * Bmi

    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    nr = length(uniqtm)

    Dm1 = zeros(p-1, p)
    for j=1:(p-1)
        Dm1[j,j] = 1
        Dm1[j,j+1] = -1
    end

    Dm2 = Dm1[1:p-2,1:p-1] * Dm1
    Dm = transpose(Dm2) * Dm2


    betam = zeros(n,p)

    for i = 1:n
        indexi = indexy.==uindex[i]
        yi = y[indexi]
        gcvvi = gcvi(yi, Bmi, Dm, nr, lamv)
        betam[i,:] = inv(BtB + lamv[argmin(gcvvi)]*Dm) * transpose(Bmi) * yi
    end

    return betam
end



function gcvi(yi::Vector, Bmi::Array, Dm::Array,
    nr::Int64,lamv::Vector)

    nlam = length(lamv)
    gcvvec = zeros(nlam)

    for j = 1:nlam
        Dmlam = lamv[j]*Dm
        Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
        Ini = diagm(0 => ones(nr))
        Ih = Ini - Hmat
        trh = tr(Ih)
        yhi = Ih * yi
        gcvvec[j] = (transpose(yhi) * yhi) * nr/(trh^2)
    end

    return gcvvec
end



##### based on BIC
function initial4(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],
    lamv::Vector = collect(range(0,10,step =1)))

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)
    BtB = transpose(Bmi) * Bmi

    p = size(Bmt, 2)
    uindex = unique(indexy) # unique id for each obs
    n = length(uindex) # number of unique ids
    nr = length(uniqtm)

    Dm1 = zeros(p-1, p)
    for j=1:(p-1)
        Dm1[j,j] = 1
        Dm1[j,j+1] = -1
    end

    Dm2 = Dm1[1:p-2,1:p-1] * Dm1
    Dm = transpose(Dm2) * Dm2


    betam = zeros(n,p)

    for i = 1:n
        indexi = indexy.==uindex[i]
        yi = y[indexi]
        bicvi = bici(yi, Bmi, Dm, nr, lamv)
        betam[i,:] = inv(BtB + lamv[argmin(bicvi)]*Dm) * transpose(Bmi) * yi
    end

    return betam
end




function bici(yi::Vector, Bmi::Array, Dm::Array,
    nr::Int64,lamv::Vector)

    nlam = length(lamv)
    bicvec = zeros(nlam)

    for j = 1:nlam
        Dmlam = lamv[j]*Dm
        Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
        Ini = diagm(0 => ones(nr))
        Ih = Ini - Hmat
        trh = tr(Hmat)
        yhi = Ih * yi
        bicvec[j] = log((transpose(yhi) * yhi)/nr)  +  log(nr)/nr*(trh)
    end

    return bicvec
end
