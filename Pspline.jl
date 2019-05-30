####### two steps algorithm ###
### the first step is to select lam1
### the first use gcv
include("header.jl")

function Psplinegcv(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],
    lamv::Vector = collect(range(0,10,step =1)))

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)
    BtB = transpose(Bmi) * Bmi
    nlam = length(lamv)

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

    gcvvec = zeros(nlam)
    for j = 1:nlam
        gcvvaluei = 0
        for i = 1:n
            indexi = indexy.==uindex[i]
            yi = y[indexi]
            gcvvaluei = gcvvaluei + gcvi(yi, Bmi, Dm, nr, lamv[j])
        end
        gcvvec[j] = gcvvaluei
    end


    lam1s = lamv[argmin(gcvvec)]

    betam = zeros(n,p)
    for i = 1:n
        indexi = indexy.==uindex[i]
        yi = y[indexi]
        betam[i,:] = inv(BtB + lam1s*Dm) * transpose(Bmi) * yi
    end

    res = (betam = betam, lam1 = lam1s)

    return res
end



function gcvi(yi::Vector, Bmi::Array, Dm::Array,
    nr::Int64,lam1::Number)


    Dmlam = lam1*Dm
    Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
    Ini = diagm(0 => ones(nr))
    Ih = Ini - Hmat
    trh = tr(Ih)
    yhi = Ih * yi
    gcvvalue = (transpose(yhi) * yhi) * nr/(trh^2)

    return gcvvalue
end



### the second use bic
function Psplinebic(indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; boundary::Vector = [0,1],
    lamv::Vector = collect(range(0,10,step =1)))

    ntotal = length(y)
    Bmt = orthogonalBsplines(tm,knots)
    uniqtm = unique(tm)
    Bmi = orthogonalBsplines(uniqtm,knots)
    BtB = transpose(Bmi) * Bmi
    nlam = length(lamv)

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

    bicvec = zeros(nlam)
    for j = 1:nlam
        bicvaluei = 0
        for i = 1:n
            indexi = indexy.==uindex[i]
            yi = y[indexi]
            bicvaluei = bicvaluei + bici(yi, Bmi, Dm, nr, lamv[j])
        end
        bicvec[j] = bicvaluei
    end


    lam1s = lamv[argmin(bicvec)]


    betam = zeros(n,p)
    for i = 1:n
        indexi = indexy.==uindex[i]
        yi = y[indexi]
        betam[i,:] = inv(BtB + lam1s*Dm) * transpose(Bmi) * yi
    end

    res = (betam = betam, lam1 = lam1s)

    return res
end



function bici(yi::Vector, Bmi::Array, Dm::Array,
    nr::Int64,lam::Number)

    Dmlam = lam*Dm
    Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
    Ini = diagm(0 => ones(nr))
    Ih = Ini - Hmat
    trh = tr(Hmat)
    yhi = Ih * yi
    bicvalue = log((transpose(yhi) * yhi)/nr)  +  trh * log(nr)/nr

    return bicvalue
end
