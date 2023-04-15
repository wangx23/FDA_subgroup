##### calculate initial value ####
#include("Bsplinestd.jl")
include("inverseb.jl")
include("orthogonalBsplines.jl")
include("header.jl")


#### initial value based on gcv penalty term is identiy matrix###

function initial5(indexy::Vector, tm::Vector, y::Vector,
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

    Imp = diagm(0 => ones(p))


    gcvvec = zeros(nlam)
    for j = 1:nlam
        gcvvaluei = 0
        for i = 1:n
            indexi = indexy.==uindex[i]
            yi = y[indexi]
            gcvvaluei = gcvvaluei + gcvi5(yi, Bmi, Imp, nr, lamv[j])
        end
        gcvvec[j] = gcvvaluei
    end


    betam = zeros(n,p)
    for i = 1:n
        indexi = indexy.==uindex[i]
        yi = y[indexi]
        betam[i,:] = inv(BtB + lamv[argmin(gcvvec)]*Imp) * transpose(Bmi) * yi
    end

    return betam
end


function gcvi5(yi::Vector, Bmi::Array, Imp::Array, nr::Int64,lam1::Number)


    Dmlam = lam1*Imp
    Hmat = Bmi*inv(transpose(Bmi) * Bmi + Dmlam) * transpose(Bmi)
    Ini = diagm(0 => ones(nr))
    Ih = Ini - Hmat
    trh = tr(Ih)
    yhi = Ih * yi
    gcvvalue = (transpose(yhi) * yhi) * nr/(trh^2)

    return gcvvalue
end




