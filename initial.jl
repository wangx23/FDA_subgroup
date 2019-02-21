##### calculate initial value ####
include("Bsplinestd.jl")
using RCall
@rlibrary Spgr

function initial(;indexy::Vector, tm::Vector, y::Vector,
    knots::Vector; g::Int = 20, boundary::Vector = [0,1],
    lam = 0.001)

    Bmt = Bsplinestd(tm,knots,g = g, boundary = boundary)

    betam0 = cal_initialrx(indexy, y, Bmt, 10)
