### Find sub groups in FDA ####
# time: observation time
# y: observation
# indexy: observation id
include("scad.jl")
include("Bsplinestd.jl")

function GrFDA(;indexy::Vector, y::Vector, time::Vector,
    betam0::Array, lam = 0.5, nu = 1, gam = 3,
    maxiter::UInt = 1000, tolabs::Number = 1e-4, tolrel::Number = 1e-2)

    
