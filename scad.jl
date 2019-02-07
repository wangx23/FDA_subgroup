#### scad function #####
using LinearAlgebra
function sfun(x::Vector,th::Real)
    xn = norm(x,2)
    thval = 1 - th/xn
    thval*(thval >0)*x
end

function scad(x::Vector, lam::Real, nu::Real, gam::Real)
    temp1 = lam/nu
    temp2 = gam * lam
    xn = norm(x,2)

    if xn <= lam + temp1
        z = sfun(x, temp1)
    elseif xn <= temp2 && xn >= lam + temp1
        z = sfun(x, temp2/((gam-1)*nu))/(1- 1/((gam - 1 )*nu))
    else
        z = x
    end
end
