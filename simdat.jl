##### simulated data ###
## two groups
function m1(x::Number)
    4*(x - 0.5)^2 + 1
end

function m2(x::Number)
    2.5*exp(-25*(x - 0.25)^2) + 2*exp(-50*(x - 0.75)^2)
end

function simdat()
