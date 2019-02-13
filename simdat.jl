##### simulated data ###
## two groups

function m1(x::Number)
4*(x - 0.5)^2 + 1
end

function m1(x::Vector)
    n = length(x)
    value = zeros(n)
    for i = 1:n
        value[i] = 4*(x[i] - 0.5)^2 + 1
    end
    return value
end

function m2(x::Number)
    2.5*exp(-25*(x - 0.25)^2) + 2*exp(-50*(x - 0.75)^2)
end


function m2(x::Vector)
    n = length(x)
    value = zeros(n)
    for i = 1:n
        value[i] = 2.5*exp(-25*(x[i] - 0.25)^2) + 2*exp(-50*(x[i] - 0.75)^2)
    end
    return(value)
end

## eigenfunctions

function psi2(x::Number)
    sqrt(2)*sin(2*pi*x)
end

function psi2(x::Vector)
    n = length(x)
    value = zeros(n)
    for i = 1:n
        value[i] = sqrt(2)*sin(2*pi*x[i])
    end
    return value
end


function psi3(x::Number)
    sqrt(2)*cos(2*pi*x)
end

function psi3(x::Vector)
    n = length(x)
    value = zeros(n)
    for i = 1:n
        value[i] = sqrt(2)*cos(2*pi*x[i])
    end
    return value
end

## equal spaced grid
## ncl: number of observations for each cluster, assume same now.
## the output is a tuple, first element is about data column information,
## the second element is about the data, the third element is the group information
function simdat(sig2::Number, lamj::Vector, m::Int, ncl::Int)
    tvec = collect(range(0, length = m, stop = 1))
