##### BIC #####
include("getgroup.jl")

## For Ind case
function BICind(obj::NamedTuple, c0::Number)
npinf = size(obj.beta)
n = npinf[1]
p = npinf[2]
group = getgroup(obj.deltam,n)
ng = size(unique(group))[1]
Cn = c0 * log(log(n*p))
value = obj.lent*log(obj.sig2) + Cn*log(n)*(ng*p)/n
return value
end


## BIC for general case, based on EM####
function BICem(obj::NamedTuple, c0::Number)
npinf = size(obj.beta)
n = npinf[1]
p = npinf[2]

group = getgroup(obj.deltam,n)
ng = size(unique(group))[1]
Cn = c0 * log(log(n*p))

value = obj.lent * log(obj.sig2) + sum(log.(obj.lamj)) +
sum(transpose(mapslices(mean,obj.postE,dims = 1))./lamj) +
Cn*log(n)*(ng*p)/n

return value
end

### BIC for iterative algorithm ####
function BIC2(obj::NamedTuple, c0::Number)
    npinf = size(obj.beta)
    n = npinf[1]
    p = npinf[2]

    group = getgroup(obj.deltam,n)
    ng = size(unique(group))[1]
    Cn = c0 * log(log(n*p))

    value = obj.lent * log(obj.sig2) + sum(obj.residsum)/obj.sig2/obj.lent/n^2 +
    sum(log.(obj.lamj)) + Cn*log(n)*(ng*p)/n

    return value
end


### BIC for GrFDAproxy ####

function BICproxy(obj::NamedTuple, c0::Number)
    npinf = size(obj.beta)
    n = npinf[1]
    p = npinf[2]

    group = getgroup(obj.deltam,n)
    ng = size(unique(group))[1]
    Cn = c0 * log(log(n*p))

    value = obj.lent * log(obj.residsum/obj.lent/n) + Cn*log(n)*(ng*p)/n

    return value
end
