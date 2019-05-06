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
function BICem(obj::NamedTuple, c0::Number = 1)
npinf = size(obj.beta)
n = npinf[1]
p = npinf[2]

P = size(obj.theta)[2]
nconstraints = P*(P+1)/2

group = getgroup(obj.deltam,n)
ng = size(unique(group))[1]
Cn = c0 * log(log(n*p))

#value = obj.lent * log(obj.sig2) + obj.residsum/n/obj.sig2 +
#Cn*log(n)*(ng*p + P*p - nconstraints)/n

#value = n*obj.lent * log(obj.residsum/n/obj.lent) +
#Cn*log(n)*(ng*p + P*p - nconstraints)

value = n *obj.lent * log(obj.residsum/n/obj.lent) +
Cn*log(n)*(ng*p) + n * (P*p - nconstraints)


#value = n *obj.lent * log(obj.sig2) + Cn*log(n)*(ng*p)

#value = obj.lent * log(obj.sig2) + obj.residsum/n/obj.sig2 +
#Cn*log(n)*(ng*p)/n

return value
end


function BICem2(obj::NamedTuple)
npinf = size(obj.beta)
n = npinf[1]
p = npinf[2]

P = size(obj.theta)[2]
nconstraints = P*(P+1)/2

group = getgroup(obj.deltam,n)
ng = size(unique(group))[1]
Cn =  log(n*p)


value = n *obj.lent * log(obj.residsum/n/obj.lent) +
Cn*log(n)*(ng*p) + n * (P*p - nconstraints)


return value
end




### BIC for iterative algorithm ####
function BIC2(obj::NamedTuple, c0::Number)
    npinf = size(obj.beta)
    n = npinf[1]
    p = npinf[2]

    P = size(obj.theta)[2]
    nconstraints = P*(P-1)/2

    group = getgroup(obj.deltam,n)
    ng = size(unique(group))[1]
    Cn = c0 * log(log(n*p))

    #value = obj.lent * log(obj.sig2) + obj.residsum/n/obj.sig2 +
     #Cn*log(n)*(ng*p + P*p - nconstraints)/n

     value = n *obj.lent * log(obj.residsum/n/obj.lent) +
     Cn*log(n)*(ng*p) + n * (P*p - nconstraints)

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

    value = n * obj.lent * log(obj.residsum/obj.lent/n) +
    Cn*log(n)*(ng*p)

    return value
end



#### GrFDAproxy use refit to calculate theta lamj, need to define a BIC to select the number of P's ###
### this BIC uses the regular BIC, since the number of groups now is a fixed number

function BICrefit(obj::NamedTuple)
    pginf = size(obj.alpm)
    p = pginf[1]

    n = length(obj.index)
    P = size(obj.theta)[2]
    nconstraints = P*(P-1)/2

    #value = obj.lent * log(obj.sig2) + sum(log.(obj.lamj)) +
    #sum(transpose(mapslices(mean,obj.postE,dims = 1))./obj.lamj) +
    #log(n)*(ng*p)/n

    value = n*obj.lent * log(obj.residsum/n/obj.lent) +
    n * (P*p - nconstraints)

end
