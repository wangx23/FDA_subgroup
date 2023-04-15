
function BICemCn(obj::NamedTuple,c0::Number = 1)
    npinf = size(obj.beta)
    n = npinf[1]
    p = npinf[2]
    lent = obj.lent
    mhat = obj.mhat
    Vi = obj.Vi
    lamj = abs.(obj.lamj)
    
    condvalue = n* sum(diag(Vi) ./lamj)
    for i = 1:n
        condvalue = condvalue + sum((mhat[i,:].^2).*(1 ./lamj))
    end

    
    P = size(obj.theta)[2]
    nconstraints = P*(P+1)/2
    
    group = getgroup(obj.deltam,n)
    ng = size(unique(group))[1]
    Cn =  c0*log(log(n))

    value  = n *lent * log(obj.sig2) + n * sum(log.(abs.(lamj))) + condvalue + 
    Cn*log(n*lent)*(ng*p + P*p - nconstraints)
    
return value
end









