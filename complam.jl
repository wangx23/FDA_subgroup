#### for comparing lamj, depends on the number of estimates components###

function complam(lamest::Vector, lam0::Vector, nPest::Int)
    nP0 = length(lam0) ### length of the truth
    lamest = sort(lamest)  ### length of the estimate
    lam0 = sort(lam0)
    rmse = 0 ### rmse of estimating lam
    if nP0 == nPest
        rmse = norm(lamest - lam0)/sqrt(nP0)
    elseif nP0 < nPest
        lam0ext = zeros(nPest)
        lam0ext[(nPest - nP0 + 1):nPest] = lam0
        rmse = norm(lamest - lam0ext)/sqrt(nPest)
    else nP0 > nPest
        lamestext = zeros(nP0)
        lamestext[(nP0 - nPest + 1):nP0] = lamest
        rmse = norm(lamestext - lam0)/sqrt(nP0)
    end
    return rmse
end
