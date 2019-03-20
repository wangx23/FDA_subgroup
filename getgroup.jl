#### getgroup

function getgroup(deltam::Array, n::Int; tol = 1e-4)
    p = size(deltam)[1]
    b2value = mapslices(norm, deltam, dims = 1)
    b2value[b2value .<= tol] .= 0

    d2 = zeros(n,n)
    for j = 1:n
      indexj1 = convert(Int,(2*n -j)*(j-1)/2 + 1)
      indexj2 = convert(Int,indexj1 + n - j - 1)
      d2[(n - indexj2 + indexj1):n,j] = b2value[indexj1:indexj2]
    end

    d2 = transpose(d2) + d2

    ngj = 1:n
    groupest = zeros(n)
    j = 0

    while (length(ngj) >0)
      global j = j + 1
      gj = (1:n)[d2[ngj[1],:] .==0]
      indexj = [x in gj for x in ngj]
      gj = ngj[indexj]
      global ngj = ngj[.!indexj]
      groupest[gj] = j * ones(length(gj))
    end

    return groupest
end
