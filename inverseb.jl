#### calculate inverse (BTB + nu A^TA)

function inverseb(indexy::Vector, Bm::Array, nu::Number,
  p::Int, n::Int, np::Int)

  uindex = unique(indexy)

  Ip = 1/nu * diagm(0 => ones(p))
  nIp = nu * n * diagm(0 => ones(p))
  Bmt = transpose(Bm)
  matinv = zeros(np,np)

  DB = zeros(p,p)
  AB = zeros(np, p)
  mati = zeros(p,p)

  idp1 = 1
  idp2 = p

  for i = 1:n
    indexi = indexy.== uindex[i]
    mati = inv(Bmt[:,indexi] * Bm[indexi,:] + nIp)
    DB = DB + mati
    AB[idp1:idp2,:] = mati
    matinv[idp1:idp2,idp1:idp2] = mati
    idp1 = idp1 + p
    idp2 = idp2 + p
  end

  IB = inv(Ip - DB)
  matinv = matinv + AB * IB * transpose(AB)
  return(matinv)

end
