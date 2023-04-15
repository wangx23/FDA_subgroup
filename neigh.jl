
function cmatfun(nr::Int64,nc::Int64)
  # input is the dimension of the grid
  N = nr*nc
  indexmat = zeros(Int64, N, 4)
  indexmat[:,1] = (1:N) .- nr
  indexmat[:,2] = (1:N) .- 1
  indexmat[:,3] = (1:N) .+ 1
  indexmat[:,4] = (1:N) .+ nr

  function fun0(x::Number)
      if x < 0
          x = 0
      elseif x >N
          x = 0
      else
          x = x
      end
      return x
  end

  indexmat = fun0.(indexmat)

  indexmat[(nr.*(1:(nr-1)).+1),2] .= 0
  indexmat[nr.*(1:nr),3] .= 0

  Cmat = zeros(Int64, N, N)
  for i = 1:N
      indexi = indexmat[i,:]
      Cmat[i, indexi[indexi.!=0]] .= 1
  end

  return Cmat
end


function calorder(index::Number, Cmat::Array)
    ord = Cmat[index,:]
    k = 0
    while sum(ord.==0)!=1
        k = k + 1
        indexknd = findall(sum(Cmat[:,findall(ord.==k)],dims = 2)[:,1] .>0)
        indexvec = vcat(index,findall(ord.!=0))
        lock = indexknd[findall([x in indexvec for x in indexknd].==false)]
        ord[lock] .= k + 1
    end
    return ord
end
