#### new OrthogonalSplineBasis ###
# knots is interior knots

function orthogonalBsplines(x::Vector, knots::Vector;boundary::Vector = [0,1],
	ord::Int = 4)
    nknots = length(knots)
    nx = length(x)

	kall = zeros(nknots + 2*ord)
	kall[1:4] = boundary[1]* ones(4)
    kall[(end-3):end] = boundary[2] * ones(4)
    kall[5:(end-4)] = knots

	n = length(kall)
	q = n - ord
	M = Mfinal(kall,ord)
	M2 = zeros(ord,q, (n - 2*ord + 1))
	for i = ord:(n-ord)
		M2[:,:,(i - ord + 1)] = M[:,:,i] * ISel(i,ord,n)
	end

	Delta = 1 ./Hankel(collect(1:(2*ord-1)),ord,ord)
	di = diff(kall[ord:(n - ord + 1)])
	d = size(M2)
	s = zeros(d[2],d[2])
	for i = 1:d[3]
		 s = s + di[i]*transpose(M2[:,:,i]) * Delta * M2[:,:,i]
	end

	s = round.(s;digits = 15)
	L = inv(cholesky(s).U)

	N = M2
	for i = 1:d[3]
		N[:,:,i] = M2[:,:,i] * L
	end

	Bmo = zeros(nx, n -ord)

	for j = 1:nx

		if x[j] == kall[(n-ord +1)]
			ind = x[j] .<= kall
		else
			ind = x[j] .< kall
		end

		if (all(ind))|(all(.!ind))
			if x[j] == kall[(n-ord+1)]
				return ones(1,ord) * N[:,:,d[3]]
			else
				return zeros(d[2])
			end
		end

		i = (1:n)[ind][1] - 1
		u = (x[j] - kall[i])/(kall[(i+1)] - kall[i])
		U = u.^(0:(ord - 1))

		Bmo[j,:] = transpose(U) * N[:,:,(i-ord+1)]

	end
	return Bmo
end


function ISel(i::Int,k::Int,n::Int)
    hcat(zeros(k,i-k),diagm(0=>ones(k)),zeros(k,n-k-i))
end

function D0(kall::Vector, i::Int, j::Int, k::Int)
	a = kall[i] - kall[j]
	b = kall[j+k-1] - kall[j]
	rtn = ifelse(b>0, a/b, 0)
	return rtn
end

function D1(kall::Vector, i::Int, j::Int, k::Int)
	a = kall[i+1] - kall[i]
	b = kall[j+k-1] - kall[j]
	rtn = ifelse(b>0, a/b, 0)
	return rtn
end


function Mf(kall::Vector, k::Int, i::Int)
	if k==1
		return [1]
	elseif (i<k)|(i>(length(kall) - k))
		return zeros(k,k)
	else
		D0v  = zeros(k-1)
		D1v  = zeros(k-1)
		for k1 = 1:(k-1)
			D0v[k1] = D0(kall,i,i-k+1+k1,k)
			D1v[k1] = D1(kall,i,i-k+1+k1,k)
		end
		return vcat(Mf(kall,k-1,i),zeros(1,k-1)) * (hcat(diagm(0=>1 .- D0v),zeros(k-1)) + hcat(zeros(k-1),diagm(0=> D0v))) +
		     vcat(zeros(1,k-1),Mf(kall,k-1,i)) * (hcat(diagm(0=> -D1v),zeros(k-1)) + hcat(zeros(k-1),diagm(0=> D1v)))
	end
end

function Mfinal(kall::Vector, k::Int)
	Matf = zeros(k,k,length(kall))
	for i = 1:length(kall)
		Matf[:,:,i] = Mf(kall,k,i)
	end
	return Matf
end

function Hankel(x::Vector, nr::Int, nc::Int)
	Z = transpose(x[1:nc]).*ones(nr, nc)
	nlx = length(x)
	for i = 1:(nr - 1)
		Z[i+1,:] = x[(i+1):nlx][1:nc]
	end
	return Z
end
