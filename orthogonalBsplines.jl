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
end
