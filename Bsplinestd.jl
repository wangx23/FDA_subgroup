#### a function to standardized Bspline matrix ###
#### knots: interior knots
#### boundary: boundary points
#### g: number of grid points
#### cubic bsplines without intercept
using LinearAlgebra

function Bsplinestd(knots::Vector, g::Int, boundary::Vector)
    nknots = length(knots)
    L = (boundary[2] - boundary[1])/(g + 1) # interval length
    # grid points
    gridx = collect(range(boundary[1] + L,length = g, stop = boundary[2] - L))

## degree 0
    Bm0 = zeros(g, nknots+ 4*2 - 1)
    xall = zeros(nknots + 2*4)
    xall[1:4] = boundary[1]* ones(4)
    xall[(end-3):end] = boundary[2] * ones(4)
    xall[5:(end-4)] = knots

    for i = 4:(nknots + 4)
        for j = 1:g
            if( xall[i] <= gridx[j] < xall[i+1] && xall[i] < xall[i+1] )
                Bm0[j,i] = 1
            end
        end
    end

    for j = 1:3
        for i = 1:(nknots + 4*2 - 1 - j)
            alpha1 = 0
            alpha2 = 0
            if xall[i+j]!= xall[i]
                alpha1 = (gridx .- xall[i])./(xall[i+j] - xall[i])
            end
            if xall[i+j+1]!= xall[i+1]
                alpha2 = (xall[i+j+1] .-  gridx)./(xall[i+j+1] - xall[i+1])
            end
            Bm0[:,i] = alpha1 .* Bm0[:,i] + alpha2 .* Bm0[:,i+1]
        end
    end

    Bm0 = Bm0[:,2:(nknots+4)] # remove intercept
    Rm = qr(Bm0).R
    Tm = sqrt(g/L) * pinv(Rm)
    Bm = Bm0 * Tm

    return Bm
end
