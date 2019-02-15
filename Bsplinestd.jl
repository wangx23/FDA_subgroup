#### a function to standardized Bspline matrix ###
#### x: predictable variable
#### knots: interior knots
#### boundary: boundary points
#### g: number of grid points used to standardize Bspine matrix
#### cubic bsplines without intercept
using LinearAlgebra

function Bsplinestd(x::Vector, knots::Vector; g::Int = 20, boundary::Vector = [0,1])
    nknots = length(knots)
    L = (boundary[2] - boundary[1])/(g + 1) # interval length
    nx = length(x)
    # grid points
    gridx = collect(range(boundary[1] + L,length = g, stop = boundary[2] - L))

## degree 3
    Bm0 = zeros(g, nknots+ 4*2 - 1)
    kall = zeros(nknots + 2*4)
    kall[1:4] = boundary[1]* ones(4)
    kall[(end-3):end] = boundary[2] * ones(4)
    kall[5:(end-4)] = knots

    for i = 4:(nknots + 4)
        for j = 1:g
            if( kall[i] <= gridx[j] < kall[i+1] && kall[i] < kall[i+1] )
                Bm0[j,i] = 1
            end
        end
    end


    for j = 1:3
        for i = 1:(nknots + 4*2 - 1 - j)
            alpha1 = 0
            alpha2 = 0
            if kall[i+j]!= kall[i]
                alpha1 = (gridx .- kall[i])./(kall[i+j] - kall[i])
            end
            if kall[i+j+1]!= kall[i+1]
                alpha2 = (kall[i+j+1] .-  gridx)./(kall[i+j+1] - kall[i+1])
            end
            Bm0[:,i] = alpha1 .* Bm0[:,i] + alpha2 .* Bm0[:,i+1]
        end
    end

    Bm0 = Bm0[:,2:(nknots+4)] # remove intercept
    Rm = qr(Bm0).R
    Tm = sqrt(g/L) * pinv(Rm)

    Bmx = zeros(nx, nknots + 4*2 - 1)
    for i = 4:(nknots + 4)
        for j = 1:nx
            if( kall[i] <= x[j] < kall[i+1] && kall[i] < kall[i+1] )
                Bmx[j,i] = 1
            end
        end
    end

    for j = 1:3
        for i = 1:(nknots + 4*2 - 1 - j)
            alpha1 = 0
            alpha2 = 0
            if kall[i+j]!= kall[i]
                alpha1 = (x .- kall[i])./(kall[i+j] - kall[i])
            end
            if kall[i+j+1]!= kall[i+1]
                alpha2 = (kall[i+j+1] .-  x)./(kall[i+j+1] - kall[i+1])
            end
            Bmx[:,i] = alpha1 .* Bmx[:,i] + alpha2 .* Bmx[:,i+1]
        end
    end

    Bmx = Bmx[:,2:(nknots+4)] # remove intercept


    Bm = Bmx * Tm

    return Bm
end
