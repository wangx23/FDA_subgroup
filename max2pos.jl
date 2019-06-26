

function max2pos(x::Vector)
    maxx = argmax(abs.(x))
    x = x .* sign(x[maxx])
    return x
end
