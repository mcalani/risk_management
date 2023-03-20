function nll_lin(theta0,datamat;vec=false)
    beta = theta0[1:(end-1)]
    sig  = abs(theta0[end])+0.000001
    y = datamat[:,1]
    x = datamat[:,2:end]
    res = y - x*beta

    nll = (-0.5)*(-log(2*pi) .- log(sig^2) .-(res.^2/(sig^2)))
    if vec 
        return nll
    else
        return sum(nll)
    end
end