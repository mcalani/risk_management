function gradp!(f,x0;h=0.00000001)
    f0 = f(x0) #should return a vector (if you want a variance estimator)
    T  = size(f0,1)
    k  = size(x0,1)

    g  = zeros(T,k)
    for j ∈ 1:k
        if x0[j]>1; # if argument is big enough, compute relative number
            delta1 = ones(k) .+ I(k).*h
            f1     = f(x0.*delta1[:,j])
            g[:,j] = (f1 - f0)./(x0[j]*h)
        else
            delta2 = I(k).*h
            f1     = f(x0 + delta2[:,j])
            g[:,j] = (f1 - f0)./h
        end
    end
    return g
end
