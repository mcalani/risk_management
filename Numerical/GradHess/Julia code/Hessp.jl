function hessp!(f,x0;h=0.000001)
    f0 = f(x0) #should return a vector (if you want a variance estimator)
    T  = size(f0,1)
    k  = size(x0,1)
    h2 = h/2

    H  = zeros(k,k) #%will contain the Hessian
    e  = I(k)
    for ii=1:k
        if abs(x0[ii])>100  #% if argument is big enough, compute relative number
          x0P     = x0.*( ones(k) .+  e[:,ii] .*h2 )
          x0N     = x0.*( ones(k) .-  e[:,ii] .*h2 )
          Deltaii = x0[ii].*h
        else
          x0P     = x0 .+  e[:,ii] .*h2
          x0N     = x0 .-  e[:,ii] .*h2
          Deltaii = h
        end

        for jj=1:ii
            if abs(x0[jj])>100 #% if argument is big enough, compute relative number
            x0PP = x0P .* ( ones(k) .+  e[:,jj] .*h2 )
            x0PN = x0P .* ( ones(k) .-  e[:,jj] .*h2 )
            x0NP = x0N .* ( ones(k) .+  e[:,jj] .*h2 )
            x0NN = x0N .* ( ones(k) .-  e[:,jj] .*h2 )
            Delta = Deltaii.*x0[jj].*h
            else
            x0PP = x0P  .+  e[:,jj] .*h2
            x0PN = x0P  .-  e[:,jj] .*h2
            x0NP = x0N  .+  e[:,jj] .*h2
            x0NN = x0N  .-  e[:,jj] .*h2
            Delta = Deltaii.*h
            end
            fPP = f(x0PP)   # % forward,forward(Wan: kiv..small steps to check on the gradient of min values)
            fPN = f(x0PN)   # % forward,backward
            fNP = f(x0NP)   # % backward,forward
            fNN = f(x0NN)   # % backward,backward

            H[ii,jj] = (sum(fPP)-sum(fPN)-sum(fNP)+sum(fNN)) ./ Delta
            H[jj,ii] = H[ii,jj]
        end
    end
    return H
end
