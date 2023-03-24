# module Discretization
# export tpm

using Distributions
using LinearAlgebra

# ---------------------------------------------------------------------------- #
#                                      MOM                                     #
# ---------------------------------------------------------------------------- #

function mom(gx,hx,varshock;J=0, method="doubling") #[sigyJ,sigxJ]
    #[sigyJ,sigxJ]=mom(gx,hx,varshock,J, method)
    # Computes the unconditional variance-covariance matrix of x(t) with x(t+J), that is sigxJ=E[x(t)*x(t+J)'], 
    #and the unconditional variance covariaance matrix of y(t) with y(t+J), that is sigyJ=E[y(t)*y(t+J)']
    # where x(t) evolves as
    # x(t+1) = hx x(t) + e(t+1)
    #and y(t) evolves according to 
    # y(t) = gx x(t)
    #where Ee(t)e(t)'=varshock
    #The parameter J can be any integer
    #method =1 : use doubling algorithm
    #method neq 1 : use algebraic method
    #(c) Stephanie Schmitt-Grohe and Martin Uribe, April 18, 1990, renewed January 24, 2000 and August 18, 2001. 

    if method=="doubling"
        #Doubling algorithm
        hx_old   = copy(hx);
        sig_old  = copy(varshock);
        sigx_old = I(size(hx,1));
        diferenz = .1;
        while diferenz > 1e-25;
            sigx     = hx_old * sigx_old * hx_old' + sig_old;

            diferenz = maximum(abs.(sigx-sigx_old));
            sig_old  = hx_old * sig_old * hx_old' + sig_old;
            hx_old   = hx_old * hx_old;
            sigx_old = copy(sigx);
        end 

    elseif method=="kronecker"
        sigx = zeros(size(hx));
        F = kron(hx,hx);
        sigx[:] = (I(size(F,1))-F)\vec(varshock);
    end


    #Get E{x(t)*x(t+J)'}
    sigxJ= (hx.^(-min(0,J))) .* sigx .* ((hx').^(max(0,J)));


    #Get E{y(t)*y(t+J)'}
    sigyJ=real(complex(gx*sigxJ*gx'));
    return [sigyJ,sigxJ]
end



# ---------------------------------------------------------------------------- #
#                                      TPM                                     #
# ---------------------------------------------------------------------------- #

#collect(range(-UB[j], step=2*UB[j] / (N[j]-1), stop=UB[j]))
#collect(range(-UB[j], UB[j], length=N[j]))

function tpm(A,omega,N;T= 1e5 , Tburn=1e+6, n_sig=sqrt(10)) # output = [Pi,S,Xvec]

    #[Pi,S] = tpm(A,omega,N,T,Tburn,  UB);
    # discretizes the AR(1) process:
    #x_t = A x_t-1  + omega e_t,    (1)
    #where x_t is m-by-1, A is m-by-m,  omega is m-by-r, and e is an r-by-1
    #white noise vector with mean zero and identity variance-covariance matrix.
    #Pi is the n-by-n transiton probability matrix of the discretized state.
    #S is an n-by-m matrix. Element (i,j) of S is the discretized value of the j-th element of x_t in state i. 
    #N is an m-by-1 vector. Element i of N indicates the number of grid points in the discretization of element i of x_t. The default value of N is all elements equal to 10. All grids contain equally spaced points. 
    #T is the length of the simulated time series from (1) used to calculate Pi. The default value of T is 1 million.
    #Tburn number of burn-in draws from simulated time series. Default 0.1 million.
    #UB is an m-by-1 vector. Element (i) of UB contains the upper bound of the
    #grid for the i-th element of the vector x_t. The default value is sqrt(10)*std(x_t(i)). 
    #For a derivation of the matrices Pi and S, see ``Finite-State Approximation Of  VAR Processes:  A Simulation Approach'' by Stephanie Schmitt-Grohé and Martín Uribe, July 11, 2010. 
    #(c) Stephanie Schmitt-Grohé and Martín Uribe, July 11, 2010. 

    #omega = [ 0.000642 0.000050 ; 0.000050 0.000167]
    #omega = [ 0.000642 0.000050 ; 0.000050 0.000168]

    m = size(A,1);
    r = size(omega,2);

    # if nargin<3
    # N = 10*ones(m,1);
    # end

    Sigg = mom(I(m),A,omega)[2]; #variance matrix of AR process

    sigg = sqrt.(diag(Sigg)) ; #Unconditional  standard deviation  of AR process

    UB = n_sig*sigg; #upper bound for grid. 

    #Grid 
    V=Any[]
    for j=1:m
        V = push!(V,collect(range(-UB[j], UB[j], length=N[j]))) ;
    end

    n = prod(N); #total number of possible values of the discretized state

    S = zeros(n,m);
    for i in 1:m
        aux1 = Int(n/(N[i]^(m+1-i)))
        aux2 = Int(n/(N[i]^(i)))
        S[:,i] = repeat(  kron(V[i], ones(aux1)), aux2  )
    end

    Pi = zeros(n,n);
    #initialize the  state
    x0 = zeros(m,1); #initialize simulated time series
    xx = repeat(x0',n,1);
    d = sum((S-xx).^2,dims=2);
    ~,i = findmin(d);

    #randn('state',0)
    Xvec = zeros(m,T+Tburn) ;
    for t=1:T+Tburn
        x = A*x0 + rand(MvNormal(zeros(m),omega));
    #Normality is not required. The command randn can be changed to any other random number generator with mean 0 and unit standard deviaiton. 
        Xvec[:,t] = x ;
        xx = repeat(x',n,1);
        d = sum( (S-xx).^2 , dims=2);
        ~,j = findmin(d);

        if t>Tburn
        Pi[i[1],j[1]] = Pi[i[1],j[1]]+1;
        end

        x0 = copy(x);
        i = deepcopy(j);
        
        if mod(t,T/10)==0
            display(t)
        end
    end

    #Trimming
    z=  Bool.(vec(sum(Pi,dims=2).>0))
    Pi = Pi[z,:];
    Pi = Pi[:,z];
    S = S[z,:];

    z   = Bool.(vec(sum(Pi,dims=2).==0))
    Pi[z,:] .= 1;
    for i=1:size(Pi,1)
        Pi[i,:] = Pi[i,:]./sum(Pi[i,:]);
    end

    Pi[z,:] .=0;
    return Pi,S,Xvec
end

# end
# ---------------------------------------------------------------------------- #
#                                    EXAMPLE                                   #
# ---------------------------------------------------------------------------- #
# using Plots

# #EXAMPLE 1
# A = [0.3244   -0.0683
#     0.2617    0.6051 ]
# omega = [0.00310298148143216 -0.000184623011507698;
#             -0.000184623011507698 0.00542375328735066]
# N = vec([8 8])
# T = 1000_000
# Tburn = 0

# Pi,S,Xvec=tpm(A,omega,N;T=T,Tburn=Tburn)

# sum(Pi,dims=2)

# plot(Pi,legend=:none)
# plot(S,legend=:none)


# # m = size(A,1)
# # gx= copy(I(m))
# # hx= copy(A)
# # varshock= copy(omega)
# #mom(gx,hx,varshock)

# #EXAMPLE 2
# T = 1000_000
# Tburn = 0

# A         = [ 0.745676 0.108413 ; 0.029050 0.739263]
# omega         = [0.000287 0.000039  ; 0.000039 0.000107]
# ny     = 12
# nr     = 12

# N=[ny;nr]

# Π,s_grid,~ = tpm(A,omega,N;T=T,Tburn=Tburn)

# sum(Π,dims=2)
# sum(Π)

# plot(Π,legend=:none)
# plot(S,legend=:none)

# n=prod(N)