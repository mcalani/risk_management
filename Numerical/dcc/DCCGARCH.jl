using Optim
using Plots, StatsPlots
using DelimitedFiles
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using Distributions
using Random
using ARCHModels
using LaTeXStrings
using ForwardDiff

# ---------------------------------------------------------------------------- #
#                                   READ DATA                                  #
# ---------------------------------------------------------------------------- #

cd(dirname(@__FILE__))

Rstock = Matrix(DataFrame(CSV.File("Rstockjl.csv")))[:,2:end]
CovQ   = Matrix(DataFrame(CSV.File("CovQjl.csv")))[:,2:end]

T,N=size(Rstock)


CovQ = reshape(CovQ, (N,N,T)) #Implied Covariances
Rstock= convert(Matrix, Rstock) #Returns

#Setting dimensions
d=vec(1:N)
k=copy(N)

# Random.setseed(1)
# j=sample(d, k, replace = false)
j=1:k

#in sample 1500
n1  = 150
y   = Rstock[1:n1,j]
iv  = CovQ[j,j,1:n1]

# ---------------------------------------------------------------------------- #
#                                Log-Likelihood                                #
# ---------------------------------------------------------------------------- #

function LogLikeGARCH(Θ::AbstractVector{Z};y=data,xreg=0,est=1) where Z

    #parameters variance
    ω0     = Θ[1]
    α0     = Θ[2]
    β0     = Θ[3]
    #initial values
    T     = length(y)
    μ     = mean(y)
    h     = zeros(Z,size(y)) ; h[1] = var(y)
    ϵ     = zeros(Z,size(y)) #; ϵ[1] = std(y)
    LLike = zeros(Z,size(y))

    # ---------------------------------------------------------------------------- #
    #                                transformation                                #
    # ---------------------------------------------------------------------------- #
    UB = 0.9998
    sigmoid(x;ν=0.75) = exp(x)^ν/(1+exp(x)^ν)
    power(x) = x^2

    #par omega
    ω = min(exp(ω0),ω0) #exp(ω0)#

    #par alpha
    β = sigmoid(β0)*UB#max(sigmoid(β0)*UB,β0)#
    
    #par beta
    UB2 = UB - β
    α = sigmoid(α0)*UB2


    # ---------------------------------------------------------------------------- #

    if xreg != 0
        @assert length(xreg)==length(y) "IV no tiene el mismo largo que los retornos >:("
        iv = copy(xreg)
        ζ0  = Θ[4]

        # transformation par zeta 
        # ζ = exp(ζ0)
        # UB3 = UB2 - β
        UB3 = UB2 - α
        ζ = sigmoid(ζ0)*UB3
    else
        iv = zeros(T)
        ζ0  = 0
        ζ  = 0
    end

    for t= 2:T        
        #variance equation ( h = σ² )
        ϵ[t] = y[t] - μ
        h[t] = ω + β*h[t-1]  + α*ϵ[t-1]^2   + ζ*iv[t-1]
        # h[t] = ω0 + β*h[t-1]  + α*ϵ[t-1]^2   + ζ*iv[t-1]
        # h[t] = ω0 + α0*ϵ[t-1]^2  + β0*h[t-1] + ζ0*iv[t-1]

        #likelihood
        LLike[t] = (1/2)*(log(2*pi) + log(h[t]) + (ϵ[t]^2)/(h[t])  )
    end
    sum_LL = sum(LLike)
    #println(Θ)
    if est==1
        return sum_LL
    else 
        xreg == 0 ? params=[ω ; β ; α ] : params=[ω ; β ; α ; ζ]
        return sum_LL, h, params
    end
end

function K_GARCHs(;add_iv= 0,tol=1e-12)
    # add_iv= 1 #incluir IV ? 
    # tol=1e-12
    LL=zeros(k)
    if add_iv == 1
        PAR = zeros(k,4)
        lower = [tol,tol,tol,-Inf]
        upper = [1-tol,1-tol,1-tol,Inf]
    elseif add_iv == 0
        PAR = zeros(k,3)
        lower = [tol,tol,tol]
        upper = [1-tol,1-tol,1-tol]
    end


    D   = zeros(k,k,n1)
    for i in 1:k
        println("GARCH for asset N°",i)

        if add_iv == 1
            xreg = (iv[i,i,:])
            par_init = [0.001, 0.01 , 0.6, 0.3]
        elseif add_iv == 0
            xreg = 0
            par_init = [0.001, 0.01 , 0.95]
        end

        ll_aux(x) = LogLikeGARCH(x;y=y[:,i],xreg=xreg,est=1)
        # opt = optimize(ll_aux, lower, upper, par_init, Fminbox(), #LBFGS(),#NelderMead(),#GradientDescent()
        #                 Optim.Options(show_trace=false,
        #                                 iterations=10_000,
        #                                 # show_every=1,
        #                                 f_tol=1e-12))
        
        opt = optimize(ll_aux, par_init, #LBFGS(),#NelderMead(),#GradientDescent()
                    Optim.Options(show_trace=false,
                                    iterations=10_000,
                                    # show_every=1,
                                    g_tol=1e-8,
                                    f_tol=1e-8))
        
        par= Optim.minimizer(opt) #Estimates

        if add_iv == 1
            ll , h, params=LogLikeGARCH(par;y=y[:,i],xreg=iv[i,i,:],est=0)#con IV
        else
            ll , h, params=LogLikeGARCH(par;y=y[:,i],est=0)#sin IV
        end

        PAR[i,:] = params
        D[i,i,:] = sqrt.(h) # h ≡ σ^2; D = diag(σ) 
        LL[i]  = ll
    end
    return LL,PAR, D
end

LL   , PAR   , D       = K_GARCHs(;add_iv= 0,tol=1e-12)
LL_iv, PAR_iv, D_iv    = K_GARCHs(;add_iv= 1,tol=1e-12)

#comparacion de parametros con libreria
i=19
println("omega , beta , alpha = ",PAR[i,:])
fit_aux=fit(GARCH{1, 1},y[:,i].-mean(y[:,i]))
println(fit_aux) # println(fit_aux.spec.coefs)
LL[i]
# PAR
# PAR_iv


#grafico
i=23
plot(y[:,i])
plot!(D[i,i,:],legend=:none, color=:red)
plot!(D_iv[i,i,:],legend=:none, color=:purple)
plot!(sqrt.(iv[i,i,:]), color=:blue)

# ---------------------------------------------------------------------------- #
#                               DCC: covariances                               #
# ---------------------------------------------------------------------------- #

# 2-stage
# Q(t) = R.*scale + a*e(t-1)'*e(t-1) +  b*Q(t-1)

# where v(t,:) = e(t,;).*(e(t,:)<0) and s = sqrt((1-sum(a)-sum(b)-gScale*sum(g))) and scale = s*s'



function LogLikeDCC( par::AbstractVector{Z} ; y = copy(y), D = copy(D) ,est=1) where Z # ::AbstractVector{Z} where Z
    # Retornos ajustados
    T,k   = size(y)
    #k*(k-1)/2

    # Matrices Q y R
    ε     = zeros(Z,T,k)#*NaN # zeros(T,k) #
    Q     = zeros(Z,k,k,T)#*NaN # zeros(k,k,T) #
    R     = zeros(Z,k,k,T)#*NaN # zeros(k,k,T) #
    LLike = zeros(Z,1,T-1)#*NaN # zeros(T) #

    
    for t = 1:T
         ε[t,:] = (y[t,:] .- mean(y[t,:]))./diag((D[:,:,t])) # D = { √h } = { σ }
    end  # ε= y / σ
    # ε = ε .- mean(ε)

    #parameters
    # par=copy(par0)
    a0        = par[1] 
    b0        = par[2]

    # ---------------------------------------------------------------------------- #
    #                                transformation                                #
    # ---------------------------------------------------------------------------- #
    UB = 0.9998
    sigmoid(x;ν=0.5) = exp(x)^ν/(1+exp(x)^ν)

    #par b
    b        = sigmoid(par[2])*UB

    #par a
    UB2      = UB - b
    a        = sigmoid(par[1])*UB2
    # ---------------------------------------------------------------------------- #

    Q_ss      = cor(ε)
    scale     = (1-a-b)# scale*scale'
    # intercept = Symmetric(Q_ss).*(scale*scale')
    intercept = Symmetric(Q_ss).*(scale)


    #initial values
    likconst = k*2*log(π)
    Q[:,:,1] = cov(ε)#*0.95
    R[:,:,1] = cor(ε)#*0.95
    #Update: @fastmath @inbounds Threads.@threads #t=2
    for t= 1:1:T-1
        Qt = Q[:,:,t]
        εtεtT = ε[t,:]*ε[t,:]'

        # Qt1= intercept  + b*Qt   + a*εtεtT   #DCC
        Qt1= Q_ss  + b*(Qt-Q_ss)   + a*(εtεtT-Q_ss)   #DCC

        q = sqrt.(diag(abs.(Qt1)))
        Rt1 = Qt1 ./ (q*q')  #correlation standarization
        Rt1[diagind(Rt1)] .=1
        Rt1 = Symmetric(Rt1)
        # ---------------------------------------------------------------------------- #
        #                              MULTIVARIATE NORMAL                             #
        # ---------------------------------------------------------------------------- #

        logdetR    = log(det(Rt1))
        invR_εtεtT = sum(diag(inv(Rt1)*εtεtT ))
        # invR_εtεtT = sum(diag((Rt1\I)*εtεtT ))

        LLike[t] = (1/2)*( likconst + logdetR  + invR_εtεtT) #2*logdetD

        Q[:,:,t+1] = Qt1
        R[:,:,t+1] = Rt1

        # ---------------------------------------------------------------------------- #
        #                            MULTIVARIATE t-student                            #
        # ---------------------------------------------------------------------------- #
    end
    Sum_LLike = sum(LLike) #optimize() -> minimization
    if est ==1
        if isnan(Sum_LLike) || ~isreal(Sum_LLike) || isinf(Sum_LLike)
            Sum_LLike = 1e7;
        end
        return Sum_LLike
    else
        H=zeros(k,k,T)
        for tt in 1:n1
            H[:,:,tt] = D[:,:,tt]*R[:,:,tt]*D[:,:,tt]
        end
        params= [a ; b]
        return Sum_LLike, H, R, Q , Q_ss, params
    end
end

function DCC_GARCH(D)
    #paràmetros iniciales
    a0   = 0.2
    b0   = 0.7
    par0 = vec([a0 b0])
    # @time LogLikeDCC(par0)

    # GRADIENTE: para calcular errores (eventualmente)
    # function Grad!(f,x)
    #     f0=f(x)
    #     if length(f(x))==1
    #         g=ForwardDiff.gradient(f,x)
    #         return g
    #     else
    #         T = length(f(x))
    #         k = length(x)
    #         g = zeros(T,k)
    #         Threads.@threads for i in 1:T
    #             f2(x)   = f(x)[i]
    #             @inbounds g[i,:] = ForwardDiff.gradient(f2,x)
    #         end
    #         return g
    #     end
    # end
    # Grad!(x -> LogLikeDCC( x ), par0)

    LL_DCC = TwiceDifferentiable(x -> LogLikeDCC( x ; y = copy(y), D = copy(D) ), par0, autodiff = :forward)


    @time opt2 = optimize(LL_DCC,     par0, #LBFGS(), #GradientDescent(), #LBFGS(), #NelderMead(),#
                            Optim.Options(f_tol        = 1e-12,
                                        g_tol          = 1e-12,
                                        extended_trace = false,
                                        iterations     = 600,
                                        show_every     = 1,
                                        store_trace    = false,
                                        show_trace     = true)) #optimizacion

    parDCC = Optim.minimizer(opt2) #Estimates

    LL, H, R, Q , Q_ss, params = LogLikeDCC(parDCC;est=0)
    return LL, H, R, Q , Q_ss, params
end

LL, H, R, Q , Q_ss, params = DCC_GARCH(D)
LL_iv, H_iv, R_iv, Q_iv , Q_ss_iv, params_iv = DCC_GARCH(D_iv)

println(params_iv)
println(params)


# ---------------------------------------------------------------------------- #
#                       Comprobaciones Gráficas: Plotting                      #
# ---------------------------------------------------------------------------- #
i=24
plot(y[:,i]) # Retornos
plot!(D[i,i,:]) # GARCH con IV ajustado
plot!(D_iv[i,i,:]) # GARCH con IV ajustado
plot!(iv[i,i,:]) # IV ajustado
# plot!(sqrt.(Q[i,i,1:n1]))
# plot!(R[i,ii,1:n1])

ii=5
jj=2
plot(H[ii,jj,1:n1])
plot!(H_iv[ii,jj,1:n1])








x
# ---------------------------------------------------------------------------- #
#                                     OTROS                                    #
# ---------------------------------------------------------------------------- #




# ---------------------------------------------------------------------------- #
#                             OPTION 2: Grid search                            #
# ---------------------------------------------------------------------------- #


# # @fastmath @inbounds 
# # begin #Threads.@threads 
# Threads.nthreads()



# # function MaxLL_MNormal(tol,step)
# tol     = 1e-12
# step    = 0.02
# λ1_grid = collect(0.1:step:0.6) #(1-tol)
# λ2_grid = collect(0.3:step:0.5) #(1-tol)


# λ_grid2D=[]
# for (i,λ1) in enumerate(λ1_grid), (j,λ2) in enumerate(λ2_grid)
#     if (λ1 + λ2) ≤ 1
#         λ_grid2D = push!(λ_grid2D, [λ1 ; λ2] )
#     else
#         nothing
#     end
# end
# N_res =length(λ_grid2D)

# sum_llikes=zeros(N_res)

# Threads.@threads for i in 1:N_res
#     println("solving point N°",i)
#     par = vec(λ_grid2D[i])
#     sum_llikes[i] = LogLikeDCC( par ; y = copy(y), D = copy(D) )
# end

# #     return 
# # end

# plot(sum_llikes)

# ind_max=findmin(sum_llikes)[2]
# println(λ_grid2D[ind_max])



# LogLikeDCC( vec(λ_grid2D[15]) ;y = copy(y), D = copy(D) )



# ---------------------------------------------------------------------------- #
#                              OPCION 3: LIBRERIA                              #
# ---------------------------------------------------------------------------- #

# y ./ diag(D)
# ε = zeros(T,k)
# for t = 1:T
#     ε[t,:] = y[t,:]./diag((D[:,:,t])) # D = { √h } = { σ }
# end  # ε= y / σ

# # example fit(DCC, DOW29)

# fit(DCC,ε)
