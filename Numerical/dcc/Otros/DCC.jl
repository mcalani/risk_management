#cd();cd("Desktop/Julia");push!(LOAD_PATH,".") ;En Julia terminal
using Optim, JuMP, CSV, LinearAlgebra, Statistics
using CSV
using DataFrames
using NLSolversBase
using Plots, StatsPlots


# ---------------------------------------------------------------------------- #
#                                   READ DATA                                  #
# ---------------------------------------------------------------------------- #
cd(dirname(@__FILE__))

Rstock = Matrix(DataFrame(CSV.File("Rstockjl.csv")))[:,2:end]
CovQ = Matrix(DataFrame(CSV.File("CovQjl.csv")))[:,2:end]

CovQ = reshape(CovQ, (24,24,3269)) #Implied Covariances
Rstock= convert(Matrix, Rstock) #Returns

#Setting dimensions
d=vec(1:24)
k=24
j=rand(d,k)

#in sample 1500
y   = Rstock[1:1500,j]
rv  = CovQ[j,j,1:1500]


# ---------------------------------------------------------------------------- #
#                                 DCC: 1st STEP                                #
# ---------------------------------------------------------------------------- #

# 1st-step: Obtaining D=√h
function LLik1(Θ;y=y,estimate=1)
    T,k = size(y)
    T   = Int16(T)
    k   = Int8(k)
    Ik  = Matrix(I, k, k) #Matriz identidad kxk
    H   = ones(k,k,T)
    D   = ones(k,k,T)
    for i in 1:T D[:,:,i]= Ik end
    for i in 1:T H[:,:,i]= Ik end
    #parameters Float32
    ω=Θ[1]
    α=Θ[2]
    β=Θ[3]
    #initial values
    H[:,:,1]=diagm(diag(cov(y))) #matriz de covarianzas
    D[:,:,1]=diagm(diag(sqrt.(cov(y))))#diagonal con desv estandar

    LLike1=0.0
    #update
    for t= 2:T
        H[:,:,t]= ω^2*Ik + α^2*diagm(diag(y[t-1,:]*y[t-1,:]')) + β^2*H[:,:,t-1] # k GARCH(1,1) univariados
        # H[:,:,t]= ω*Ik + α*diagm(diag(y[t-1,:]*y[t-1,:]')) + β*H[:,:,t-1] # k GARCH(1,1) univariados

        D[:,:,t]= sqrt.(H[:,:,t]) # matriz diagonal con k desv estandar
        LLike1 += .5*(2*logdet(D[:,:,t])    +    y[T,:]'*inv(D[:,:,t]^2)*y[t,:])
    end
    if estimate==1
        return LLike1
    else
        return LLike1, D
    end
end

#paràmetros iniciales
w=0.1; a=0.1; b=0.8
par0 = vec([w a b])

#Optimizacion
LL1   = TwiceDifferentiable(x -> LLik1(x;y=y,estimate=1), par0, autodiff = :forward)

@time opt = optimize(LLik1,par0, #GradientDescent(), #LBFGS(), #NelderMead(),
                Optim.Options(g_tol = 1e-12,
                        #iterations = 10,
                        show_every  = 5,
                        store_trace = false,
                        show_trace = true)) #optimizacion
par = Optim.minimizer(opt) #Estimates
println(par)
# println(par.^2)


#check it out: in-sample results
~,D=LLik1(par;estimate=0)
size(D)
i=5
plot(y[:,i])
plot!(D[i,i,:])


#---------------- 2nd-step: DCC -----------------
# Retornos ajustados
ε=zeros(t,k)
for T = 1:t ε[T,:]= y[T,:]./diag(D[:,:,T]) end  #ε=D^-1*y
# Matrices Q y R
Q=zeros(k,k,t)
R=zeros(k,k,t)
# ML Function
function LLik2(Θ2)
    #parameters
    a=Θ2[1]
    b=Θ2[2]
    #initial values
    Q[:,:,1]=cov(ε)
    R[:,:,1]=cor(ε)
    LLike2=0.0
    #Update:
    for T= 2:t
        # Q[:,:,T]=cov(ε)*(1-a^2-b^2)+a^2*ε[T-1,:]ε[T-1,:]'+b^2*Q[:,:,T-1] #DCC
        Q[:,:,T]=cov(ε)*(1-a-b)+a*ε[T-1,:]ε[T-1,:]'+b*Q[:,:,T-1] #DCC

        R[:,:,T]=inv(sqrt.(diagm(diag(abs.(Q[:,:,T])))))*abs.(Q[:,:,T])*inv(sqrt.(diagm(diag(abs.(Q[:,:,T]))))) #cov2corr
        LLike2 += .5*(log(abs.(det(R[:,:,T])))+ ε[T,:]'*inv(abs.(R[:,:,T]))*ε[T,:])
    end
    return LLike2
end

#paràmetros iniciales
a0= 0.1
b0= 0.7
par0 = Vector{Float64}(vec([a0 b0]))

#Optimizacion
@time opt2 = optimize(LLik2,par0,GradientDescent(),
                        Optim.Options(g_tol = 1e-12,
                                #iterations = 10,
                                show_every  = 10,
                                store_trace = false,
                                show_trace = true)) #optimizacion
parDCC = Optim.minimizer(opt2) #Estimates
println(parDCC)

#Comprobaciones Gráficas: Plotting
i=2
plot(Q[i,i,1:20])

plot(ε[:,i]) #residuos e= y/√D
plot(Q[i,i,:])
plot(R[i,1,:])

plot(y[:,i]) # Retornos
plot!(D[i,i,:]) # GARCH ajustado
