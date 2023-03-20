#cd();cd("Desktop/Julia");push!(LOAD_PATH,".") ;En Julia terminal
using Optim, JuMP, CSV, LinearAlgebra, Statistics
using CSV
using DataFrames
using NLSolversBase
using Plots, StatsPlots

####################     Loading data     ########################
@time Rstock = CSV.read("/Users/marco/Desktop/Julia/Rstockjl.csv", copycols=false)
@time CovQ = CSV.read("/Users/marco/Desktop/Julia/CovQjl.csv", copycols=false)
#typeof(rv)
#typeof(y)

#Dataframe to Arrays
CovQ=convert(Matrix, CovQ[2:78457]) #deprecated
CovQ = reshape(CovQ, (24,24,3269)) #Implied Covariances
Rstock= convert(Matrix, Rstock[2:25]) #Returns

#Setting dimensions
d=vec(1:24)
k=3
j=rand(d,k)

#in sample 1000
y=Rstock[1:1000,j]
rv=CovQ[j,j,1:1000]
#for i = 1:k, j=1:t y[j,i]=y[j,i]-mean(y[:,i]) end #substract mean

t,k=size(y)
t=Int16(t)
k=Int8(k)

Ik=Matrix{Int8}(I, k, k) #Matriz identidad kxk

#Matrices H y D
H=ones(k,k,t)
D=ones(k,k,t)
# Retornos ajustados
ε=zeros(t,k)
# Matrices Q y R
R=zeros(k,k,t)
Q=zeros(k,k,t)
################      DCC estimation    ###############

function LLik(Θ)
    #parameters Float32
    ω=Float32(Θ[1])
    α=Float32(Θ[2])
    β=Float32(Θ[3])
    #parameters
    a=Float64(Θ[4])
    b=Float64(Θ[5])
    c=Float64(Θ[6])
    #initial values
    #initial values
    H[:,:,1]=diagm(diag(cov(y)))
    D[:,:,1]=diagm(diag(sqrt.(cov(y))))#diagonal con desv estandar
    Q[:,:,1]=cov(ε)
    R[:,:,1]=cor(ε)
    LLike1=0.0
    LLike2=0.0
    LLike=LLike1+LLike2
    #update
    for T= 2:t
        #1st step
        H[:,:,T]= ω^2*Ik+α^2*diagm(diag(y[T-1,:]*y[T-1,:]'))+β^2*H[:,:,T-1] # k GARCH(1,1) univariados
        D[:,:,T]=sqrt.(H[:,:,T]) # matriz diagonal con k desv estandar
        ε[T,:]= y[T,:]./diag(D[:,:,T]) #retornos ajustados
        #2nd step
        Q[:,:,T]=cov(ε)*(1-a^2-b^2)+a^2*ε[T-1,:]ε[T-1,:]'+b^2*Q[:,:,T-1] #DCC
        R[:,:,T]=inv(sqrt.(diagm(diag(abs.(Q[:,:,T])))))*abs.(Q[:,:,T])*inv(sqrt.(diagm(diag(abs.(Q[:,:,T]))))) #cov2corr
        LLike1 += .5*(2*logdet(D[:,:,T])+y[T,:]'*inv(D[:,:,T]^2)*y[T,:])
        LLike2 += .5*(log(abs.(det(R[:,:,T])))+ε[T,:]'*inv(R[:,:,T])*ε[T,:])
        LLike=LLike1+LLike2
    end
    return LLike
end

#paràmetros iniciales GARCH
ω0=0.01
α0=0.1
β0= 0.8
#paràmetros DCCX
a0=0.1
b0=0.9
c0=0.2

par0=[ω0,α0,β0, a0, b0, c0]

#Optimizacion: poner restricciones w,a,b>0
@time opt = optimize(LLik,par0; iterations=1000) #optimizacion
parDCC = Optim.minimizer(opt) #Estimates
println(parDCC)

plot(y[:,1])
D[3,3,1:50]
