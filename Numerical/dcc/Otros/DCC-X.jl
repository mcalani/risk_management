#cd();cd("Desktop/Julia");push!(LOAD_PATH,".") ;En Julia terminal
using Optim, JuMP, CSV, LinearAlgebra, Statistics
using CSV
using DataFrames
using NLSolversBase
using Plots, StatsPlots

#Loading data
@time Rstock = CSV.read("/Users/marco/Desktop/Julia/Rstockjl.csv", copycols=false)
@time CovQ = CSV.read("/Users/marco/Desktop/Julia/CovQjl.csv", copycols=false)
#typeof(rv)
#typeof(y)

#Dataframe to Arrays (IN-SAMPLE t=2000 aprox)
CovQ=convert(Matrix, CovQ[2:78457]) #deprecated
CovQ = reshape(CovQ, (24,24,3269)) #Implied Covariances
Rstock= convert(Matrix, Rstock[2:25]) #Returns

#Setting dimensions
d=vec(1:24)
k=4
j=rand(d,k)

#in sample 1000
y=Rstock[1:1000,j]
rv=CovQ[j,j,1:1000]

t,k=size(y)
t=Int16(t)
k=Int8(k)

Ik=Matrix{Int8}(I, k, k) #Matriz identidad kxk

H=ones(k,k,t)
D=ones(k,k,t)
for i in 1:t D[:,:,i]= Ik end
for i in 1:t H[:,:,i]= Ik end
############################## DCC-X estimation ################################
#----------------------------  1st-step: Obtaining D=√h ------------------------
function LLik1(Θ1)
    #parameters Float32
    ω=Float32(Θ1[1])
    α=Float32(Θ1[2])
    β=Float32(Θ1[3])
    #initial values
    H[:,:,1]=diagm(diag(cov(y))) #matriz de covarianzas
    D[:,:,1]=diagm(diag(sqrt.(cov(y))))#diagonal con desv estandar
    LLike1=0.0
    #update
    for T= 2:t
        H[:,:,T]= ω^2*Ik + α^2*diagm(diag(y[T-1,:]*y[T-1,:]')) + β^2*H[:,:,T-1] # k GARCH(1,1) univariados
        D[:,:,T]=sqrt.(H[:,:,T]) # matriz diagonal con k desv estandar
        LLike1 += .5*(2*logdet(D[:,:,T])+y[T,:]'*inv(D[:,:,T]^2)*y[T,:])
    end
    return LLike1
end

#paràmetros iniciales Float64
w=0.1
a=0.1
b=0.9
par0 = Vector{Float64}(vec([w a b])) #Float64= 13segs

#Optimizacion: poner restricciones w,a,b>0
@time opt = optimize(LLik1,par0) #optimizacion
par = Optim.minimizer(opt) #Estimates
println(par)

plot(D[1,1,:])

#----------------------  2nd-step: DCC-X ---------------------------------------
#Matriz media de RV escalado por 1/100
mean_rv=zeros(k,k)
for i= 1:k, j=1:k mean_rv[i,j] = mean(rv[i,j,:]) end
mean_rv

# Retornos ajustados
ε=zeros(t,k)
for T = 1:t ε[T,:]= y[T,:]./diag(D[:,:,T]) end  #ε = D^-1*y
plot(y[:,1])
plot!(ε[:,1])

# Matrices Q y R
R=zeros(k,k,t)
Q=zeros(k,k,t)
# ML Function
function LLik2(Θ2)
    #parameters
    a=Θ2[1]
    b=Θ2[2]
    c=Θ2[3]
    #initial values
    Q[:,:,1]=cov(ε)
    R[:,:,1]=cor(ε)
    #R[1,:]=vec(cor(ε))
    LLike2=0.0
    #update
    for T= 2:t
    Q[:,:,T]= (cov(ε)*(1+a^2+b^2)+c^2*cov(ε))+a^2*ε[T-1,:]ε[T-1,:]'+b^2*Q[:,:,T-1] #DCCX
    R[:,:,T]=inv(sqrt.(diagm((diag(Q[:,:,T])))))*(Q[:,:,T])*inv(sqrt.(diagm((diag(Q[:,:,T]))))) #cov2corr
    LLike2 += .5*(log((det(R[:,:,T])))+ε[T,:]'*inv(R[:,:,T])*ε[T,:])
    end
    return LLike2
end

#paràmetros iniciales
a0=0.1
b0=0.8
c0=0.2
par0 = Vector{Float64}(vec([a0 b0 c0]))

#Optimizacion: PONER RESTRICCIONES PARA QUE Q SEA POSTIVA SEMIDEFINIDA
@time opt2 = optimize(LLik2,par0;iterations=10000) #optimizacion
parDCCX = Optim.minimizer(opt2) #Estimates
println(parDCCX)

#Plotting
plot(y[:,2].^2)
plot!(rv[2,2,:])

plot(y[:,1].^2)
plot(Q[1,1,:])
plot(R[1,3,:])

#cholesky decomposition
cov(ε)
CH=cholesky(cov(ε))
CH.L*CH.U

CH.L.*ones(k,k)
F=CH.L*Ik
F*F'

LowerTriangular(cov(ε))*Ik

#---------- Posibles soluciones
# 1.Modelar R y no Q?
# 2.Constrained Optim
# 3. Notar que RV puede explicar Y, y tal vez no a e=Y/√D
