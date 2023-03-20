#------------------------------------------------------------------------------#
#         Deterministic growth model
#------------------------------------------------------------------------------#
using LinearAlgebra, Statistics
using Parameters, Plots, QuantEcon
using Interpolations

#Initial Values
σ=2;
β=0.95;
δ=0.1;
α=0.33;

#Grid
kstar= ( α ./ ((1/β)-(1-δ)) ) .^ (1/(1-α));
kmin= 0.25*kstar;
kmax=1.75*kstar;
kgrid= 120;
N=kgrid;
kmat=Array(range(kmin,kmax,length= kgrid));

#Preallocation
v1=zeros(kgrid,1);  # maximo de la funcion de valor t
v0=zeros(kgrid,1);  # maximo de la funcion de valor t-1
c=zeros(kgrid,1);   # consumo
k=zeros(kgrid,1);   # capital
k0=zeros(kgrid,1);
K=zeros(kgrid,1);   # capital optimo
C=zeros(kgrid,1);   # consumo optimo
val=ones(kgrid,1)* (-1); # value function
ind= CartesianIndex();
#Loop parameters
tol= 100;
iter=0;
maxit=400;

while  tol> 1e-5 && iter< maxit
    global tol, iter, v0, v1, val,c, k,C,K, maxit, ind

    for i in 1:N, j in 1:N
        k0=kmat[i,1]
        c[j,1]= k0^α -kmat[j,1]+ (1-δ)*k0;

        c[j,1] <=0 ?
        val[j,1]= -999_999 :
        val[j,1]= (1/(1-σ))*(c[j,1] ^(1-σ)-1)+β*v0[j,1];

    v1[i,1],ind = findmax(val)
    K[i,1]= kmat[ind] #encontrar el k que corresponde a cada max val
    C[i,1]=kmat[i,1].^α+(1-δ).*kmat[i,1]-K[i,1]; #consumo optimo
    end
        println("convergence:$tol ; iteration:$iter")

    tol= norm(v1-v0)
    iter+=1;
    v0[:,1]=v1[:,1];
end

p1= plot(v1,kmat, label="Value Function")  #Value function
p2= plot(K,kmat, label= "Policy Function K")    #policy function
plot!(kmat,kmat, label="45º line")
p3= plot(C,kmat,label= "Policy Function C")    #policy function

fg=plot(p1,p2,p3,layout=3)

#savefig("/Users/marco/Desktop/det_growth_model.png")

#------------------------------------------------------------------------------#
#        interpolation: example
#------------------------------------------------------------------------------#

f(x)= 2 .* cos.(6x) .+ sin.(14x) .+ 2.5;
c_grid= 0:.1:1;
f_grid= collect(range(0,1, length=150));

Af = LinearInterpolation(c_grid, f(c_grid));
plt=plot(xlim=(0,1), ylim=(0,6))
plot!(plt, f, f_grid, color= :blue,
 lw=2, alpha= 0.8, label= "true function")
plot!(plt, f_grid, Af.(f_grid), color= :green,lw=2,
 alpha=0.8, label="linear approximation")
plot!(plt, f, c_grid, seriestype= :sticks,
 linestyle= :dash, linewidth=2, alpha=0.8)
plot!(plt, legend= :top)

#------------------------------------------------------------------------------#
#  value function iteration with interpolation (work in progress)
#------------------------------------------------------------------------------#

plot(v1,K)
plot(v1)
plot(K)
V_int= LinearInterpolation(K,v1)

while  tol> 1e-5 && iter< maxit
    global tol, iter, v0, v1, val,c, k,C,K, maxit, ind

    for i in 1:N, j in 1:N
        k0=kmat[i,1]
        c[j,1]= kmat[j,1]^α -kmat[j,1]+ (1-δ)*k0;

        c[j,1] <=0 ?
        val[j,1]= -999_999 :
        val[j,1]= (1/(1-σ))*(c[j,1] ^(1-σ)-1)+β*v0[j,1];

    v1[i,1],ind = findmax(val)
    K[i,1]= kmat[ind] #encontrar el k que corresponde a cada max val
    C[i,1]=kmat[i,1].^α+(1-δ).*kmat[i,1]-K[i,1]; #consumo optimo
    end
        println("convergence:$tol ; iteration:$iter")

    tol= norm(v1-v0)
    iter+=1;
    v0[:,1]=v1[:,1];
end
