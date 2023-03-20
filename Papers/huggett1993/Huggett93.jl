#-------------------------------------------------------------------------------
#   Huggett 1993
#-------------------------------------------------------------------------------
using Parameters, Plots, QuantEcon
using LinearAlgebra, Statistics
using Optim
using Interpolations

#Each agent wants max E{  sum_{t=0}^{Inf} β^t u(c_{t}) }
# st: c+ a'q ≤ a + e
# where:
# c: consumption
# a: credit balance
# e: endowment
# q: prince of next-credit balance a'

#Household parameters
β=0.9;
σ=1.5;
u(c)= σ==1 ? log.(c) : (c^(1-σ) -1)/(1-σ);
#endowment
e=[1 .1];
e_h=e[1];
e_l=e[2];
#Transition between states
π=[0.925 0.075; 0.5 0.5];
#grid on assets (underbar& overbar)
a̅= 4;
a̲= 1e-10;
n_a=200;
a= Array(range(a̲,a̅, length=n_a));
#price of a'
q=1;
#Loop parameters
c1=zeros(n_a);    #consumo
v1=zeros(n_a);
v01=zeros(n_a);   #value function
val_h=zeros(n_a);
c2=zeros(n_a);
v2=zeros(n_a);
v02=zeros(n_a);
val_l=zeros(n_a);
tol=1e-5;
conv1=100;
conv2=100;
conv3=100;
#Excess of demand
ExcessDemand=Float64[];
Q=Float64[];
Q=push!(Q,0.1);
t=2;
qhi=1.2;
qlo=β/2;
Qhi=Float64[];
Qlo=Float64[];

################################################################################
################################################################################

@time while conv3 > tol
    global qhi, qlo, t, Qhi, Qlo
    global q= (qhi+qlo)/2
################################################################################
#Value Function Iteration
    global conv1=100;
   while conv1 > tol
       global iter, maxit, ind_h, ind_l
       global v1_i=Int[];
       global v2_i=Int[];
       for i in 1:n_a
       a_t=a[i]
           for j in 1:n_a
           a_t_1=a[j]
           c1[j]=a_t + e_h - q*a_t_1
           c2[j]=a_t + e_l - q*a_t_1
           c1[j]>0 ?
           val_h[j]= u(c1[j]) + β*(π[1,1] *v01[j] + π[1,2] *v02[j]) :
           val_h[j]= -999_999
           c2[j]>0 ?
           val_l[j]= u(c2[j]) + β*(π[2,1] *v01[j] + π[2,2] *v02[j]) :
           val_l[j]=-999_999
           end
          #optimal vf
          v1[i],ind_h=findmax(val_h)
          v2[i],ind_l=findmax(val_l)
          #indices of controls
          global v1_i=push!(v1_i,ind_h)
          global v2_i=push!(v2_i,ind_l)
          end
          global conv1=norm(v1-v01)
          println("VF Convergence: $conv1")
       v01[:]=v1[:]
       v02[:]=v2[:]
   end
################################################################################
#Distribution Iteration
   pos_h=zeros(n_a,n_a);
   pos_l=zeros(n_a,n_a);
#coordinates pos(state, control)
   for i in 1:n_a
       pos_h[i,v1_i[i]]=1; pos_l[i,v2_i[i]]=1
   end
#transition probabilities from state i to control j, given e ϵ[e_h,e_l]
   P= [π[1,1]*pos_h  π[1,2]*pos_l
       π[2,1]*pos_h  π[2,2]*pos_l];
#Space distribution
   global ψ= (ones(2*n_a)/(2*n_a))';
   global conv2=100;
   while conv2 > tol
       global tol, ψ, conv
       ψ1= ψ * P
       global conv2= maximum(abs.(ψ1-ψ))
       println("Distribution convergence: $conv2")
       ψ[:,:]=ψ1[:,:]
   end
################################################################################
#Convergence of price q
   global ExDD= ψ*[a ; a];
   global ExcessDemand=push!(ExcessDemand,ExDD)
   global Q=push!(Q,q)
   global conv3=abs.(Q[t] - Q[t-1])
         if t>1 && conv3 <tol
         break
         end
         ExDD<0 ?
         qhi = (qhi + qlo) / 2 :
         qlo = (qhi + qlo) / 2 ;
         Qhi  =push!(Qhi,qhi);
         Qlo  =push!(Qlo,qlo);
         t+=1;
   println("Excess of Demand: $conv3")
end
################################################################################
################################################################################

#optimal value functions
plot(v1)
plot!(v2)

#optimal decision rule
σ_h=zeros(n_a);
σ_l=zeros(n_a);
for i in 1:n_a
    σ_h[i]= a[v1_i[i]]
    σ_l[i]= a[v2_i[i]]
end
plot(σ_h)
plot!(σ_l)

#price convergence
plot(Q)
#Excess Demand converge
plot(ExcessDemand)

#distributions
plot(ψ[1:400])
Ψ_h=cumsum(ψ[1:200]);
Ψ_l=cumsum(ψ[201:400]);

plot(a,Ψ_h)
plot!(a,Ψ_l)
