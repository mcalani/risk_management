#-------------------------------------------------------------------------------
#   Huggett 1993 : Discrete states and continuous control
#-------------------------------------------------------------------------------
using Parameters, Plots, QuantEcon
using LinearAlgebra, Statistics
using Optim
using Interpolations

# Each agent wants max E{ sum_{t=0}^{Inf} β^t u(c_{t}) }
# st: c+ a'q ≤ a + e
# where:
# c: consumption
# a: credit balance
# e: endowment
# q: prince of next-credit balance a'

##################################### Parametros ###############################
#Household parameters
β=0.9;
γ=1.5;
function u(c)
    if c > 0.0
    γ==1 ? log(c) : (c^(1-γ) -1)/(1-γ)
    else
    0.0
    end
end
#endowment
e=[1 ; .1];
e_h=e[1];
e_l=e[2];
#Transition between states
π=[0.925 0.075; 0.5 0.5];
#grid on assets (underbar& overbar)
a̅= 4;
a̲= 1e-10;
n_a=200;
a_grid= range(a̲,a̅, length=n_a);
#price of a'
qh=1.2;
ql=0.8;
q=mean([qh ql])
#funcion auxiliar:
#Find the position of the nearest element of v to c (Float)
function findnearest(v::AbstractArray, c::Real)
    findmin(abs.(v .- c ))[2]
end

# V(s)=max_{a'} {U(s-a) + β* Σ V(s')Q(s,a,s')}
# ct= at + e - at+1

################################# Bellman operator #############################
function T(v01,v02,q,a_grid ; compute_policy=false)
    v1_intp=LinearInterpolation(a_grid,v01)
    v2_intp=LinearInterpolation(a_grid,v02)
    Tv1=zeros(n_a)
    Tv2=zeros(n_a)
    Tv= [Tv1 ; Tv2]
    σ1=zeros(n_a)
    σ2=zeros(n_a)
    σ= [σ1, σ2]
    #la maximización se hace desde i, "j" es continuo
    for i in 1:n_a
        obj1(at_1) = u(a_grid[i]+e_h-at_1) +β*(π[1,1]*v1_intp((at_1))
                                              +π[1,2]*v2_intp(at_1))

        obj2(at_1) = u(a_grid[i]+e_l-at_1) +β*(π[2,1]*v1_intp(at_1)
                                              +π[2,2]*v2_intp(at_1))

        res1=maximize(obj1, a̲,(a_grid[i]+e_h))
        res2=maximize(obj2, a̲,(a_grid[i]+e_l))
        Tv1[i]=Optim.maximum(res1)
        Tv2[i]=Optim.maximum(res2)
        if compute_policy
            σ1[i]=Optim.maximizer(res1)
            σ2[i]=Optim.maximizer(res2)
        end
    end
    Tv=[Tv1 ; Tv2]
    if compute_policy
        σ=[σ1 ; σ2]
        return (Tv , σ)
    end
    return Tv
end
v01=log.(a_grid);
v02=log.(a_grid);
q=1.0;
w=T(v01,v02,q,a_grid;compute_policy=false)

################################# Value Function Iteration #####################
function VFI(v01,v02,q,a_grid,T)
    iter=1
    N_iter=500
    conv=100
    tol=1e-5
    v_init=[v01 ; v02]
    V=zeros(n_a)
    σ=zeros(n_a)
    @time while iter < N_iter && conv > tol
        #Bellman operator
        v_next,σ=T(v_init[1:200],v_init[201:400],q,a_grid;compute_policy=true)
        #converce criterion
        conv  = norm(v_next - v_init)
        iter +=1
        println("VF convergence: $conv : Iteration $iter")
        #loop
        v_init[:] = v_next[:]
        V=v_next
    end
    return V , σ
end

V,σ=VFI(v01,v02,q,a_grid,T);
σ_h,σ_l = σ[1:200] , σ[201:400];

################################# Distribution Iteration #######################
collect(a_grid)
[a_grid σ_h]
function stationary_dist(σ_h,σ_l,π)
    conv=100
    tol=1e-5
    pos_h=zeros(n_a,n_a)
    pos_l=zeros(n_a,n_a)
    for i in 1:n_a
        pos_h[i, findnearest(a_grid,σ_h[i])] =1 # 1{a : a'=g(s,a)}
        pos_l[i, findnearest(a_grid,σ_l[i])] =1
    end
    P= [π[1,1]*pos_h  π[1,2]*pos_l
        π[2,1]*pos_h  π[2,2]*pos_l]
    ψ= (ones(2*n_a)/(2*n_a))' #equally spaced distribution
    while conv > tol
        #global conv,tol
        ψ1     = ψ * P
        conv   = norm(ψ1 - ψ)
        println("Distribution convergence: $conv")
        ψ[:,:] =ψ1[:,:]
    end
    return ψ
end

ψ=stationary_dist(σ_h,σ_l,π)
################################### Iteration of price q #######################
#Convergence of price q
function price_q(ψ, qh0, ql0)
   q              = mean([qh0 ql0])
   qh             = Float64[]
   ql             = Float64[]
   ExDD           = ψ*[a_grid ; a_grid]
   ExcessDemand   = push!(Float64[],ExDD)
   if ExDD < 0
       qh = mean([qh0 ql0]) ; ql = ql0
   else
       qh = qh0 ; ql = mean([qh0 ql0])
   end
   return qh, ql
end
qh
ql
price_q(ψ,qh,ql)

################################# Huggett model ################################

v01=log.(a_grid);
v02=log.(a_grid);
qh0=1.2
ql0=0.8

function Huggett_model(v01, v02, qh0, ql0, a_grid)
    q      = mean([qh0 ql0])
    Q      = push!(Float64[], q)
    conv   = 100
    tol    = 1e-5
    iter   = 2
    N_iter = 500
 @time while conv > tol && iter < N_iter
        V, σ       = VFI(v01,v02,q,a_grid,T)
        σ_h, σ_l   = σ[1:200] , σ[201:400]
        ψ          = stationary_dist(σ_h,σ_l,π)
        qh, ql     = price_q(ψ,qh0,ql0)
        q          = mean([qh ql])
        Q          = push!(Q,q)
        conv       = abs.(Q[iter] - Q[iter-1])
        iter      +=1
        println("Convergence of price: $conv")
       end
      return V, σ, Q, ψ
end


############################# Run the model ####################################
qh0=1.2;
ql0=0.8;
v01=log.(a_grid);
v02=log.(a_grid);

V,σ,Q,ψ=Huggett_model(v01,v02, qh0, ql0, a_grid)

plot(V[1:200])
plot!(V[201:400])

plot(σ[1:200])
plot!(σ[201:400])

Ψ_h = cumsum(ψ[1:200]);
Ψ_l = cumsum(ψ[201:400]);

plot(a_grid,Ψ_h)
plot!(a_grid,Ψ_l)
