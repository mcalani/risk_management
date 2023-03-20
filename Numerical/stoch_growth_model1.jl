#------------------------------------------------------------------------------#
# QuantEcon Cap.29: Optimal Growth I: The stochastic optimal growth model
#------------------------------------------------------------------------------#
# y_t+1 := f(k_t+1)ξ_t+1, ξ_t+1~ϕ
# constraints: k_t+1+c_t =< y_t
using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version="0.5.0")
using LinearAlgebra, Statistics, Optim, Random
using Parameters, Plots, QuantEcon, Interpolations, NLsolve

#Bellman operator:
function T(w, grid, β, u, f, shocks, Tw=similar(w);compute_policy=false)
    w_func = LinearInterpolation(grid,w) #w(grid)
    #objective for each grid point
    objectives= (c -> u(c) +β*mean(w_func.(f(y-c).* shocks))
    for y in grid_y)
        results= maximize.(objectives, 1e-10, grid_y)
        Tw=Optim.maximum.(results)
        if compute_policy
            σ= Optim.maximizer.(results)
            return Tw, σ
        end
    return Tw
end
#-------------------------------------------------------#
#Example:
#-------------------------------------------------------#

# y=f(k)=(k^α) * exp(μ+σζ) ,    exp(μ+σζ)~ ϕ, ζ~N(0,1)
# u(c)=ln(c)
# k_t+1+c_t =< y_t

#initial values
α=0.4;
β=0.96;
μ=0;
s=0.1;  # s= σ en exp(μ+σζ)
#-------------------------------------------------------#
#True values
c1=log(1-α*β)/(1-β);
c2=(μ +α*log(α*β))/(1-α);
c3= 1/(1-β);
c4= 1/(1-α*β);
#True optimal policy
c_star(y)=(1-α*β)*y;
#True value function
v_star(y)= c1+c2*(c3-c4)+c4*log(y);

#-------------------------------------------------------#

#Utility
u(c)= log(c);
∂u∂c(c)= 1/c;

#deterministic part of production function
f(k)= k^α;
deriv_f(k)= α*k^(α-1);

grid_max=4;
grid_size= 200;
shock_size=250;

grid_y= range(1e-5, grid_max, length=grid_size);
shocks= exp.(μ .+ s*randn(shock_size));

#Policy function
w=T(v_star.(grid_y), grid_y, β, log, k->k^α, shocks); #w-greedy

plt= plot(ylim= (-35,-24))
plot!(plt, grid_y, w, linewidth=2, alpha=0.6, label="T(v_star)")
plot!(plt, v_star, grid_y, linewidth=2, alpha=0.6, label="v_star")

#-------------------------------------------------------#
#Plot
global w= 5*log.(grid_y);
n=35;

plot(xlim=(extrema(grid_y)), ylim=(-50,10))
lb="initial condition"
plt=plot(grid_y,w,color=:green, linewidth=2,alpha=0.8, label=lb)

for i in 1:n
    global w
    w=T(w, grid_y, β, log, k -> k^α, shocks)
    plot!(grid_y,w, color=RGBA(i/n, 0,1-i/n,0.8),
    linewidth=2, alpha=0.6, label="")
end

lb="true value function";
plot!(plt, v_star, grid_y, color= :black, linewidth=2, alpha=0.8,label=lb)
plot!(plt, legend= :bottomright)
