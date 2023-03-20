#------------------------------------------------------------------------------#
# QuantEcon Cap.29: Optimal Growth I: The stochastic optimal growth model
#------------------------------------------------------------------------------#
using LinearAlgebra, Statistics, Optim, Random
using Parameters, Plots, QuantEcon, Interpolations, NLsolve
Random.seed!(42)

#Bellman operator: with monte carlo shocks
function T(w, grid, β, u, f, shocks, Tw=similar(w);compute_policy=false)
    w_func = LinearInterpolation(grid,w) #w(grid)
    #objective for each grid point
    #the expectation is computed via monte carlo (with mean function)
    objectives= (c -> u(c)+β*mean(w_func.(f(y-c).* shocks))
    for y in grid_y)
        results= maximize.(objectives, 1e-10, grid_y)
        Tw=Optim.maximum.(results)
        if compute_policy
            σ= Optim.maximizer.(results)
            return Tw, σ
        end
    return Tw
end

#initial values
α=0.4;
β=0.96;
μ=0.0;
σ_ϵ=0.1;  # σ en exp(μ+σζ)

#grid for y=f(k)
grid_max=4;      #maximum value in the grid
grid_size= 200;
shock_size=250;
grid_y = range(1e-5, grid_max, length=grid_size);
shocks = exp.(μ .+ σ_ϵ * randn(shock_size));

#initial guess w(y)= 5ln(y) (arbitrary)
w0=5*log.(grid_y);
#Bellman operator and policy function σ_c
w1, σ_c=T(w0, grid_y, β, log, k->k^α, shocks, compute_policy=true);

#-------------------------------------------------------#
#Comparision with true Values
#-------------------------------------------------------#

#True values: from recursive macroeconomics section 3.1.2
c1=log(1-α*β)/(1-β);
c2=(μ +α*log(α*β))/(1-α);
c3= 1/(1-β);
c4= 1/(1-α*β);
#True optimal policy
c_star(y)=(1-α*β)*y;
#True value function
v_star(y)= c1+c2*(c3-c4)+c4*log(y);

w2=T(v_star.(grid_y), grid_y, β, log, k->k^α, shocks);

#we expect the same function (practically)
plot(grid_y,w1, label="Bellman Operator")
plot!(grid_y,w0, label="Initial Value")
plot!(v_star,grid_y, label= "True Value Function")

#convergencia de la funcion de valor
global w= 5*log.(grid_y);
n=40;
for i in 1:n
    # notar que este loop está aplicando el operador de Bellman
    # al value function w hasta que converge al true valor v*
    global w
    w=T(w, grid_y, β, log, k -> k^α, shocks)
    plot!(grid_y,w, color=RGBA(i/n, 0,1-i/n,0.8),
    linewidth=1, alpha=0.6, label="")
end
plot!(legend= :bottomright)

#-------------------------------------------------------#
#      Compute fixed point
#-------------------------------------------------------#

u(c)=log(c);
f(k)=k^α;

#funcion para encontrar punto fijo
function solve_optgrowth(initial_w; tol = 1e-6, max_iter = 500)
Tw = similar(grid_y);
v_star_approx = fixedpoint(w -> T(w, grid_y, β, u, f, shocks, Tw),
initial_w).zero; # gets returned
end

#valor inicial y value function
w0 = 5 * log.(grid_y);
v_star_approx = solve_optgrowth(w0)


plt = plot(ylim = (-35, -24))
plot!(plt, grid_y, v_star_approx, linewidth = 2, alpha = 0.6,
label = "approximate value function")
plot!(plt, v_star, grid_y, linewidth = 2, alpha = 0.6, label = "true￿
↪value function")
plot!(plt, legend = :bottomright)
