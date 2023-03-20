
# ---------------------------------------------------------------------------- #
#                             Parameter Uncertainty                            #
# ---------------------------------------------------------------------------- #
using Random
using Statistics
using Distributions
using LinearAlgebra
using Plots

#Utility function
u(W;A=2) = W^(1-A) /(1-A)

#Wealth
r_f = 0.01
W(ω,T,R_T) = (1-ω) .*exp.(r_f.*T) .+ ω .*exp.(r_f.*T .+ R_T) 


#Time structure (years)
time_bounds = [0 ; 70]
grid_points = 1000
dt = (diff(time_bounds)/grid_points)[]

t = collect(range( time_bounds[1], time_bounds[2], step=dt))
T = maximum(t)

# ---------------------------------------------------------------------------- #
#                                    Wiener                                    #
# ---------------------------------------------------------------------------- #
Np  = 20    #n.processes
Nt  = 252*70  #n.periods
dt  = 1/Nt
μ   = 0.02
σ   = 0.3
S0  = 100

#preallocation
S      =  zeros(Nt,Np); S[1,:] .=S0
d_logS =  zeros(Nt,Np)

for i in 2:Nt
    dW          = rand(Normal(0,1),Np) .* sqrt(dt)
    d_logS[i,:] = ((μ - 0.5*σ^2)*dt .+ σ .*dW) #return
    S[i,:]      = S[i-1,:] .* exp.(d_logS[i,:]) #Continuously compounded
    # S[i,:] = S[i-1,:] .* (1 .+ d_logS)   # Discrete composition
end

#option1
ret_d       = diff(log.(S);dims=1)
sigma_hat   = std(ret_d)/sqrt(dt)
mu_hat      = mean(ret_d)/T + sigma_hat^2/2

#option2
sigma_hat   = std(d_logS)/sqrt(dt)
mu_hat      = mean(d_logS)/T + sigma_hat^2/2




println(mu_hat)
println(sigma_hat)


mean(log.(S[end,:])-log.(S[1,:]))

plot(cumsum(d_logS,dims=1),legend=:none)
plot(S,legend=:none)
plot(d_logS)


gr()
using StatsPlots
density(S[end,:])

log(100)