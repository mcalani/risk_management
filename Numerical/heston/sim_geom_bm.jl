# ---------------------------------------------------------------------------- #
#                      Simulate Geometric Brownian Motion                      #
# ---------------------------------------------------------------------------- #
using Random
using Statistics
using Distributions
using LinearAlgebra
using Plots


# ---------------------------------------------------------------------------- #
#                               Binomial process                               #
# ---------------------------------------------------------------------------- #
p   = 0.5
Np  = 1    #n.processes
Nt  = 25200  #n.periods
dt  = 1/Nt
μ   = 0.06
σ   = 0.4
S0  = 100

binomial_zt = sign.(p .- rand(Binomial(), (Np,Nt+1)))
binomial_zt[:,1] .= 0

daily_returns = μ*dt .+ σ*sqrt(dt) .* binomial_zt
cum_returns = cumsum(daily_returns, dims=2)
S = S0 .* exp.(cum_returns)

plot(S')

mu_hat      = mean(daily_returns)/Nt # - std(daily_returns)^2/2
sigma_hat   = std(daily_returns)/sqrt(dt)
# plot(daily_returns')




println(mu_hat)
println(sigma_hat)


# ---------------------------------------------------------------------------- #
#                                    Wiener                                    #
# ---------------------------------------------------------------------------- #
Np  = 1    #n.processes
Nt  = 2520  #n.periods
dt  = 1/Nt
μ   = 0.1
σ   = 0.4
S0  = 100

S   =   zeros(Nt,Np); S[1,:] .=S0
d_logS =  zeros(Nt,Np)
for i in 2:Nt
    dW = rand(Normal(0,1),Np) .* sqrt(dt) 
    d_logS[1,:] = ((μ - σ^2 /2)*dt .+ σ .*dW)[] #return
    S[i,:] = S[i-1,:] .* exp.(d_logS[i,:]) #Continuously compounded
    # S[i,:] = S[i-1,:] .* (1 .+ d_logS)   # Discrete composition
end


(mean(d_logS) - std(d_logS)^2/2)*dt
std(d_logS)/sqrt(dt)



ret_d       = diff(log.(S);dims=1)
mu_hat      = (mean(ret_d)-std(ret_d)^2/2)*dt
sigma_hat   = std(ret_d)/sqrt(dt)

println(mu_hat)
println(sigma_hat)


plot(S)
plot(ret_d)
plot(diff(S;dims=1))
