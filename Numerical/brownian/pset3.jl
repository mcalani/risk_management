using Distributions

p = 0.5
Nt = 252
dt = 1/Nt
μ = 0.06
σ = 0.4
S0 = 100
Np = 10000
binom_dist = Binomial()

binomial_zt = sign.(p .- rand(binom_dist, (Np,Nt+1)))
binomial_zt[:,1] .= 0

daily_returns = μ*dt .+ (σ*sqrt(dt) .* binomial_zt)
cum_returns = cumsum(daily_returns, dims=2)
S = S0 .* exp.(cum_returns)

mu_hat, sigma_hat = mean(S), std(S)
println(mu_hat)
println(sigma_hat)