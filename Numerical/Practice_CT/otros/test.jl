using Distributions
using Plots
using Random

# ---------------------------------------------------------------------------- #
#                                Simulate Prices                               #
# ---------------------------------------------------------------------------- #
Random.seed!(3)

#R.W : Y_t = z_1 +...+ z_T , where z_t ~ N(μ,σ^2)
#mean(Y) = μ
#var(Y)  = σ^2

ns= 100  # Number of stocks
T = 252      # Days (252 in a year)
dt= 1/T      # aily step

#R.W moments (in a year)
μ = .15
σ = 0.3 
z=rand(Normal(dt*μ,sqrt(dt)*σ),ns,T)

y=zeros(size(z))
for t in 2:T
    y[:,t] = y[:,t-1] + z[:,t-1]
end 

#check
mean(y[:,T]) # -> μ in the limit ns -> ∞
std(y[:,T])  # -> σ in the limit ns -> ∞

P_0 = 100
P_t = P_0 .* exp.(y) 

plot(P_t',legend=:none,title="Simulated Prices",yaxis="Price",xaxis="Date")


# ---------------------------------------------------------------------------- #
#                                    Returns                                   #
# ---------------------------------------------------------------------------- #

# (Equally weighted) portfolios' profit between two dates t1 and t2 (t ∈ [1, T])

function Return(t0,t1,P)
    println("Return between periods $t0 and $t1")
    
    #period_ret = mean( P_t[:,t1] ./P_t[:,t0] .-1 )
    period_ret = mean(log.(P_t[:,t1] ./P_t[:,1]) )
    an_ret = (  (period_ret + 1)^((T-1)/(t1-t0))  -1   )*100
    println("Annual Return:",round(an_ret,digits=4),"%")
    println("Period return:",round(period_ret,digits=4),"%")
end


#Arbitrary dates
t0= 1
t1= 3
Return(t0,t1,P_t)

#All period (1 year)
Return(1,T,P_t)

### Note that annual returns are -> μ in the limit ns -> ∞