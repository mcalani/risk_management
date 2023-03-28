using Distributions
using Plots

θ = 0.95; # Vary persistence to check behavior
μ = 1;
σ = 0.007;

# Numerical parameters
tb = [0 1];  # Time bounds 
N  = 10000;  # 10000 grid points

# Components simulation
dt = (tb[2]-tb[1])./N;
t  = collect(range(tb[1],tb[2],length=N)); # From t0-->t1 with N points
y  = zeros(N);              # 1xN Matrix of zeros

# Simulated paths of N points

y[1] = μ ;  # Initial y condition (converge to the mean=1)


i=1
for i = 2:N
    dW   = rand(Normal(0,sqrt(dt)),1)[]
    y[i] = y[i-1] + θ*(μ-y[i-1])*dt + σ*dW
end

plot(y)