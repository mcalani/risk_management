# using JuMP
# using TimeSeries
# using Dates
# using XLSX
# using NLSolversBase
using Optim
using Plots, StatsPlots
using DelimitedFiles
using LinearAlgebra
using Statistics
using CSV
using DataFrames
using Distributions


include("utilities.jl");

# log(pdf( Normal(0,1.5), 0.3))
# (-1/2)*(log(2*pi) + log(1.5^2) + (0.3)^2/(1.5^2)  )


# ---------------------------------------------------------------------------- #
#                                   LOAD DATA                                  #
# ---------------------------------------------------------------------------- #

data_csv = CSV.File("tc24febrero.csv") |> Tables.matrix
data     = na_omit!(data_csv)
tc       = Float64.(data[:,2])
ret_tc   = diff(log.(tc))
#mean(ret_tc)

# ---------------------------------------------------------------------------- #
#                                    LOGLIKE                                   #
# ---------------------------------------------------------------------------- #

function LogLike(Θ,y)
    #parameters mean
    # ϕ0=Θ[1]
    # ϕ1=Θ[2]
    #parameters variance
    ω     = Θ[1]
    α     = Θ[2]
    β     = Θ[3]
    #initial values
    μ     = mean(y)
    h     = zeros(size(y)) ; h[1] = var(y)
    ϵ     = zeros(size(y)) #; ϵ[1] = std(y)
    y_hat = zeros(size(y))
    LLike = zeros(size(y))
    for t= 2:(length(y))
        #mean equation
        # y_hat[t] = ϕ0 + ϕ1*y[t-1]
        # ϵ[t] = y[t] -  y_hat[t]
        #variance equation ( h = σ² )
        ϵ[t] = y[t] - μ
        h[t] = ω + α*ϵ[t-1]^2  + β*h[t-1]

        #likelihood
        LLike[t] = (-1/2)*(log(2*pi) + log(h[t]) + (ϵ[t]^2)/(h[t])  )
    end
    sum_LL = sum(LLike)
    #println(Θ)
    return sum_LL
end

max(ll)
min(-ll)

y = copy(ret_tc)
par_init = [0.001, 0.1 , 0.85]
LogLike(par_init,y)



#Optimizacion: poner restricciones w,a,b>0?
@time opt = optimize(x -> LogLike(x,y),par_init,
                        #NelderMead(),
                        LBFGS(),
                        Optim.Options(show_trace=true,
                                        iterations=10_000,
                                        show_every=1,
                                        f_tol=1e-15))

par = Optim.minimizer(opt) #Estimates
println(par)





#advanced optimization 
LL   = TwiceDifferentiable(x -> LogLike(x,ret_tc), par_init, autodiff = :forward)
res  = optimize(LL, par_init, LBFGS(),Optim.Options(show_trace=true,
                                                iterations=100_000,
                                                show_every=10,
                                                f_tol=1e-8))
opt = res.minimizer









# ---------------------------------------------------------------------------- #
#                                   EXAMPLE 1                                  #
# ---------------------------------------------------------------------------- #
using DelimitedFiles
data=readdlm("Example_FEX.txt")
y = data[:,1]
# plot(y)

#paràmetros iniciales Float64
par_init = [0.3, 0.1 , 0.99]
LogLike(par_init,y)

#Optimizacion: poner restricciones w,a,b>0?
@time opt = optimize(x -> LogLike(x,y),par_init,
                        LBFGS(),
                        Optim.Options(show_trace=true,
                                    iterations=100_000,
                                    show_every=1,
                                    f_tol=1e-20))

par = Optim.minimizer(opt) #Estimates
println(par)

# ---------------------------------------------------------------------------- #
#                                   EXAMPLE 2                                  #
# ---------------------------------------------------------------------------- #
using CSV
using DataFrames

ret_stocks = CSV.File("Rstockjl.csv") |> Tables.matrix
#ret_stocks = DataFrame(CSV.File("Rstockjl.csv"))

y2= ret_stocks[:,4] .- mean(ret_stocks[:,4])
#plot(y2)
mean(y2)

#paràmetros iniciales
par_init = [0.3, 0.1 , 0.99]
LogLike(par_init,y2)

# LL   = TwiceDifferentiable(x -> LogLike(x,y2), par_init, autodiff = :forward)
#Optimizacion: poner restricciones w,a,b>0?
@time opt2 = optimize(x->LogLike(x,y2),par_init, #[1e-5 0.2], [0.5 0.999],
                        LBFGS(),
                        Optim.Options(show_trace=true,
                                    iterations=100_000,
                                    show_every=1,
                                    #h_tol=1e-16,
                                    g_tol=1e-16,
                                    f_tol=1e-16))

par2 = Optim.minimizer(opt2) #Estimates
println(par2)

sum(par2)