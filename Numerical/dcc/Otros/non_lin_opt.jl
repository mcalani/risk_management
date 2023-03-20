using JuMP
using Ipopt #GLKP
using Plots, StatsPlots
using DelimitedFiles
using LinearAlgebra
using Statistics

data=readdlm("Example_FEX.txt")
y = data[:,1] .- mean(data[:,1])

function LogLike(Θ::T,y) where {T<:Real} 
    #parameters variance
    α0 = Θ[1]
    α1 = Θ[2]
    β  = Θ[3]
    #initial values
    h     = zeros(T, size(y)) ; h[1] = var(y)
    ϵ     = zeros(T, size(y)) ; ϵ[1] = std(y)
    y_hat = zeros(T, size(y))
    LLike = zeros(T, size(y))
    for t= 2:(length(y))
        #mean equation
        # y_hat[t] = ϕ0 + ϕ1*y[t-1]
        # ϵ[t] = y[t] -  y_hat[t]
        #variance equation ( h = σ² )
        ϵ[t] = y[t]
        h[t] = α0 + α1*ϵ[t-1]^2  + β*h[t-1]

        #likelihood
        LLike[t] = (-1/2)*(log(2*pi) + log(h[t]) + (ϵ[t]^2)/(h[t])  )
    end
    sum_LL = sum(LLike)
    return sum_LL
end

model = Model(Ipopt.Optimizer)

LL(Θ::Real) = LogLike(Θ,y)
register(model, :LLike, 3, LL; autodiff = true)

#ATTEMPT 1
@variable(model, Θ[1:3])
@constraint(model, Θ[1] >= 0)
@constraint(model, Θ[2] >= 0)
@constraint(model, 0 <=  Θ[3] <= 0.9999)
@constraint(model, Θ[2] + Θ[3] <= 0.99999)

@NLobjective( model, Min, LLike(Θ[1], Θ[2], Θ[3]) )

print(model)
optimize!(model)
@show objective_value(model)
@show value(Θ[1])
@show value(Θ[2])
@show value(Θ[3])



# ---------------------------------------------------------------------------- #
#                                    EXAMPLE                                   #
# ---------------------------------------------------------------------------- #


