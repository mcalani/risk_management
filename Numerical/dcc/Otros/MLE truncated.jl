using Distributions, Optim

# hard-coded data\observations
odr=[0.10,0.20,0.15,0.22,0.15,0.10,0.08,0.09,0.12]
Q_t = quantile.(Normal(0,1), odr)

# return a function that accepts `[mu, sigma]` as parameter
function neglik_tn(Q_t)
    maxx = maximum(Q_t)
    f(μσ) = -sum(logpdf.(Truncated(Normal(μσ[1],μσ[2]), -Inf, maxx), Q_t))
    f
end

neglikfn = neglik_tn(Q_t)

# optimize!
# start searching
@time res = optimize(neglikfn, [mean(Q_t), std(Q_t)]) # 0.8 seconds
@time res = optimize(neglikfn, [mean(Q_t), std(Q_t)]) # 0.000137 seconds

# the \mu and \sigma estimates
Optim.minimizer(res) # [-1.0733250637041452,0.2537450497038758] # or

# use `fieldnames(res)` to see the list of field names that can be referenced via . (dot)
res.minimizer # [-1.0733250637041452,0.2537450497038758]
println(res.minimizer)
