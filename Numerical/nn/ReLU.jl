# ---------------------------------------------------------------------------- #
#                                 ReLU function                                #
# ---------------------------------------------------------------------------- #

using LinearAlgebra    #norm function (compute norm of a matrix)
using Plots            #plots
using Optim
using Random         #seed
# using Parameters     
# using LaTeXStrings     #to write LaTeX in legends
# using Statistics       #kernel density
# using JLD              #to save and load results
# using Distributions  
# using NLsolve          
# using Roots            

# ---------------------------------------------------------------------------- #
#                            Function to aproximate                            #
# ---------------------------------------------------------------------------- #
#Example1:
# f(x) = x.^3 .+ x.^2 .- x -1
#Example2: 
f(x) = x.^3 
#Example3: 
# f(x) = sin.(x)

x_min = -5
x_max = 5
x_space = x_min:0.01:x_max


#generate random points from the space defined (sample)
Random.seed!(3)

draw = 20
x_grid=sort(unique(rand(x_space,draw)))

plot(x_space, f.(x_space),label="f(x)",legend=:none)
plot!(x_grid, f.(x_grid),markershape = :circle,label="Sample ($draw)",legend=:topleft)


# ---------------------------------------------------------------------------- #
#                              Activation function                             #
# ---------------------------------------------------------------------------- #
#ReLU

ReLU(x) = max(0,x)

function nnet(par)
    C= zeros(9) ; B = zeros(9)
    C=par[1:length(C)]
    B=par[length(C)+1:length(C)*2]

    app= ReLU.(C[1] .+ B[1].*x_grid) .-
            ReLU.(C[2] .+ B[2].*x_grid) .-
                ReLU.(C[3] .+ B[3].*x_grid) .-
                    ReLU.(C[4] .+ B[4].*x_grid) .-
                        ReLU.(C[5] .+ B[5].*x_grid).-
                            ReLU.(C[6] .+ B[6].*x_grid).-
                                ReLU.(C[7] .+ B[7].*x_grid).-
                                    ReLU.(C[8] .+ B[8].*x_grid).+
                                        ReLU.(C[9] .+ B[9].*x_grid)

    sum((app.-f.(x_grid)).^2)
end

# ---------------------------------------------------------------------------- #
#                                    Fitting                                   #
# ---------------------------------------------------------------------------- #
Random.seed!(3)
par_init=randn(18)
res=optimize(nnet,par_init,Optim.Options(g_tol = 1e-11,f_tol=1e-11,x_tol=1e-11,
                                        iterations = 25_000,
                                        store_trace = true,
                                        show_every = 150,
                                        show_trace = true))

C_opt=Optim.minimizer(res)[1:9]
B_opt=Optim.minimizer(res)[10:18]



approx = ReLU.(C_opt[1] .+ B_opt[1].*x_space) .-
            ReLU.(C_opt[2] .+ B_opt[2].*x_space) .-
                ReLU.(C_opt[3] .+ B_opt[3].*x_space) .-
                    ReLU.(C_opt[4] .+ B_opt[4].*x_space) .-
                        ReLU.(C_opt[5] .+ B_opt[5].*x_space).-
                            ReLU.(C_opt[6] .+ B_opt[6].*x_space).-
                                ReLU.(C_opt[7] .+ B_opt[7].*x_space).-
                                    ReLU.(C_opt[8] .+ B_opt[8].*x_space).+
                                        ReLU.(C_opt[9] .+ B_opt[9].*x_space)




plot(x_space, f.(x_space),linewidth=2)
plot!(x_space, approx,linewidth=2,label="Prediction (sample=$draw)",legend=:topleft)
plot!(x_grid, f.(x_grid),linealpha=0,markershape = :circle,label="Sample points",legend=:topleft)

# ---------------------------------------------------------------------------- #
#                                     Error                                    #
# ---------------------------------------------------------------------------- #
trace = Optim.trace(res)

trace_iter = []
trace_err  = []
trace_time = []
for i in 1:length(trace)
    append!(trace_err, parse(Float64, split(string(trace[i]))[2]))
    append!(trace_iter, parse(Float64, split(string(trace[i]))[1]))
    append!(trace_time, parse(Float64, split(string(trace[i]))[end]))
end
p1=plot(log10.(trace_time), log10.(trace_err/trace_err[end]),xlabel="log(Time)",ylabel="log(Error)");
p2=plot(log10.(trace_iter), log10.(trace_err/trace_err[end]),xlabel="log(Iterations)",ylabel="log(Error)");

plot(p1,p2)