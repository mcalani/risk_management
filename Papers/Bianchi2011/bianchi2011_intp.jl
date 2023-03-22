# ---------------------------------------------------------------------------- #
#                         Replication of Bianchi (2011)                        #
# ---------------------------------------------------------------------------- #
using Distributions
using Interpolations #interpolations
using Optim          #optimize
using Roots          #find roots
using ProgressMeter  #to show progress
using LinearAlgebra  #norm function (compute norm of a matrix)
using JLD            #to save and load results
using Plots          #plots
using LaTeXStrings   #to write LaTeX in legends
using QuantEcon      #stationary distributions
using StatsPlots     #kernel density
using Random         #seed
cd(dirname(@__FILE__))
#Load stochastic structure
include("calibration.jl")
include("aux_funcs.jl")

#Parameters
r  = 0.04             #interest rate
σ  = 2                #risk aversion
η  = 1/0.83-1         #elasticity of substitution
ω  = 0.307000693252802             #weight on tradables
β  = 0.906             #discount factor
κ  = 0.3235
n  = 4                #nodes
ny = n*n

#grid for b
b_min  = -1.02
b_max  = -0.4

#OPTION 1
# nb     =  800
# b_grid = collect(range(b_min,b_max,length=nb))

#OPTION 2
Random.seed!(3) ; x=rand(Normal(mean([b_max b_min]),1),500)
b_grid= sort(unique(clamp.(x,b_min,b_max)))
nb    =  length(b_grid)

density(b_grid, title="Grid points distribution")

#grid for y
YN = repeat(yN_grid,1,nb)
YT = repeat(yT_grid,1,nb)
B  = repeat(b_grid,1, ny)'

#initial values (matrices)
bp = repeat(b_grid,1,ny)'
cT = ones(ny,nb)
pN = ones(ny,nb)
μ  = zeros(ny,nb)
Eμ = zeros(ny,nb)

#Marginal utility
mgu(ct,yn)=( ω*ct^(-η) + (1-ω)*yn^(-η) )^(σ/η-1/η-1) * ω * ct^(-η-1)
# ---------------------------------------------------------------------------- #
#                            Competitive Equilibrium                           #
# ---------------------------------------------------------------------------- #
# Technical parameters
updt    = 0.2  ; #Updating rule: Must be slow, important!
n_iter  = 700  ;
counter = 0    ;
tol     = 1e-8 ;
show_every= 10 ;

while counter < n_iter
    global updt,counter,n_iter, μ_in ,Eμ,cT,bp,pN,cT_bind,bp_bind

    #save old values
    bp_old = copy(bp)
    cT_old = copy(cT)
    pN_old = copy(pN)

    #Next binding values
    bp_bind = -κ.*(pN.*YN+YT)
    cT_bind = (1+r).*B .+ YT .- bp_bind

    #marginal utility
    λ = ( ω.*cT.^(-η) .+ (1-ω).*YN.^(-η) ).^(σ/η-1/η-1) .* ω .* cT.^(-η-1)

    #Expected marginal utility
    for i in 1:nb
        for j in 1:ny
            Eμ[j,i] = β*(1+r)*T[j,:]'*interp1(b_grid,λ,bp[j,i])
        end
    end

    #Compute Euler's residual assuming that constraint binds
    λ_bind = ( ω.*cT_bind.^(-η) .+ (1-ω).*YN.^(-η) ).^(σ/η-1/η-1) .* ω .* cT_bind.^(-η-1)
    μ      = 1 .- Eμ./λ_bind

    #system of equations
    bp[μ .≥ tol] = -κ.*(pN.*YN+YT)[μ .≥ tol]
    cT[μ .≥ tol] = ((1+r).*B .+ YT .- bp)[μ .≥ tol]
    pN[μ .≥ tol] = ((1-ω)/ω .* (cT./YN).^(1+η))[μ .≥ tol]

    for i in 1:nb
        for j in 1:ny
            if  μ[j,i] < -tol
                u(cc)   = mgu(cc,YN[j,i]) .- Eμ[j,i] # 1 .- Eμ[j,i]/mgu(cc,YN[j,i])
                cT[j,i] = fzero( u , 0.1 ) ; μ[j,i] = u(cT[j,i])  # (check μ ≥ 0)
                bp[j,i] = (B*(1+r) + YT - cT)[j,i]
                pN[j,i] = ((1-ω)/ω .* (cT./YN).^(1+η))[j,i]
            end
        end
    end
    
    #=check collateral constraint
    bp[ bp .< -κ*(pN.*YN+YT) ] = (-κ.*(pN.*YN+YT))[  bp .< -κ*(pN.*YN+YT) ]
    cT = ((1+r).*B .+ YT .- bp)
    pN = (1-ω)/ω .* (cT./YN).^(1+η) # =#

    #updating rule
    bp =  bp_old + updt.*(bp - bp_old)
    cT =  cT_old + updt.*(cT - cT_old)
    pN =  pN_old + updt.*(pN - pN_old)

    #convergence criterion
    diff=maximum( [ maximum(abs.(bp_old .- bp)),
                    maximum(abs.(cT_old .- cT)),
                    maximum(abs.(pN_old .- pN))  ] )

    diff2=round(diff,digits=10)

    if sum(counter .==collect(0:show_every:n_iter)) != 0
        println("Diff = $diff2 ; Iter= $counter")
    else
        nothing
    end

    if diff<tol
        println("Convergence reached :D! convergence=$diff2")
        global cT_out=copy(cT)
        global bp_out=copy(bp)
        global pN_out=copy(pN)
        global μ_out=copy(μ)
        break
    elseif counter ≥ n_iter
        println("Failed convergence D:! iteration=$max_iter")
        break
    end
    counter += 1
end

# sum(bp .< -κ*(pN.*YN .+YT)) #Collateral constraint
# sum(μ .< -tol)
# plot(μ')

# #Slackness condition: 
# TC= μ.*(bp .< -κ*(pN.*YN .+YT))
# plot(TC', title="Slackness Condition")

# #Policies
# plt_b=plot(b_grid,[bp[1,:] bp[end,:] ],
#   title="Policy function for debt b'", legend=:none, color= [:red :black],
#    xaxis="Current debt b" );

# plt_p=plot(b_grid,[pN[1,:] pN[end,:] ],
#    label=["yN low, yT low" "yN low, yT high"],
#     legend=:bottomright, title="Policy function for price pN", color= [:red :black],
#      xaxis="Current debt b" );

# plot(plt_b, plt_p, layout=(2,1))

# ---------------------------------------------------------------------------- #
#                            Social Planner Problem                            #
# ---------------------------------------------------------------------------- #
#initial values
μsp  = copy(μ)
Eμsp = copy(Eμ) 
cTsp = copy(cT)
bpsp = copy(bp)
pNsp = copy(pN)

#Technical parameters
counter=0
updt= 0.02
tol=1e-8
n_iter=1300
while counter < n_iter
    global updt,counter,n_iter,μsp ,Eμsp,cTsp,bpsp,pNsp,cTsp_bind,bpsp_bind

    #save old values
    bpsp_old = copy(bpsp)
    cTsp_old = copy(cTsp)
    pNsp_old = copy(pNsp)

    #Next binding values
    bpsp_bind = -κ.*(pNsp.*YN+YT)
    cTsp_bind = (1+r).*B .+ YT .- bpsp_bind

    #marginal utility
    #ψ   = κ*(1+η)*pNsp.*(YN./cTsp)
    ψ   = κ*(1+η)*(1-ω)/ω .* (cTsp./YN).^η

    λsp = ( ω.*cTsp.^(-η) .+ (1-ω).*YN.^(-η) ).^(σ/η-1/η-1) .* ω .* cTsp.^(-η-1) ./(1 .- μsp.*ψ)
    
    #Expected marginal utility
    for i in 1:nb
        for j in 1:ny
            Eμsp[j,i] = β*(1+r)*T[j,:]'*interp1(b_grid,λsp,bpsp[j,i])
        end
    end

    #Compute Euler's residual assuming that constraint binds
    λsp_bind = ( ω.*cTsp_bind.^(-η) .+ (1-ω).*YN.^(-η) ).^(σ/η-1/η-1) .* ω .* cTsp_bind.^(-η-1) ./(1 .- μsp.*ψ)
    μsp = 1 .- Eμsp ./ λsp_bind

    #system of equations
    bpsp[μsp .≥ tol] = -κ.*(pNsp.*YN+YT)[μsp .≥ tol]
    cTsp[μsp .≥ tol] = ((1+r).*B .+ YT .- bpsp)[μsp .≥ tol]
    pNsp[μsp .≥ tol] = ((1-ω)/ω .* (cTsp./YN).^(1+η))[μsp .≥ tol]
    
    for i in 1:nb
        for j in 1:ny
            if μsp[j,i] < -tol #then constraint binds
                u(cc)   = mgu(cc,YN[j,i]) - Eμsp[j,i] 
                cTsp[j,i] = find_zero( u , 0.1 ) ; μsp[j,i] = u(cTsp[j,i]) # (check μ ≥ 0)
                bpsp[j,i] = B[j,i]*(1+r) + YT[j,i] - cTsp[j,i]
            end
        end
    end

    #updating rule
    bpsp =  bpsp_old + updt.*(bpsp - bpsp_old)
    cTsp =  cTsp_old + updt.*(cTsp - cTsp_old)
    pNsp =  pNsp_old + updt.*(pNsp - pNsp_old)

    #convergence criterion
    diff=maximum( [ norm(bpsp_old .- bpsp),
                    norm(cTsp_old .- cTsp),
                    norm(pNsp_old .- pNsp)  ] )

    diff2=round(diff,digits=10)

    if sum(counter .==collect(0:show_every:n_iter)) != 0
        println("Diff = $diff2 ; Iter= $counter")
    else
        nothing
    end

    if diff<tol
        println("Convergence reached :D! convergence=$diff2")
        global cTsp_out = copy(cTsp)
        global bpsp_out = copy(bpsp)
        global pNsp_out = copy(pNsp)
        global μsp_out  = copy(μsp)
        break
    elseif counter ≥ n_iter
        println("Failed convergence D:! iteration=$max_iter")
        break
    end
    counter += 1
end


# sum(bpsp .< -κ*(pNsp.*YN .+YT)) #Collateral constraint
# sum(μsp .< -tol)

# #Slackness condition: 
# plot(μsp[ bpsp .> -κ*(pNsp.*YN .+YT)],
#     title="Slackness Condition", 
#     label= L"$\mu_t$ when the constraint is slack")

# plt_bsp=plot(b_grid,[bpsp_out[1,:] bpsp_out[end,:] ],
#   title="Policy function for debt b'", legend=:none, color= [:red :black],
#    xaxis="Current debt b" );

# plt_psp=plot(b_grid,[pNsp_out[1,:] pNsp_out[end,:] ],
#    label=["yN low, yT low" "yN low, yT high"],
#     legend=:bottomright, title="Policy function for price pN", color= [:red :black],
#      xaxis="Current debt b" );

# plot(plt_bsp, plt_psp, layout=(2,1) )

# ---------------------------------------------------------------------------- #
#                                  Simulations                                 #
# ---------------------------------------------------------------------------- #

pyplot(size=(900,800))
n_sims=800_000
bp_sim, bpsp_sim =Simulate!(n_sims,bp_out,bpsp_out);
b_den=density([bp_sim, bpsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
        xaxis=(L"Current debt ($b_t$)" ,[-1 ,-0.7]),color=[:red :blue],
        label=["Competitive Equilibrium" "Social Planner Solution"],
        title= L"Density of next period debt (b$_{t+1}$)")

cd(dirname(@__FILE__))
savefig(b_den, "./b_density.png")