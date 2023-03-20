using Distributions
using ProgressMeter  #to show progress
using LinearAlgebra  #norm function (compute norm of a matrix)
using JLD            #to save and load results
using Plots          #plots
using LaTeXStrings   #to write LaTeX in legends
using QuantEcon      #stationary distributions
using StatsPlots     #kernel density
using Random         #seed
using Parameters
using BenchmarkTools #test 

Threads.nthreads()
include("tauchen2.jl"); #it allows the use of real numbers in the tauchens' sigma (not just integers)

# ---------------------------------------------------------------------------- #
#                                  Calibration                                 #
# ---------------------------------------------------------------------------- #
#include("tauchenHussey.jl"); #tauchenHussey() function
#include("tauchen_var.jl"); #tauchen_var() function (correlated shocks)
#findnearest(A::AbstractArray,t) = findmin(abs.(A.-t))
# r      = 0.04             #interest rate
# σ      = 2                #risk aversion
# η      = 1/0.83-1         #elasticity of substitution
# ω      = 0.3296             #weight on tradables
# β      = 0.94             #discount factor
# b_min  = -1.2533
# b_max  = .0 
# nb     = 400
# #INCOME
# ny     = 20
# ρ_yT   = 0.2407          #1st order autocorrel of tradable
# σ_yT   = 0.0341         #sd of tradable output
# #KAPPA
# nκ     = 17
# ρ_κ    = 0.8242     #1st order autocorrel of nontradable **  0.82
# σ_κ    = 0.1064          #sd of nontradable output *
# μ_κ    = exp(-0.7836)

function Economy(;  r      = 0.04,             #interest rate
                    σ      = 2,                #risk aversion
                    η      = 1/0.83-1,         #elasticity of substitution
                    ω      = 0.3296,             #weight on tradables
                    β      = 0.94,             #discount factor
                    b_min  = -1.2533,
                    b_max  = .0, 
                    nb     = 400,
                    #INCOME
                    ny     = 20,
                    ρ_yT   = 0.2407,          #1st order autocorrel of tradable
                    σ_yT   = 0.0341,         #sd of tradable output
                    #KAPPA
                    nκ     = 17,
                    ρ_κ    = 0.8242,     #1st order autocorrel of nontradable **  0.82
                    σ_κ    = 0.1064,          #sd of nontradable output *
                    μ_κ    = exp(-0.7836)
                    )          # mean of κ AR(1) process
                    

    R         = (1+r)
    q         = 1/R 
    ns        = ny*nκ #state-space dimension
    # ---------------------------------------------------------------------------- #
    #                                   DEBT GRID                                  #
    # ---------------------------------------------------------------------------- #
    # DEBT-GRID OPTION 1: Uniform
    b_grid = collect(range(b_min,b_max,length=nb))

    # DEBT-GRID OPTION 2: Random
    # Random.seed!(3) ; x=rand(Normal(mean([b_max b_min]),1),500)
    # b_grid= sort(unique(clamp.(x,b_min,b_max)))
    # nb    =  length(b_grid)
    # density(b_grid, title="Grid points distribution")

    # ---------------------------------------------------------------------------- #
    #                             SHOCK DISCRETIZATION                             #
    # ---------------------------------------------------------------------------- #
    # OPTION A: ADiscretize shocks: Tauchen (1986)
    Πy        =      tauchen(ny,ρ_yT,σ_yT,0,2).p
    yT_grid   = exp.(tauchen(ny,ρ_yT,σ_yT,0,2).state_values) #states for nontradable
    Πκ        =      tauchen2(nκ,ρ_κ,σ_κ,log(μ_κ),1.5).p
    κ_grid    = exp.(tauchen2(nκ,ρ_κ,σ_κ,(1-ρ_κ)*log(μ_κ),1.5).state_values) #states for fin. shock
    Π         = kron(Πy,Πκ)         #transition matrix ( y move slower )
    s_grid    = gridmake(κ_grid,yT_grid) #grid for shocks ( y move slower ) =#
    

    #OPTION B: Use Tauchen & Hussey (1991) 
    # #yT process
    # Πy        = tauchenHussey(ny,ρ_yT,σ_yT)[2],#trans matrix of AR(1) for tradable
    # yT_grid = exp.(tauchenHussey(ny,ρ_yT,σ_yT)[1]), #states for tradable
    # #κT process
    # Πκ        = tauchenHussey(nκ,ρ_κ,σ_κ;mmu=μ_κ)[2],#trans matrix of AR(1) for nontradable
    # κ_grid  = (tauchenHussey(nκ,ρ_κ,σ_κ;mmu=μ_κ)[1]), #states for nontradable
    # #combining processes
    # Π         = kron(Πy,Πκ),         #transition matrix #
    # s_grid    = gridmake(κ_grid,yT_grid), #grid for shocks

    #OPTION C: tpm.jl (CORRELATED SHOCKS)
    #=
    Π,s_grid,~ = tpm(A,Σ,[ny , nκ];T=1_000_000,Tburn=100_000)
    ns = size(Π,1)
    s_grid[:,1] = exp.(s_grid[:,1])
    s_grid[:,2] = exp.(s_grid[:,2]) #.+ μ_κ

    s_grid_i  = collect(1:size(s_grid,1))
    yT_grid   = s_grid[:,1]
    κ_grid    = s_grid[:,2] #index # =#

    #grid for y and b
    K  = repeat(s_grid[:,1]',nb,1)
    YT = repeat(s_grid[:,2]',nb,1)
    B  = repeat(b_grid,1,ns)

    # ---------------------------------------------------------------------------- #
    #                                ALLOCATE MEMORY                               #
    # ---------------------------------------------------------------------------- #
    #initial values: COMPETITIVE
    BP = repeat(b_grid,1,ns)
    CT = ones(nb,ns)
    C  = ones(nb,ns)
    PN = ones(nb,ns)
    BC = ones(nb,ns)
    V  = zeros(nb,ns)
    μ  = zeros(nb,ns)
    λ  = zeros(nb,ns)

    #initial values: PLANNER
    BPsp = repeat(b_grid,1,ns)
    CTsp = ones(nb,ns)
    Csp  = ones(nb,ns)
    PNsp = ones(nb,ns)
    BCsp = ones(nb,ns)
    Vsp  = zeros(nb,ns)
    μsp  = zeros(nb,ns)
    λsp  = zeros(nb,ns)

    #initial values: PLANNER
    BPir = repeat(b_grid,1,ns)
    CTir = ones(nb,ns)
    Cir  = ones(nb,ns)
    PNir = ones(nb,ns)
    BCir = ones(nb,ns)
    Vir  = zeros(nb,ns)
    μir  = zeros(nb,ns)
    λir  = zeros(nb,ns)
    APir = repeat(b_grid,1,ns)
    TAX  = zeros(nb,ns)


    # ---------------------------------------------------------------------------- #
    #                                   Functions                                  #
    # ---------------------------------------------------------------------------- #
    u(ec,c) =( c^(1-ec.σ) -1 ) / (1 - ec.σ)
    #mgu(ec,ct,yn)= ct > 0 ? ( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(ec.σ/ec.η-1/ec.η-1) * ec.ω * ct^(-ec.η-1) : 1e15
    mgu(ec,ct,yn)=( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(ec.σ/ec.η-1/ec.η-1) * ec.ω * ct^(-ec.η-1)
    CES(ec,ct,yn) = ( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(-1/ec.η)

    return (β = β, σ = σ, r = r,R=R, q=q, η = η , ω=ω, nκ=nκ,  ny = ny, nb = nb, ns=ns,
            ρ_κ=ρ_κ,  σ_κ=σ_κ,  μ_κ=μ_κ,  ρ_yT=ρ_yT,  σ_yT=σ_yT,
            Π = Π, u=u, mgu=mgu, CES=CES,
            κ_grid=κ_grid,yT_grid=yT_grid,b_grid=b_grid,s_grid=s_grid, K = K,YT = YT, B = B,
            BP = BP    ,BC=BC    , CT = CT    , PN = PN,
            BPsp = BPsp,BCsp=BCsp, CTsp = CTsp, PNsp = PNsp,
            BPir = BPir,BCir=BCir, CTir = CTir, PNir = PNir, APir=APir,TAX=TAX,
            V = V,Vsp = Vsp,Vir = Vir, μ=μ, μsp=μsp , μir=μir,λ=λ, λsp=λsp, λir=λir)
end

#∂u∂ct(ct,yn) = ct > 0 ? ( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(ec.σ/ec.η-1/ec.η-1) * ec.ω * ct^(-ec.η-1) : 1e15
#∂u∂ct(ct,yn) = ct > 0 ? ( ω*ct^(-η) + (1-ω)*yn^(-η) )^(σ/η-1/η-1) * ω * ct^(-η-1) : 1e15

# ---------------------------------------------------------------------------- #
#                            Competitive Equilibrium                           #
# ---------------------------------------------------------------------------- #

# ec = Economy()
# ec.q
# ec.K
# mean(ec.K)
# Eλ = ec.λ*ec.Π'
# Policy_operator!(ec,Eλ)
# is=11
# ib=1
# ib_next=1
function Policy_operator!(ec,Eλ)
    # unpack stuff
    @unpack β, σ, q, η, ω, ny, nb, ns, r = ec
    @unpack PN, BP, CT, BC, V, λ, μ = ec
    @unpack κ_grid, yT_grid,s_grid ,b_grid,K, YT, B, Π = ec
    @unpack u, CES,mgu = ec
    euler_res = zeros(nb)
    slackness = zeros(nb)
    # pol_ind   = Any
    @fastmath @inbounds begin #Threads.@threads 
        for ib in 1:nb
            b=b_grid[ib]
            for is in 1:ns
                κ,yt = s_grid[is,:]
                yn   = 1
                pn   = PN[ib,is]
                borr_const = -κ*(pn*yn + yt)/q #Threads.@threads 
                Threads.@threads for ib_next in 1:nb

                    #Borrowing Constraint
                    if b_grid[ib_next] > borr_const #not binding
                        bp  = copy(b_grid[ib_next])
                    else #binding
                        bp  = copy(borr_const)
                    end

                    #Consumption > 0 
                    if  yt + b - q*bp > 1e-15 
                        ct  = yt + b - q*bp
                    else
                        ct  = copy(1e-15)
                    end
                    euler_res[ib_next] = mgu(ec,ct,yn) - β*Eλ[ib_next,is]/q
                    # slackness[ib_next] = euler_res[ib_next]*(bp - borr_const)
                    # #slackness[ib_next] = bp

                    #Opt2
                    # bp  = max(b_grid[ib_next],borr_const)
                    # ct  = max(yt + b - q*bp , 1e-15)
                    # euler_res[ib_next] = mgu(ec,ct,yn) - β*Eλ[ib_next,is]/q
                end
                #update
                #pol_ind      = min(searchsortedfirst(euler_res,0),nb)
                #pol_ind      = min(searchsortedfirst(slackness,0),nb)
                pol_ind      = findmin(abs.(euler_res))[2] #find zero
                #pol_ind      = findmin(abs.(slackness))[2] #find zero

                ec.BP[ib,is] = max(b_grid[pol_ind],borr_const)
                ec.CT[ib,is] = max(yt + b - ec.BP[ib,is]*q ,1e-15)
                ec.PN[ib,is] = (1-ω)/ω * ( ec.CT[ib,is]/yn ) ^(η+1)
                ec.BC[ib,is] = borr_const
                ec.λ[ib,is]  = mgu(ec,CT[ib,is],yn)
                ec.μ[ib,is]  = euler_res[pol_ind]
            end
        end
    end
end

function PFI!(ec;tol=1e-5,maxit=300,show_iter=10,updt=0.1)
    BP_0 = copy(ec.BP)
    CT_0 = copy(ec.CT)
    PN_0 = copy(ec.PN)
    DIFF = Any[]
    @unpack Π,q,β = ec
    #@fastmath @inbounds begin
        for it in 1:maxit
            copyto!(BP_0, ec.BP + updt.*(BP_0 - ec.BP))
            copyto!(CT_0, ec.CT + updt.*(CT_0 - ec.CT))
            copyto!(PN_0, ec.PN + updt.*(PN_0 - ec.PN))

            Eλ = ec.λ * Π'
            Policy_operator!(ec,Eλ)
            
            #dist = maximum(abs(x - y) for (x, y) in zip([BP_0 CT_0 PN_0], [ec.BP ec.CT ec.PN]))
            dist=maximum( [ maximum(abs.(BP_0 .- ec.BP)),
                            maximum(abs.(CT_0 .- ec.CT)),
                            maximum(abs.(PN_0 .- ec.PN))  ] )

            dist < tol ? break : nothing
            #check if iteration are converging (if not try a smaller step)
            DIFF=push!(DIFF,dist)
            it > 1 ?  (DIFF[it] > DIFF[it-1] ? updt=updt*(1-updt) : nothing) : nothing

            it>30 ? (DIFF[it-30] ≤ DIFF[it] ? break : nothing ) : nothing 

            it % show_iter == 0 ? println("Iteration $(it) ; Dist $(round(dist,digits=11)) ; Step $updt") : nothing
        end
    #end
end


# ---------------------------------------------------------------------------- #
#                                Social Planner                                #
# ---------------------------------------------------------------------------- #

function T_operator!(ec,EVsp)
    @unpack β, σ, q, η, ω, ny, nb,ns = ec
    @unpack PNsp, BPsp, CTsp, BCsp, Vsp = ec
    @unpack κ_grid, yT_grid,s_grid ,b_grid,K, YT, B, Π = ec
    @unpack u, CES = ec

    Vprime = zeros(size(b_grid))
    @fastmath @inbounds begin #Threads.@threads 
        for ib in 1:nb
            b = b_grid[ib]
            for is in 1:ns
                κ,yt = s_grid[is,:]
                #yt,κ = s_grid[is,:]
                yn = 1
                pn = PNsp[ib,is]
                #borr_const = κ*(pn*yn + yt)/q
                borr_const = -κ*( (1-ω)/ω *( ec.CTsp[ib,is]/yn ) ^(η+1) *yn + yt)/q
                Threads.@threads for ib_next in 1:nb
                    bp  = max(b_grid[ib_next],borr_const)
                    ct  = max(yt + b - q*bp,1e-15) 
                    c   = max(CES(ec,ct,yn),1e-15)
                    Vprime[ib_next] = u(ec, c) + β*EVsp[ib_next,is]
                end
                #save values
                ec.Vsp[ib,is], pol_ind  = findmax(Vprime)
                ec.BPsp[ib,is]          = max(b_grid[pol_ind],borr_const)
                ec.CTsp[ib,is]          = max(yt + b - ec.BPsp[ib,is]*q ,1e-15) 
                ec.PNsp[ib,is]          = ((1-ω)/ω) *( ec.CTsp[ib,is]/yn ) ^(η+1)
                ec.BCsp[ib,is]          = borr_const
            end
        end
        #Reserves Economy
        #copyto!(ec.BPir, -ec.K.*(ec.PNsp .+ ec.YT)/q)
        #copyto!(ec.BPir, -ec.K.*( (1-ω)/ω .*( ec.CTsp ) .^(η+1) .+ ec.YT)/q)
        copyto!(ec.BPir, ec.BCsp)  
        copyto!(ec.TAX , (ec.BPsp - ec.BPir)*q)
        copyto!(ec.APir, ec.TAX / q)
    end
end

function VFI!(ec;tol=1e-8,maxit=500,show_iter=10)
    Vsp_0 = similar(ec.Vsp)
    @unpack Π = ec
    #@fastmath @inbounds begin #
        for it in 1:maxit

            copyto!(Vsp_0, ec.Vsp) #save old value

            EVsp = ec.Vsp*Π'
            T_operator!(ec,EVsp) #update value
            
            dist = maximum(abs(x - y) for (x, y) in zip(Vsp_0, ec.Vsp))
            dist < tol ? break : nothing
            it % show_iter == 0 ? println("Finished iteration $(it) with dist of $(round(dist,digits=11))") : nothing
        end
    #end
end

# ---------------------------------------------------------------------------- #
#                                Solve the model                               #
# ---------------------------------------------------------------------------- #

#ec = Economy(nb=800,ny=14,nκ=14,b_max  = 1.5,σ_κ=0.011,ρ_κ=0.0,μ_κ=0.25)
ec = Economy(nb=400)

# 1. Solve competitive equilibrium
@time PFI!(ec;tol=1e-2,maxit=200,show_iter=1,updt=0.1)

# 2. solving planner's problem
@time VFI!(ec;tol=1e-5,maxit=300,show_iter=1)

# ---------------------------------------------------------------------------- #
#                             Plot Policy function                             #
# ---------------------------------------------------------------------------- #
plot(-ec.b_grid,-ec.BP/mean(ec.PN + ec.YT),legend=:none,yaxis=[0.05 ,0.9])
plot(-ec.b_grid,-ec.BPsp/mean(ec.PNsp   + ec.YT),legend=:none,yaxis=[0.05 ,0.9])
plot(-ec.b_grid,-ec.BPir/mean(ec.PNsp + ec.YT),legend=:none,yaxis=[0.05 ,0.9])

# Check collateral constraint
sum(ec.BC .≤ ec.BP)
sum(ec.BCsp .≤ ec.BPsp)
sum(ec.BCsp .≈ ec.BPir) == prod([ec.nb ec.ns])

#Density
# bp_sim=Simulate1!(ec,1000_000,ec.BP)
# density(bp_sim)

# bpsp_sim=Simulate1!(ec,1000_000,ec.BPsp)
# density(bpsp_sim)

#output=deepcopy(ec)
#save(string(pwd(), "/solution_shock_k.jld"), "output", output)

# ---------------------------------------------------------------------------- #
#                                   Simulate                                   #
# ---------------------------------------------------------------------------- #


function Simulate!(ec;GDP=0,n_sims=1_000_000)
    cut=Int(floor(n_sims*0.2))

    function mc_sample_path(P; init = 1, sample_size = 1000)
        @assert size(P)[1] == size(P)[2] # square required
        N = size(P)[1] # should be square
    
        # create vector of discrete RVs for each row
        dists = [Categorical(P[i, :]) for i in 1:N]
    
        # setup the simulation
        X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
        X[1] = init # set the initial state
    
        for t in 2:sample_size
            dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
            X[t] = rand(dist) # draw new value
        end
        return X
    end

    s_shock = mc_sample_path(ec.Π;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(s_shock))
    bpsp_sim = zeros(length(s_shock))
    bpir_sim = zeros(length(s_shock))

    b_t   = mean(ec.b_grid)
    bsp_t = mean(ec.b_grid)
    bir_t = mean(ec.b_grid)
    
    gdp_sp     = mean(ec.PN   + ec.YT)
    gdp_ce     = mean(ec.PNsp + ec.YT)
    gdp_ir     = mean(ec.PNsp + ec.YT)
    @showprogress for i in 1:length(s_shock)
        #Simulate y-shocks
        is  = s_shock[i]
        ib   = clamp(searchsortedfirst(ec.b_grid,b_t),1,ec.nb)
        ibsp = clamp(searchsortedfirst(ec.b_grid,bsp_t),1,ec.nb)
        ibir = clamp(searchsortedfirst(ec.b_grid,bir_t),1,ec.nb)

        #Competitive Equilibrium Simulations
        bp_sim[i]     = ec.BP[ib,is] #./((ec.PN[ib,is] + ec.YT[ib,is]))
        #Social Planner Simulations
        bpsp_sim[i]   = ec.BPsp[ibsp,is] #./((ec.PNsp[ibsp,is] + ec.YT[ibsp,is]))
        #Reserves Economy
        bpir_sim[i]   = ec.BPir[ibsp,is] #./((ec.PNir[ibir,is] + ec.YT[ibir,is]))

        b_t   = copy(bp_sim[i])
        bsp_t = copy(bpsp_sim[i] )
    end
    if GDP == 0
        return bp_sim[cut:end], bpsp_sim[cut:end],bpir_sim[cut:end]
    elseif GDP == 1
        return bp_sim[cut:end]/gdp_ce, bpsp_sim[cut:end]/gdp_sp, bpir_sim[cut:end]/gdp_ir
    end
end


pyplot(size=(900,800))
n_sims=1000_000

bp_sim, bpsp_sim,bpir_sim =Simulate!(ec;n_sims,GDP=1);
b_den=density(-ec.b_grid,[-bp_sim, -bpsp_sim, -bpir_sim], yaxis= ("Frequency"),lw=2.,
        xaxis=(L"Debt to GDP ($b_t/(pN_t yN_t + yT_t)$)",[0.1, 0.8]),color=[:red :blue :green],
        label=["Competitive Equilibrium" "Social Planner Solution" "Reserves Economy"],
        title= L"Density of next period debt (b$_{t+1}/(pN_t*yN_t + yT_t)$)")

# ---------------------------------------------------------------------------- #
#                                    OTHERS                                    #
# ---------------------------------------------------------------------------- #

# function Simulate1!(ec,n_sims;cut=Int(floor(n_sims*0.2)))
#     function mc_sample_path(P; init = 1, sample_size = 1000)
#         @assert size(P)[1] == size(P)[2] # square required
#         N = size(P)[1] # should be square
    
#         # create vector of discrete RVs for each row
#         dists = [Categorical(P[i, :]) for i in 1:N]
    
#         # setup the simulation
#         X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
#         X[1] = init # set the initial state
    
#         for t in 2:sample_size
#             dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
#             X[t] = rand(dist) # draw new value
#         end
#         return X
#     end

#     y_shock = mc_sample_path(ec.Π;sample_size = n_sims+cut)#Simulations
#     bp_sim  = zeros(length(y_shock))
#     b_t     = mean(ec.b_grid)
#     gdp     = mean(ec.PN + ec.YT)
#     @showprogress for i in 2:length(y_shock)
#         #Simulate y-shocks
#         iy=y_shock[i]
#         ib   = searchsortedfirst(ec.b_grid,b_t)
#         bp_sim[i]   = ec.BP[ib,iy]

#         b_t   = copy(bp_sim[i])
#     end
#     return bp_sim[cut:end]
# end

# bp_sim=Simulate1!(ec,900_000)
# gdp     = mean(ec.PN + ec.YT)
# density(ec.b_grid,-bp_sim/gdp,xaxis=[0.1 ,0.8],color=:blue,)

# function stat_dist(ec)
#     λ=eigvals(ec.Π)[size(ec.Π,1)]
#     Reig=eigvecs(ec.Π) #right eigenvalues
#     Leig=inv(Reig)  #left  eigenvalues
#     π_bar = real.(vec(Leig[size(ec.Π,1),:]))
#     #P*v = λ*v; P*v = v (λ=1)
#     @assert vec(π_bar'*ec.Π) ≈ λ*π_bar
#     @assert sum(π_bar./sum(π_bar)) ≈ 1
#     sdist=π_bar./sum(π_bar)
#     return sdist
# end