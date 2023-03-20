# ---------------------------------------------------------------------------- #
#                         Replication of Bianchi (2011)                        #
# ---------------------------------------------------------------------------- #
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
using BenchmarkTools

#Load stochastic structure

Threads.nthreads()
# ---------------------------------------------------------------------------- #
#                                  Calibration                                 #
# ---------------------------------------------------------------------------- #
include("calibration.jl")

function Economy(;  r  = 0.04,             #interest rate
                    σ  = 2,                #risk aversion
                    η  = 1/0.83-1,         #elasticity of substitution
                    ω  = 0.307000693252802,             #weight on tradables
                    β  = 0.906,             #discount factor
                    κ  = 0.3235,
                    b_min  = -1.02,
                    b_max  = -0.4, 
                    nb     =  400,
                    yT_grid=yT_grid,
                    yN_grid=yN_grid,
                    Π = T)

    #OPTION 1
    b_grid = collect(range(b_min,b_max,length=nb))

    #OPTION 2
    # Random.seed!(3) ; x=rand(Normal(mean([b_max b_min]),1),500)
    # b_grid= sort(unique(clamp.(x,b_min,b_max)))
    # nb    =  length(b_grid)
    # density(b_grid, title="Grid points distribution")

    #grid for y and b
    ny = length(yN_grid)
    YN = repeat(yN_grid',nb,1)
    YT = repeat(yT_grid',nb,1)
    B  = repeat(b_grid,1,ny)

    #initial values: COMPETITIVE
    BP = repeat(b_grid,1,ny)
    CT = ones(nb,ny)
    C  = ones(nb,ny)
    PN = ones(nb,ny)
    BC = ones(nb,ny)
    V  = zeros(nb,ny)
    μ  = zeros(nb,ny)
    λ  = zeros(nb,ny)

    #initial values: PLANNER
    BPsp = repeat(b_grid,1,ny)
    CTsp = ones(nb,ny)
    Csp  = ones(nb,ny)
    PNsp = ones(nb,ny)
    BCsp = ones(nb,ny)
    Vsp  = zeros(nb,ny)
    μsp  = zeros(nb,ny)
    λsp  = zeros(nb,ny)


    #Functions
    u(ec,c) =( c^(1-ec.σ) -1 ) / (1 - ec.σ)
    mgu(ec,ct,yn)=( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(ec.σ/ec.η-1/ec.η-1) * ec.ω * ct^(-ec.η-1)
    CES(ec,ct,yn) = ( ec.ω*ct^(-ec.η) + (1-ec.ω)*yn^(-ec.η) )^(-1/ec.η)

    return (β = β, σ = σ, r = r, η = η ,κ=κ , ω=ω,  ny = ny, nb = nb,
            Π = Π, u=u, mgu=mgu, CES=CES,
            yN_grid=yN_grid,yT_grid=yT_grid,b_grid=b_grid, YN = YN,YT = YT, B = B,
            BP = BP,BC=BC, CT = CT, PN = PN,
            BPsp = BPsp,BCsp=BCsp, CTsp = CTsp, PNsp = PNsp,
            V = V,Vsp = Vsp, μ=μ, μsp=μsp ,λ=λ, λsp=λsp)
end

# ---------------------------------------------------------------------------- #
#                            Competitive Equilibrium                           #
# ---------------------------------------------------------------------------- #
# ib=200
# iy=1
function Policy_operator!(ec,Eμ)
    # unpack stuff
    @unpack β, σ, r, η, κ, ω, ny, nb = ec
    @unpack YN, YT, B, Π, PN, BP, CT, BC, V, λ, μ = ec
    @unpack yN_grid, yT_grid ,b_grid = ec
    @unpack u, CES,mgu = ec

    policy_candidates = zeros(nb)

    #@fastmath @inbounds begin #Threads.@threads 
        for ib in 1:nb
            b=b_grid[ib]
            for iy in 1:ny
                yt = yT_grid[iy]
                yn = yN_grid[iy]
                pn = PN[ib,iy]
                #ct = CT[ib,iy]
                borr_const = -κ*(pn*yn + yt)
                for ib_next in 1:nb
                    bp = max(b_grid[ib_next],borr_const)
                    ct  = yt + b*(1+r) - bp 
                    policy_candidates[ib_next] = mgu(ec,ct,yn) - Eμ[ib_next,iy]
                    #policy_candidates[ib_next] = ct
                end
                #update
                #pol_ind      = min(searchsortedfirst(policy_candidates,0),nb)
                pol_ind = findmin(abs.(policy_candidates))[2]

                ec.BP[ib,iy] = max(b_grid[pol_ind],borr_const)
                ec.CT[ib,iy] = yt + b*(1+r) - ec.BP[ib,iy]
                ec.PN[ib,iy] = (1-ω)/ω * ( max(ec.CT[ib,iy],0)/yn ) ^(η+1)
                ec.BC[ib,iy] = borr_const
                ec.λ[ib,iy]  = mgu(ec,CT[ib,iy],yn)
                ec.μ[ib,iy] = policy_candidates[pol_ind]
            end
        end
    #end
end


#plot(abs.(policy_candidates))

function PFI!(ec;tol=1e-4,maxit=300,show_iter=10,updt=0.1)
    BP_0 = copy(ec.BP)
    CT_0 = copy(ec.CT)
    PN_0 = copy(ec.PN)
    DIFF = Any[]
    @unpack Π,r,β = ec
    #@fastmath @inbounds begin
    for it in 1:maxit
        copyto!(BP_0, ec.BP + updt.*(BP_0 - ec.BP))
        copyto!(CT_0, ec.CT + updt.*(CT_0 - ec.CT))
        copyto!(PN_0, ec.PN + updt.*(PN_0 - ec.PN))

        Eμ = β*(1+r)*ec.λ*Π'
        Policy_operator!(ec,Eμ)
        
        #dist = maximum(abs(x - y) for (x, y) in zip([BP_0 CT_0 PN_0], [ec.BP ec.CT ec.PN]))
        dist=maximum( [ maximum(abs.(BP_0 .- ec.BP)),
                        maximum(abs.(CT_0 .- ec.CT)),
                        maximum(abs.(PN_0 .- ec.PN))  ] )

        dist < tol ? break : nothing
        #check if iteration are converging (if not try a smaller step)
        DIFF=push!(DIFF,dist)
        it > 1 ?  (DIFF[it] > DIFF[it-1] ? updt=updt*(1-updt) : nothing) : nothing
        it % show_iter == 0 ? println("Iteration $(it) with dist of $(round(dist,digits=11)) and a step of $updt") : nothing
    end
    #end
end

# ---------------------------------------------------------------------------- #
#                                Social Planner                                #
# ---------------------------------------------------------------------------- #

#initialize with optimal values
# copyto!(ec.BPsp,ec.BP)
# copyto!(ec.CTsp,ec.CT)
# copyto!(ec.PNsp,ec.PN)
# copyto!(ec.BCsp,ec.BC)
# copyto!(ec.Vsp,ec.V)


function T_operatorSP!(ec,EVsp)
# unpack stuff
    @unpack β, σ, r, η, κ, ω, ny, nb = ec
    @unpack PNsp, BPsp, CTsp, BCsp, Vsp = ec
    @unpack yN_grid, yT_grid ,b_grid,YN, YT, B, Π = ec
    @unpack u, CES = ec

    policy_candidates = zeros(size(b_grid))
    @fastmath @inbounds begin
        for ib in 1:nb
            b=b_grid[ib]
            for iy in 1:ny
                yt = yT_grid[iy]
                yn = yN_grid[iy]
                pn = PNsp[ib,iy]
                #borr_const = -κ*(pn   *yn + yt)
                borr_const = -κ*( (1-ω)/ω *( max(ec.CTsp[ib,iy],0)/yn ) ^(η+1) *yn + yt)
                BCsp[ib,iy] = borr_const
                Threads.@threads for ib_next in 1:length(b_grid)
                    bp = b_grid[ib_next] ≥ borr_const ? b_grid[ib_next] : borr_const
                    ct  = yt + b*(1+r) - bp 
                    c = CES(ec,ct,yn) > 1e-15 ? CES(ec,ct,yn) : 1e-15
                    policy_candidates[ib_next] = u(ec, c) + β*EVsp[ib_next,iy]
                end
                ec.Vsp[ib,iy], pol_ind  = findmax(policy_candidates)
                ec.BPsp[ib,iy] = max(b_grid[pol_ind],borr_const)
                ec.CTsp[ib,iy] = yt + b*(1+r) - ec.BPsp[ib,iy]
                ec.PNsp[ib,iy] = ((1-ω)/ω) *( max(ec.CTsp[ib,iy],0)/yn ) ^(η+1)
            end
        end
    end
end
    
function VFIsp!(ec;tol=1e-5,maxit=300,show_iter=10)
    Vsp_0 = similar(ec.Vsp)
    @unpack Π = ec
    #@fastmath @inbounds begin
        for it in 1:maxit

            copyto!(Vsp_0, ec.Vsp)
            EVsp = ec.Vsp*Π'
            T_operatorSP!(ec,EVsp)
            
            dist = maximum(abs(x - y) for (x, y) in zip(Vsp_0, ec.Vsp))
            dist < tol ? break : nothing
            it % show_iter == 0 ? println("Finished iteration $(it) with dist of $(round(dist,digits=11))") : nothing

        end
    #end
end

# ---------------------------------------------------------------------------- #
#                                    Solving                                   #
# ---------------------------------------------------------------------------- #

ec = Economy(nb=400)

#Solution: COMPETITIVE
@btime PFI!(ec;tol=1e-3,maxit=200,show_iter=1,updt=0.5)

#Solution: PLANNER
@btime VFIsp!(ec;tol=1e-8,maxit=1000,show_iter=1)

#Checking Collateral Constraint
sum(ec.BP .< ec.BC)
sum(ec.BPsp .< ec.BCsp)

#Plot policies
bp    = ec.BP
bpsp  = ec.BPsp
b_grid= ec.b_grid

# neg_shock_yT = 0.95*mean(ec.yT_grid)
# iyt_shock =searchsortedfirst(yT_grid,neg_shock_yT)
plot(b_grid,bp, legend= :none)

plot(b_grid,[bp[:,1] bpsp[:,1] ], label = ["Competitive" "Planner"])

# ---------------------------------------------------------------------------- #
#                                  Simulations                                 #
# ---------------------------------------------------------------------------- #



function Simulate!(ec,n_sims;cut=Int(floor(n_sims*0.2)))
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

    y_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(y_shock))
    bpsp_sim = zeros(length(y_shock))

    b_t   = mean(ec.b_grid)
    bsp_t = mean(ec.b_grid)
    
    @showprogress for i in 2:length(y_shock)
        #Simulate y-shocks
        iy=y_shock[i]
        # ib   = searchsortedfirst(ec.b_grid,b_t)
        # ibsp = searchsortedfirst(ec.b_grid,bsp_t)

        ib   = clamp(searchsortedfirst(ec.b_grid,b_t),1,ec.nb)
        ibsp = clamp(searchsortedfirst(ec.b_grid,bsp_t),1,ec.nb)

    
        #Competitive Equilibrium Simulations
        bp_sim[i]   = ec.BP[ib,iy] / (4*(ec.PN[ib,iy]*ec.YN[ib,iy] + ec.YT[ib,iy]))

        #Social Planner Simulations
        bpsp_sim[i]   = ec.BPsp[ibsp,iy] / (4*(ec.PNsp[ibsp,iy]*ec.YN[ibsp,iy] + ec.YT[ibsp,iy]))

        b_t   = copy(bp_sim[i])
        bsp_t = copy(bpsp_sim[i] )
    end
    return bp_sim[cut:end], bpsp_sim[cut:end]
end

pyplot(size=(900,800))
n_sims=800_000
bp_sim, bpsp_sim =Simulate!(ec,n_sims);
b_den=density(ec.b_grid,[bp_sim, bpsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
        xaxis=(L"Current debt ($b_t$)"),color=[:red :blue],
        label=["Competitive Equilibrium" "Social Planner Solution"],
        title= L"Density of next period debt (b$_{t+1}$)")

#cd(dirname(@__FILE__))
#savefig(b_den, "./b_density.png")
# ---------------------------------------------------------------------------- #
#                               Simulate 1 policy                              #
# ---------------------------------------------------------------------------- #
function Simulate1!(ec,n_sims,bp,;cut=Int(floor(n_sims*0.2)))
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

    y_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    bp_sim  = zeros(length(y_shock))
    b_t     = mean(ec.b_grid)
    
    @showprogress for i in 2:length(y_shock)
        #Simulate y-shocks
        iy=y_shock[i]
        ib   = searchsortedfirst(ec.b_grid,b_t)
        bp_sim[i]   = bp[ib,iy]

        b_t   = copy(bp_sim[i])
    end
    return bp_sim[cut:end]
end
