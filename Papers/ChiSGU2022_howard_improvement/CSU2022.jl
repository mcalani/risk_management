# ---------------------------------------------------------------------------- #
#                       Chi, Schmidt-Grohe & Uribe (2022)                      #
# ---------------------------------------------------------------------------- #
using Roots
using NLsolve
using Plots
using ProgressMeter  #to show progress
using LaTeXStrings   #to write LaTeX in legends
using StatsPlots     #kernel density
using LinearAlgebra
using SparseArrays   #solucion del planner
using Optim
using Random
using JLD            #to save and load results

cd(dirname(@__FILE__))
include("tpm.jl")
include("aux_funcs.jl")

# Threads.nthreads()

# NOTA: El equilibrio competitivo está restringido en el nivel de loans and Reserves
#       segun el nivel depósitos D (dado D, decide L y R)
#       El planner no está restringido en Reserves porque puede establecer i_r.
#       por lo que primero determinar su nivel de L y R, y luego D

# ---------------------------------------------------------------------------- #
#                                  CALIBRATION                                 #
# ---------------------------------------------------------------------------- #

#calibration: structural parameters
const σ       = 2
const a       = 0.31
const ξ       = 0.83
const i_ast   = 0.04
const β       = 0.8667
const κ       = 0.3205
const A       = 0.0089
const α       = 1.8104
const ϕ       = 6.7983
const r_bar   = 0.5848
const B       = 2.6852

#calibration: disctretization space
const n_yT        = 6 #13
const n_yN        = 6 #13
const ny          = n_yT*n_yN
const d_min,d_max = [0.4,1.2] # {0.4 1.05} está bien para el competitive, pero el planner necesita {0.4 1.2}

#OPTION 1: UNIFORM GRID
const nd          = 200 #800
const d_grid   = collect(range(d_min,d_max,length=nd))
# const dsp_grid   = collect(range(d_min,1.3,length=nd))

#OPTION 2: NOT UNIFORM GRID
# Random.seed!(3) ; x=rand(Normal(mean([d_max d_min]),1),1000)
# d_grid= sort(unique(clamp.(x,d_min,d_max)))
# nd    =  length(d_grid)
# density(d_grid, title="Grid points distribution")


#calibration: transition matrix (bianchi2011 parameters)
const ρ = [ 0.901 0.495 ; -0.453 0.225]           ;
const Σ = [0.00219 0.00162  ; 0.00162 0.00167]    ;
const N = vec([n_yT n_yN])
const T = 1_000_000
const Tburn = Int(T*0.1)

const Π,S,~ = tpm(ρ,Σ,N;T=T, n_sig=1.5,Tburn=Tburn)

const yN_grid  = exp.(S[:,1])
const yT_grid  = exp.(S[:,2])
#state-variables
const YN = repeat(yN_grid',nd,1)
const YT = repeat(yT_grid',nd,1)
const D  = repeat(d_grid,1,ny)
const DD = repeat(d_grid,1,nd) #for computing expectations

# natural debt limit
const NdL = minimum(yT_grid).*(1+i_ast)./i_ast

# ---------------------------------------------------------------------------- #
#                            I. Unregulated Economy                            #
# ---------------------------------------------------------------------------- #

function Unregulated_Economy(;n_iter=5)
    i_d = copy(i_ast)
    i_r = 0.

    #initial values: COMPETITIVE
    DP = repeat(d_grid,1,ny)
    CT = YT .- D.*(1+i_ast) + DP
    PN = (1-a)/a .* ( max.(CT,0)./YN ) .^(1/ξ)
    BC = κ.*(PN.*YN .+ YT)
    μ  = zeros(nd,ny)
    λ  = ∂u∂cT.(CT,YN)#zeros(nd,ny)

    L  = zeros(Float64,nd,ny)
    R  = zeros(Float64,nd,ny)
    IL = zeros(Float64,nd,ny)

    Eλ = λ*Π' #init value

    #iteration parameters
    DP_0 = copy(DP)
    CT_0 = copy(CT)
    PN_0 = copy(PN)

    DIFF   = Float64[]
    n_iter = n_iter
    updt   = 0.0
    tol    = 1e-4

    @fastmath @inbounds begin 
        for iter in 1:n_iter #Threads.@threads 

            next_period_euler_residuals = zeros(nd)
            #begin
            Threads.@threads for id in 1:nd #@fastmath @inbounds @simd Threads.@threads 
                d=d_grid[id]
                Threads.@threads for iy in 1:ny
                    yt = yT_grid[iy]
                    yn = yN_grid[iy]
                    
                    pn         = PN[id,iy]
                    borr_const = κ*(pn*yn + yt)
                    BC[id,iy]  = borr_const
                    i_l        = IL[id,iy]

                    # dp_g  = min(DP[id,iy],borr_const) # dont need min() here
                    dp_g  = DP[id,iy] # dont need min() here
                    lt,rt = get_L_R_CE(dp_g)

                    ##### ##### ##### Expectation ##### ##### #####
                    # OPTION 1 (element-wise)
                    # for id_next in 1:nd #Threads.@threads 
                    #     dp_g  = min(d_grid[id_next],borr_const)
                    #     lt,rt = get_L_R_CE(dp_g)
                    #     if lt ≤ borr_const
                    #         l_grid[id_next]   = lt
                    #         r_grid[id_next]   = rt
                    #     elseif lt > borr_const
                    #         lt = copy(borr_const)
                    #         find_rt(r_x) = dp_g - (lt + r_x + Γ(lt,r_x)) 
                    #         rt  = fzero( find_rt , 0.1 )
                    #         l_grid[id_next]   = lt
                    #         r_grid[id_next]   = rt
                    #     end
                    #     ct  = max(yt - d*(1+i_ast) + dp_g - Γ(lt,rt) - Γr(rt),1e-15)
                    #     next_period_euler_residuals[id_next] = 1- β*(1+i_l)*Eλ[id_next,iy]/∂u∂cT(ct,yn)
                    # end
                    
                    # # OPTION 2 (using vectors->faster)

                    # dp_vec = copy(d_grid) # TEMP
                    dp_vec = min.(d_grid,borr_const)
                    # lt,rt  = Get_L_R_CE(dp_vec)
                    ct_vec = max.(yt - d*(1+i_ast) .+ dp_vec .- Γ.(lt,rt) .- Γr.(rt),1e-15)
                    next_period_euler_residuals = 1 .- β*(1+i_l)*Eλ[:,iy]./ ∂u∂cT.(ct_vec,yn)
                    ##### ##### ##### ##### ##### ##### ##### ##### 

                    #update
                    pol_ind   = findmin((next_period_euler_residuals).^2)[2]
                    μ[id,iy]  = next_period_euler_residuals[pol_ind]

                    # {L,R} or D first?  : COMPETITIVE DECIDE D, THEN L,R
                    DP[id,iy] = min(d_grid[pol_ind],borr_const) #need this
                    # DP[id,iy] = d_grid[pol_ind] 
                    # L[id,iy]  = l_grid[pol_ind]
                    # R[id,iy]  = r_grid[pol_ind]
                end
            end
    
            #Compute optimal values ({L,R} or D first?)
            # DP = L + R + Γ.(L,R)
            L,R = GET_L_R_CE(DP)  


            CT  = YT - D.*(1+i_ast) + DP - Γ.(L,R) - Γr.(R)
            PN  = (1-a)/a .* ( max.(CT,zeros(nd,ny))./YN ) .^(1/ξ)
            λ   = ∂u∂cT.(CT,YN)
            IL  = (1+i_d) .* Γ_l.(L,R) .+ i_d

            #convergence
            dist=maximum( [ maximum(abs.(DP_0 .- DP)),
                            maximum(abs.(CT_0 .- CT)),
                            maximum(abs.(PN_0 .- PN))  ] )

            dist < tol ? break : println("Iter: $iter ; Conv: $dist")
            # ProgressMeter.update!(prog, dist)
            # ProgressMeter.next!(prog, spinner="🌑🌒🌓🌔🌕🌖🌗🌘")

            #check if iteration are converging (if not try a smaller step)
            DIFF=push!(DIFF,dist)

            #again
            copyto!(DP_0, DP + updt.*(DP_0 - DP))
            copyto!(CT_0, CT + updt.*(CT_0 - CT))
            copyto!(PN_0, PN + updt.*(PN_0 - PN))

            #update expectations
            Eλ = λ*Π'
        end 
    end
    return DP,CT,PN,L,R,IL,BC,μ,i_d,i_r
end

# using BenchmarkTools
@time DP,CT,PN,L,R,IL,BC,μ,i_d,i_r=Unregulated_Economy(;n_iter=10) #16 sec

#save
# output_ce=[DP,CT,PN,L,R,IL,BC,μ,i_d,i_r]
# save(string(pwd(), "/output_ce.jld"), "output_ce", output_ce)

plot(d_grid,DP, legend=:none)
plot(d_grid,CT, legend=:none)
plot(d_grid,μ , legend=:none)
plot(d_grid,L , legend=:none)
plot(d_grid,R , legend=:none)
plot(d_grid,IL, legend=:none)
plot(d_grid,BC, legend=:none)
plot(d_grid, [L[:,1] DP[:,1]], labels=["Loans" "Deposits"],legend=:topleft)
# plot(DIFF)

#check conditions (each of which should be zero)
#Eq. 14: Loans spread
sum( ((IL .- i_d)./(1 + i_d) .- Γ_l.(L,R)).*L .> 1e-8)

#Eq. 15: Reserves spread
sum( ((i_r - i_d)/(1 + i_d) .- Γ_r.(L,R)).*R .> 1e-8)

#Eq. 16: Balance sheet
sum( (L + R + Γ.(L,R) .- DP) .> 1e-8 ) 

#Eq. 22: Collateral constraint
sum(μ .< -1e-2 )

#Eq. 23: Collateral constraint
sum(L.> BC)
sum(DP.> BC)

#Eq. 24: Collateral constraint
sum(μ.*(κ.*(YT + PN.*YN) - L) .< -1e-2 )

# ---------------------------------------------------------------------------- #
#                        II. Social Planner Problem                            #
# ---------------------------------------------------------------------------- #
# i_d = copy(i_ast)
# i_r = 0

function Optimal_constrained(;n_iter=2,pfi=1)
    i_d  = copy(i_ast)

    # Planner needs a wider grid cause dp is not constraint
    dsp_grid = copy(d_grid)
    # dsp_grid = collect(range(d_min,1.2,length=nd))
    Dsp  = repeat(dsp_grid,1,ny)
    DDsp = repeat(dsp_grid,1,nd) #for computing expectations
    DPsp = repeat(dsp_grid,1,ny)

    # DPsp = copy(DP)
    CTsp = YT + D*(1+i_ast) - DPsp 
    Csp  = CES.(CTsp,YN)
    PNsp = ones(nd,ny)
    BCsp = ones(nd,ny)
    Vsp  = CRRA.(Csp)
    μsp  = zeros(nd,ny)
    λsp  = ∂u∂cT.(CTsp,YN)#zeros(nd,ny)

    Lsp  = zeros(nd,ny)
    Rsp  = zeros(nd,ny)
    ILsp = zeros(nd,ny)
    # IDsp = zeros(nd,ny) # para cuando control de capitales
    IRsp = zeros(nd,ny)

    #howard's improvement algorithm
    OPT  = zeros(nd,ny)
    util = zeros(nd,ny)
    Q    = sparse(zeros(nd*ny,nd*ny))

    l_grid = zeros(nd) 
    r_grid = zeros(nd) 

    #init value
    EVsp   = Vsp*Π' 

    n_iter = copy(n_iter)
    tol = 1e-4

    #method
    pfi = copy(pfi)
    @fastmath @inbounds begin 
        for iter in 1:n_iter
            DIFFsp = Float64[]

            if pfi == 1
                #### PFI
                DPsp_0 = copy(DPsp)
                CTsp_0 = copy(CTsp)
                PNsp_0 = copy(PNsp)
            else
                #### VFI
                Vsp_0= copy(Vsp) # VALUE FUNCTION ITERATION
            end


            policy_candidates = zeros(size(dsp_grid))
            # @fastmath @inbounds @simd for id in 1:nd
            for id in 1:nd    
                d=dsp_grid[id]
                # @fastmath @inbounds @simd for iy in 1:ny
                    for iy in 1:ny
                    yt = yT_grid[iy]
                    yn = yN_grid[iy]
                    pn = PNsp[id,iy]

                    # borr_const = κ*( (1-a)/a *( max(CTsp[id,iy],0)/yn ) ^(1/ξ) *yn + yt)
                    borr_const   = κ*(pn   *yn + yt) #ALTERNATIVE
                    BCsp[id,iy]  = borr_const
                    
                    #### #### #### #### #### Expectation #### #### #### #### #### 
                    # #OPTION1 (WORKING)
                    # for id_next in 1:length(dsp_grid) #@fastmath @inbounds @simd 
                    #     dp_g      = dsp_grid[id_next]  #guess on dp ##
                    #     # dp_g    = DPsp[id_next,iy]   #guess on dp ## dont haha

                    #     lt,rt      = get_L_R_SP(dp_g) # {lt, rt} using (7) and f.o.c that maximizes c^T                    
                    #     if lt ≤ borr_const
                    #         l_grid[id_next] = lt 
                    #         r_grid[id_next] = rt 

                    #     elseif lt > borr_const
                    #         lt = copy(borr_const)
                    #         find_rt(r_x)   = dp_g - (lt + r_x + Γ(lt,r_x)) 
                    #         rt  = fzero( find_rt , 0.1 )

                    #         l_grid[id_next] = lt 
                    #         r_grid[id_next] = rt 
                    #     end
                    #     ct  = max( yt - d*(1+i_ast) + dp_g - Γ.(lt,rt) - Γr.(rt), 1e-15)

                    #     c   = CES(ct,yn) # > 1e-15 ? CES(ct,yn) : 1e-15
                    #     policy_candidates[id_next] = CRRA(c) + β*EVsp[id_next,iy]
                    # end

                    #OPTION2 (ARREGLAR ESTO PARA QUE CORRA MAS JUERTE)
                    dp_vec = copy(dsp_grid)
                    lt,rt  = Get_L_R_SP(dp_vec) 
                    # donde lt ≤ borr_const:
                    l_grid[lt .≤ borr_const] = lt[lt .≤ borr_const] 
                    r_grid[lt .≤ borr_const] = rt[lt .≤ borr_const]
                    # donde lt > borr_const:
                    l_grid[lt .> borr_const] .= copy(borr_const)
                    #find optimal rt, given lt=constraint
                    n_binds      = sum(lt .> borr_const)            #preallocation
                    lt_bind      = l_grid[lt .> borr_const]    #preallocation
                    rt_bind      = r_grid[lt .> borr_const]    #preallocation
                    Threads.@threads for i in 1:n_binds
                        dp_vec_bind  = dp_vec[lt .> borr_const]
                        find_rt(r_x) = dp_vec_bind[i] - (lt_bind[i] + r_x + Γ(lt_bind[i],r_x)) 
                        rt_bind[i] = fzero( find_rt , 0.1 )
                    end
                    r_grid[lt .> borr_const] .= rt_bind
                    # end

                    ct_vec = yt .- d*(1+i_ast) .+ dp_vec - Γ.(l_grid,r_grid) - Γr.(r_grid)
                    c      = max.(CES.(ct_vec,yn),1e-15)
                    policy_candidates = CRRA.(c) + β*EVsp[:,iy]

                    #### #### #### #### #### #### #### #### #### #### #### #### #### #### #### 
                    Vsp[id,iy], pol_ind  = findmax(policy_candidates)
                    OPT[id,iy]           = pol_ind # PFI method
                    Lsp[id,iy]           = l_grid[pol_ind]
                    Rsp[id,iy]           = r_grid[pol_ind]
                end
            end
            ILsp = (1+i_d).*Γ_l.(Lsp,Rsp) .+ i_d
            IRsp = (1+i_d).*Γ_r.(Lsp,Rsp) .+ i_d

            DPsp = Lsp + Rsp + Γ.(Lsp,Rsp)
            CTsp = YT + DPsp - Dsp.*(1+i_ast)  - Γ.(Lsp,Rsp) - Γr.(Rsp)
            PNsp = ((1-a)/a) .*( max.(CTsp,0)./YN ) .^(1/ξ)
            
            if pfi == 1
                ####################################################
                ########## HOWARD'S IMPROVEMENT ALGORITHM ##########
                ####################################################
                Q = sparse(zeros(nd*ny,nd*ny))
                for jy=1:ny
                    yn  = yN_grid[jy]
                    ct  = CTsp[:,jy]
                    c   = CES.(ct,yn)

                    #update the value
                    util[:,jy] = CRRA.(c)

                    Q0 = sparse(zeros(nd,nd))
                    for jd=1:nd
                        Q0[jd,Int.(OPT[jd,jy])] = 1
                    end
                    Q[(jy-1)*nd+1:jy*nd,:] = kron(Π[jy,:]',Q0);
                end
                TVsp = (sparse(I(nd*ny)).-β*Q)\vec(util)
                Vsp = reshape(TVsp,nd,ny)

                dist=maximum( [ maximum(abs.(DPsp_0 .- DPsp)),
                                maximum(abs.(CTsp_0 .- CTsp)),
                                maximum(abs.(PNsp_0 .- PNsp))  ] ) # POLICY FUNCTION ITERATION
                        
                ####################################################
                ####################################################
                ####################################################
            else
                #### VALUE FUNCTION ITERATION
                ####check convergence
                dist=maximum( [ maximum(abs.(Vsp_0 .- Vsp))] )   
            end
            # ---------------------------------------------------------------------------- #
            #                                    UPDATE                                    #
            # ---------------------------------------------------------------------------- #
            DIFFsp=push!(DIFFsp,dist)

            println("Iter: $iter ; Conv: $dist")
            dist < tol ?  break : nothing

            EVsp = Vsp*Π'
        end
    end
    return DPsp,CTsp,PNsp,Lsp,Rsp,ILsp,IRsp,BCsp#,IDsp
end

@time DPsp,CTsp,PNsp,Lsp,Rsp,ILsp,IRsp,BCsp=Optimal_constrained(;n_iter=10,pfi=1)

#save
# output_sp=[DPsp,CTsp,PNsp,Lsp,Rsp,ILsp,IRsp,BCsp]
# save(string(pwd(), "/output_sp.jld"), "output_sp", output_sp)


# dsp_grid = collect(range(d_min,1.2,length=nd))
dsp_grid = copy(d_grid)

plot(dsp_grid,DPsp, legend=:none)
plot(dsp_grid,CTsp, legend=:none)
# plot(d_grid,μsp, legend=:none)
plot(dsp_grid,Lsp, legend=:none)
plot(dsp_grid,Rsp, legend=:none)
plot(dsp_grid,ILsp, legend=:none)
plot(dsp_grid,IRsp, legend=:none)
plot(dsp_grid,BCsp, legend=:none)
plot(dsp_grid, [Lsp[:,1] DPsp[:,1]], labels=["Loans" "Deposits"],legend=:topleft)




#Eq. 23: Collateral constraint
sum(Lsp.== BCsp)
sum(Lsp.> BCsp)
sum(DPsp.> BCsp)


# maximum(d_grid)
# maximum(DPsp)

# ---------------------------------------------------------------------------- #
#                              III. Simulate                                   #
# ---------------------------------------------------------------------------- #

function Simulate2!(n_sims,T,DP,L,R;cut=Int(floor(n_sims*0.2)))
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
    dp_sim   = zeros(length(y_shock))
    lt_sim   = zeros(length(y_shock))
    rt_sim   = zeros(length(y_shock))

    # bpsp_sim = zeros(length(y_shock))

    d_t   = mean(d_grid)
    # dp_sim[1] = mean(d_grid)
    # bsp_t = mean(b_grid)
    
    @showprogress for i in 2:length(y_shock)
        #Simulate y-shocks
        iy   = y_shock[i]
        id   = clamp(searchsortedfirst(d_grid,d_t),1,nd)
    
        #Competitive Equilibrium Simulations
        dp_sim[i]   = DP[id,iy] #/ (PN[id,iy]*YN[id,iy] + YT[id,iy]) *100
        lt_sim[i]   = L[id,iy]  #/ (PN[id,iy]*YN[id,iy] + YT[id,iy])*100
        rt_sim[i]   = R[id,iy]  #/ (PN[id,iy]*YN[id,iy] + YT[id,iy])*100

        #Social Planner Simulations
        # bpsp_sim[i]   = BPsp[ibsp,iy]# / (4*(ec.PNsp[ibsp,iy]*ec.YN[ibsp,iy] + ec.YT[ibsp,iy]))

        d_t   = copy(dp_sim[i])
        # bsp_t = copy(bpsp_sim[i] )
    end
    return dp_sim[cut:end], lt_sim[cut:end], rt_sim[cut:end]
end

pyplot(size = (900,800))
n_sims      = 1_000_000

# I. Simulate competitive
dp_sim  , lt_sim  , rt_sim   = Simulate2!(n_sims,Π,DP  ,L  ,R)
dpsp_sim, ltsp_sim, rtsp_sim = Simulate2!(n_sims,Π,DPsp,Lsp,Rsp)


fig1=density(d_grid,[dp_sim ,dpsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
        xaxis=(L"Deposits ($d_t$)"),color=[:red :blue],linestyle=[:dash :solid],
        label=["Competitive Equilibrium" "Optimal Constrained"],
        title= L"Density of next period debt (d$_{t+1}$)")
cd(dirname(@__FILE__))
savefig(fig1, "fig1.png")

# dsp_grid = collect(range(d_min,1.2,length=nd))

# density(d_grid,[dp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Deposits ($d_t$)"),color=[:red :blue],linestyle=[:dash :solid],
#         label=["Competitive Equilibrium"],
#         title= L"Density of next period debt (d$_{t+1}$)")
# density(dsp_grid,[dpsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Deposits ($d_t$)"),color=[:blue],linestyle=[:dash :solid],
#         label="Optimal Constrained",
#         title= L"Density of next period debt (d$_{t+1}$)")

r_den=density(d_grid,[rt_sim ,rtsp_sim], yaxis= ("Frequency", [0., 90]),lw=2.,
        xaxis=(L"Reserves ($r_t$)" , [0., 0.3]),color=[:red :blue],linestyle=[:dash :solid],
        label=["Competitive Equilibrium" "Optimal Constrained"],
        title= L"Reserves (r$_{t}$)");
# density(d_grid,[rt_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Reserves ($d_t$)"),color=[:red :blue],linestyle=[:dash :solid],
#         label="Competitive Equilibrium",
#         title= L"Reserves (r$_{t}$)")
# density!(dsp_grid,[rtsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Reserves ($d_t$)"),color=[:blue],linestyle=[:dash :solid],
#         label="Optimal Constrained",
#         title= L"Reserves (r$_{t}$)")

        
l_den=density(d_grid,[lt_sim,ltsp_sim ], yaxis= ("Frequency", [0., 90]),lw=2.,
        xaxis=(L"Loans ($l_t$)" , [0.7, 1.07]),color=[:red :blue],linestyle=[:dash :solid],
        label=["Competitive Equilibrium" "Optimal Constrained"],
        title= L"Loans (l$_{t}$)");

# histogram(d_grid,[lt_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Loans ($d_t$)"),color=[:red :blue],linestyle=[:dash :solid],
#         label="Competitive Equilibrium",
#         title= L"Loans (l$_{t}$)")
# density(dsp_grid,[ltsp_sim], yaxis= ("Frequency",[0. 60.]),lw=2.,
#         xaxis=(L"Loans ($d_t$)"),color=[:blue],linestyle=[:dash :solid],
#         label="Optimal Constrained",
#         title= L"Loans (l$_{t}$)")
       

fig2=plot(l_den,r_den)
cd(dirname(@__FILE__))
savefig(fig2, "fig2.png")

# histogram(d_grid,rt_sim)
# histogram(d_grid,lt_sim)