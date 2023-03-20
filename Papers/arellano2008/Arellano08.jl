# ---------------------------------------------------------------------------- #
#                                 Arellano 2008                                #
# ---------------------------------------------------------------------------- #
using Parameters, Plots, QuantEcon
using LinearAlgebra, Statistics
using Optim
using Interpolations
using ProgressMeter
#using BenchmarkTools

#Parameters
β   = 0.953
γ   = 2.0
r   = 0.017
ρ   = 0.945
η   = 0.025
θ   = 0.282 #reentering international credit market probability
ny  = 10#21
nB  = 100 #251

#Discretization
Π        = tauchen(ny, ρ, η).p
y_grid   = exp.(tauchen(ny, ρ, η).state_values)
ydefgrid = min.(.969 * mean(y_grid), y_grid)
B_grid   = collect(range(-.4, .4, length = nB))

#Initial values
vf = zeros(nB, ny)
vd = zeros(1, ny)
vc = zeros(nB, ny)
policy = zeros(nB, ny)
q = ones(nB, ny) .* (1 / (1 + r))
defprob = zeros(nB, ny)

function u(c;γ=2)
    if c > 0.0
        γ==1 ? log(c) : (c^(1-γ) -1)/(1-γ)
    else
        -Inf
    end
end

#Problem:
# E{ Σ_{t=0}^{∞} β^t * u(c_t)  }
# c = y + B - q(B',y)B'
# B ≥ -Z 

#Equilibrium: 
#Default:
# V_d = u(y - qB' + B) + β*E{  θ*V(0,y') + (1-θ)V_d } 

#Paying:
# V_c = max{  u(y - qB' + B) + β*E{V}  β*E{ V } }

#V = max{ V_c, V_d }

#Default prob:
# δ := ∫ 1{V_c < V_d}p(y,y')dy'
# q= (1-δ)/(1+r)

# ---------------------------------------------------------------------------- #
#                               Bellman Operator                               #
# ---------------------------------------------------------------------------- #


function T_operator(vd,vc,vf,q)
    #Preallocations
    Tvd = zeros(1, ny)
    Tvc = zeros(nB, ny)
    Tvf = zeros(nB, ny)
    defs= zeros(nB, ny)
    Bpol= zeros(nB, ny)
    policy_candidates = zeros(nB)
   
    #Expectation
    EVd = vd* Π'
    EVc = vc* Π'
    EVf = vf* Π'

    #Find index where B=0
    zero_ind = searchsortedfirst(B_grid,0.)
    for iy in 1:ny
        y    = y_grid[iy]
        ydef = ydefgrid[iy]

        #Default Value
        Tvd[1,iy] = u(ydef) + β*(θ*EVc[zero_ind, iy] + (1-θ)*EVd[1, iy])

        for ib in 1:nB
            B = B_grid[ib]
            #Paying debt Value
            for ibp in 1:nB #try every possible next value and keep the maximum
                Bp = B_grid[ibp]
                c  = max(y + B - q[ibp,iy]*Bp, 1e-14)
                policy_candidates[ibp] = u(c) + β*EVf[ibp,iy]
            end
            Tvc[ib,iy],Bpol[ib,iy]=findmax(policy_candidates)

            #Value Function
            Tvf[ib,iy]   = max(Tvc[ib,iy],Tvd[1,iy])
        end
    end
    return Tvd, Tvc, Tvf, Bpol
end

# # #test
# @time Tvd, Tvc, Tvf, Bpol = T_operator(vd,vc,vf,q)
# #Tvd, Tvc, Tvf,Bp, defs = T_operator(vd,vc,vf,q)
# plot(B_grid,Tvc)
# plot(B_grid,Tvf)
# plot(Tvd')

# ---------------------------------------------------------------------------- #
#                                      VFI                                     #
# ---------------------------------------------------------------------------- #

function VFI(vd,vc,vf,q,T_operator;show_every=10)
    iter  = 0
    N_iter= 800
    conv  = 100
    tol   = 1e-8
    δ     = zeros(nB, ny)
    Bpol  = zeros(nB, ny)
    for i in 1:N_iter
        iter += 1

        #Bellman operator
        vd1, vc1, vf1, Bpol= T_operator(vd,vc,vf,q)        
        #converce criterion
        conv  = norm( vf1 - vf )

        if iter % show_every == 0
            println("Diff = $(round(conv,digits=10)) ; Iter= $iter")
        end
        if conv < tol
            break
            println("Convergence achieved :D !")
        end

        #update bond price
        default_states = repeat(vd1,nB) .> vc1
        δ = default_states*Π'
        q = (1 .- δ) / (1 + r)

        #Update values
        vd = copy(vd1)
        vc = copy(vc1)
        vf = copy(vf1)
    end
    return vd,vc,vf, Bpol, q, δ
end

#test
@time vd,vc,vf, Bpol, q, δ = VFI(vd,vc,vf,q,T_operator)


# ---------------------------------------------------------------------------- #
#                                     Plot                                     #
# ---------------------------------------------------------------------------- #

# create "Y High" and "Y Low" values as 5% devs from mean
high, low = 1.05 * mean(y_grid), 0.95 * mean(y_grid)
iy_high, iy_low = map(x -> searchsortedfirst(y_grid, x), (high, low))

plot(B_grid,[ q[:,iy_high] q[:,iy_low] ],xaxis=[-0.4,0.])


heatmap(B_grid[1:end-1],y_grid[2:end], xaxis=[-0.4, 0.05],
    reshape(clamp.(vec(δ[1:end - 1, 1:end - 1]), 0, 1), nB-1, ny-1)')
plot!(xlabel = "B'", ylabel = "y", title = "Probability of default",
    legend = :topleft)

