#-------------------------------------------------------------------------------
#   Aiyagari 1994
#-------------------------------------------------------------------------------
using Parameters, Plots, QuantEcon
using LinearAlgebra, Statistics
using Optim

# ---------------------------------------------------------------------------- #
#                                  HOUSEHOLDS                                  #
# ---------------------------------------------------------------------------- #
#utility function
β = 0.96
μ = 3
u(x;σ=μ) = σ==1 ? log.(Float64(x)) : (Float64(x).^(1-μ)-1)./(1-μ)
r = 0.03
w = 0.94
# ---------------------------------------------------------------------------- #
#                                     FIRMS                                    #
# ---------------------------------------------------------------------------- #
#production function
α = 0.36
δ = 0.08
f(x) = x.^α .- (1-δ) .* x

# ---------------------------------------------------------------------------- #
#                                State Variables                               #
# ---------------------------------------------------------------------------- #
#labor endowment shocks
ρ = 0.6
σ_l = 0.2
mc = tauchen(2,ρ,σ_l);
π  = mc.p
states_z = ℯ.^collect(mc.state_values)
Nz   =length(states_z)

#grid for a
Na    = 200
a_min = 1e-10
a_max = 20
Δa    = (a_max-a_min)/(Na -1)
a_grid= collect(a_min:Δa:a_max)

#grid: states (a,z)
s_grid   = [vec(repeat(a_grid,Nz,1)) kron(states_z, repeat([1], Na))] # Ns x Number of state variables
s_grid_i = [vec(repeat(1:1:Na,Nz,1)) kron(1:1:Nz, repeat([1], Na))]



# ---------------------------------------------------------------------------- #
#                           COMPUTE HOUSEHOLD PROBLEM                          #
# ---------------------------------------------------------------------------- #


function Households(;tol=1e-3,w=w,r=r,print=true)
    # Preallocation
    R       = fill(-Inf, Na*Nz,Na)  # reward function: (state, action)
    Q       = zeros(Na*Nz,Na,Na*Nz) # transition probability function: (state, action, next_state)

    #Rewards
    function compute_R!(R, s_grid,a_grid) #para cada estado, retorna todos los rewards posibles
        #loop over control
        for i in 1:size(R,1) #state
            a,z = s_grid[i,:]
            for j in 1:size(R,2) #action 
                ap = a_grid[j]
                c  = w*z+(1+r)*a - ap
                if c > 0
                    R[i,j] = u(c)
                end
            end
        end
        return R
    end
    
    #Transition
    function compute_Q!(Q, s_grid_i, Π)
        for si in 1:size(Q,1) #state  s:=(a,z)
            ~,zi =  s_grid_i[si,:] #prob t
            for  σi in 1:size(Q,2) #action a' = σ(a)
                for sii in 1:size(Q,3) #next_state s':=(a',z')
                    aii,zii =  s_grid_i[sii,:] #prob t+1 of state aii
                    σi == aii ? Q[si,aii,sii] = Π[zi,zii] : nothing
                end
            end
        end
        return Q
    end
    
    function max_n_ind!(vals,TV,σ_greedy)
        nrow,ncol=size(vals) #(s,a)
        # σ      = zeros(nrow)
        # TV   = zeros(nrow)
        σ_i    = zeros(nrow)
        TV_i = zeros(nrow)
        for ir in 1:nrow
            TV[ir],TV_i[ir] = findmax(vals[ir,:])
        end 
        σ_greedy[:] = s_grid[Int.(TV_i),1]
        TV,σ_greedy
    end
    
    function BellmanOp(R,Q,β,v,Tv,σ,w,r)
        vals = R + β * (Q * v)
        # s_wise_max!(vals, Tv, σ)
        max_n_ind!(vals, Tv, σ)
    
        Tv, σ
    end
    

    R = compute_R!(R, s_grid,a_grid)
    Q = compute_Q!(Q,s_grid_i,π)
    
    v  = zeros(size(R,1))
    Tv = zeros(size(R,1))
    σ  = zeros(size(R,1))

    dif = 100
    # for i in 1:250
    while dif > tol
        # println(i)
        Tv,σ = BellmanOp(R,Q,β,v,Tv,σ,w,r)

        dif = norm(v.-Tv)
        print== true ? println(dif) : nothing
    
        v[:]=Tv[:]
    end    
    return v,σ
end

V,σ = Households(w=w,r=r;tol=1e-3)

pyplot()
pl_sigma = plot(a_grid,[σ[1:Na] σ[(Na+1):end]])
pl_value = plot(a_grid,[V[1:Na] V[(Na+1):end]])

# plot([σ[1:200] σ[201:end]])





# ---------------------------------------------------------------------------- #
#                                     FIRMS                                    #
# ---------------------------------------------------------------------------- #

# Firms' parameters
A = 1
N = 1

#production function
α = 0.36
δ = 0.08
# f(x) = x.^α .- (1-δ) .* x


#f.o.c : 
r_demand(K) = A * α * (N / K) ^ (1 - α) - δ
new_w(r) = A * (1 - α) * (A * α / (r + δ)) ^ (α / (1 - α))


# ---------------------------------------------------------------------------- #
#                              Compute Equilibrium                             #
# ---------------------------------------------------------------------------- #
Π_ss = stationary_distributions(MarkovChain(π))[1]

# eigen_ind = findall(eigen(π).values .== 1)
# Π_ss = vec(eigen(π).vectors[:,eigen_ind]./sum(eigen(π).vectors[:,eigen_ind]))



r_grid = collect(range(0.005, 0.08, length = 30))
K      = zeros(length(r_grid))

ai = 1 #fixing a level of capital
# V_iter = copy(V)

@time for (i,ri) in enumerate(r_grid)
  println("iter=$i ; IntRate =$ri")
  wi   = new_w(ri)
  ~,σ  = Households(w=wi, r=ri;tol=1e-2,print=false)
  K[i] = dot(Matrix([σ[ai] σ[Na+ai-1]] ), Π_ss) #supply of capital
end


pyplot()
pl_eq = plot(K,[r_demand.(K) r_grid])

summ = plot(pl_eq,pl_sigma,pl_value)
cd(dirname(@__FILE__))
savefig(pl_eq,"froot_hast.png")