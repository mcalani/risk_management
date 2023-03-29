# ---------------------------------------------------------------------------- #
#                                 AIYAGARI 1994                                #
# ---------------------------------------------------------------------------- #
using LinearAlgebra, Statistics
using Parameters, Plots, QuantEcon

#Initial parameters
const r = 0.03
const w = .956
const σ = 1.0
const β = 0.96

#Grid
const a_min  = 1e-10
const a_max  = 20.0
const Na     = 200
const a_grid = collect(range(a_min, a_max, length=Na))

#Stochastic process
const z_chain      = MarkovChain([0.95 0.05; 0.1 0.9],[0.1;1.0]) # 2 states: 0.1 or 1
const z_grid       = z_chain.state_values
const Π            = z_chain.p
const Nz           = length(z_chain.state_values)
const Ns           = Na*Nz
const s_grid       = gridmake(a_grid, z_grid) 
const si_grid      = gridmake(1:Na, 1:Nz)          #indices
const u(c; σ=σ)    = σ == 1 ? log(c) : (c^(1-σ)-1)/(1-σ)

const Π_ss = stationary_distributions(z_chain)[1]

# ---------------------------------------------------------------------------- #
#                                  HOUSEHOLDS                                  #
# ---------------------------------------------------------------------------- #


function Households(;r=r,w=w,a_grid=a_grid,z_grid=z_grid,tol=1e-2,V_init=false,print=false)
  #parameters

  #initial values
  aux_positive(x) = x > 0 ? x : 1e-25
  a_guess  = repeat(a_grid,1,Nz)
  c_guess  = w*z_grid' .+ (1+r).*a_grid - a_guess
  U        = u.(aux_positive.(c_guess))
  V_init == false ? V0 = copy(U) : V0 = copy(V_init)
  
  #first iteration
  V_max_i = zeros(Na,Nz) #V_max index
  V_max   = zeros(Na,Nz)
  V_val   = zeros(Na,Nz,Na)

  diff = norm(V0 .- V_max)
  #states -> (A,Z,A')
  # for zzz in 1:n_iter
  while diff > tol #1e-3
    EV = Π * V0' #[Int.(V_max_i)]
    for (i,a) in enumerate(a_grid)
      for (ii,ap) in enumerate(a_grid) 
        C = w.*z_grid .+ (1+r)*a .- ap
        C[C.<0] .= 0
        V_val[i,:,ii] = u.(C) .+ β .* EV[ii] 
      end
      for zi in 1:length(z_grid)
        V_max[i,zi],V_max_i[i,zi] = findmax(V_val[i,zi,:])
      end
      # V_max[i,:]   .= vec(findmax(V_val[i,:,:],dims=2)[1])
      # V_max_i       = push!(V_max_i,findmax(V_val[i,:,:],dims=2)[2])
    end
    diff = norm(V0 .- V_max)
    print== true ? println(diff) : nothing 
    V0 = copy(V_max)
  end
  return V_max , V_max_i   
end

V , Vi= Households(;r=r,w=w,a_grid=a_grid,z_grid=z_grid,tol=1e-3,print=true) #,V_init=V

pyplot()
plot(a_grid,V, size=(900,600),legend=:bottomright, title="Value Function")
plot(a_grid,a_grid[Int.(Vi)],legend=:bottomright, title="Policy Function")
plot!(a_grid,a_grid,color=:black,linestyle=:dash,legend=:bottomright, title="Policy Function")



# ---------------------------------------------------------------------------- #
#                                     FIRMS                                    #
# ---------------------------------------------------------------------------- #

# Firms' parameters
const A = 1
const N = 1
const α = 0.33
# const β = 0.96
const δ = 0.05



#f.o.c : 
r_demand(K) = A * α * (N / K) ^ (1 - α) - δ
new_w(r) = A * (1 - α) * (A * α / (r + δ)) ^ (α / (1 - α))


#Iteration
r_grid = collect(range(0.02, 0.06, length = 30))
K      = zeros(length(r_grid))

ai = 120 #fixing a level of capital
# V_iter = copy(V)

@time for (i,ri) in enumerate(r_grid)
  println("iter=$i ; IntRate =$ri")
  wi = new_w(ri)
  ~ , Vi_iter = Households(;w=wi, r=ri,tol=1e-5,print=false)
  K[i] = dot(a_grid[Int.(Vi_iter)][ai,:] , Π_ss) #supply of capital
end

plot(K,[r_demand.(K) r_grid])