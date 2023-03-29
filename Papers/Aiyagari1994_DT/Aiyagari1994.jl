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
μ = 1
u(x;σ=μ) = σ==1 ? log.(Float64(x)) : (Float64(x).^(1-μ)-1)./(1-μ)
r = 0.01
w = 1#0.94
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
# ρ = 0.6
# σ_l = 0.2
# mc = tauchen(2,ρ,σ_l);
# π  = mc.p
# states_z = ℯ.^collect(mc.state_values)
# Nz   =length(states_z)
π  = [0.95 0.05; 0.1 0.9]
states_z = [0.1;1.0]
Nz   =length(states_z)

#grid for a
Na    = 800
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
function compute_R!(R, s_grid,a_grid;r=r,w=w) #para cada estado, retorna todos los rewards posibles
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

function BellmanOp(R,Q,β,v,Tv,σ)
    vals = R + β * (Q * v)
    # s_wise_max!(vals, Tv, σ)
    max_n_ind!(vals, Tv, σ)

    Tv, σ
end


function Households(;tol=1e-3,w=w,r=r,print=true)
    # Preallocation
    R       = fill(-Inf, Na*Nz,Na)  # reward function: (state, action)
    Q       = zeros(Na*Nz,Na,Na*Nz) # transition probability function: (state, action, next_state)
    #Rewards    

    R = compute_R!(R, s_grid,a_grid;r=r,w=w)
    Q = compute_Q!(Q,s_grid_i,π)
    
    v  = zeros(size(R,1))
    Tv = zeros(size(R,1))
    σ  = zeros(size(R,1))

    dif = 100
    # for i in 1:250
    while dif > tol
        # println(i)
        Tv,σ = BellmanOp(R,Q,β,v,Tv,σ)

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




# ---------------------------------------------------------------------------- #
#                             Stationary Equilibrum                            #
# ---------------------------------------------------------------------------- #

function stationary_equilibrum(;r=r,w=w,σ=σ)
    R = fill(-Inf, Na*Nz,Na)  # reward function: (state, action)
    Q = zeros(Na*Nz,Na,Na*Nz) # transition probability function: (state, action, next_state)

    R = compute_R!(R, s_grid,a_grid;r=r,w=w)
    Q = compute_Q!(Q, s_grid_i, π)

    #recover optimal indices
    σ1 = σ[1:Na] 
    σ2 = σ[Na+1:end] 
    σ_i = zeros(size(σ))
    for i in 1:size(σ_i,1)
        if i ≤ Na
            σ_i[i] = Int.(findfirst(σ1[i] .== a_grid ))
        else
            σ_i[i] = Int.(findfirst(σ2[i-Na] .== a_grid ))
        end
    end

    #recover rewards from R
    R_vec = zeros(size(σ_i,1))
    for i in 1:size(σ_i,1)
            R_vec[i] = R[i,Int.(σ_i[i])]
    end
    # plot(R_vec[1:Na])
    # plot!(R_vec[(Na+1):end])


    #recover transition matriz from Q
    Q_mat = zeros(size(σ_i,1),size(σ_i,1))
    for i in 1:size(Q_mat,1)
        for j in 1:size(Q_mat,2)
            Q_mat[i,j] = Q[i,Int.(σ_i[i]),j]
        end
    end

    # Q_mat

    #stationary distribution

    function sdist(T ; method="Eigenvector")
        function findnearest(v::AbstractArray, c::Real)
            findmin(abs.(v .- c ))[2]
        end;
        rows =[ sum(T[i,:])  for i in 1:size(T,1)]
        @assert size(T,1) == size(T,2)      #squared matrix
        @assert sum(rows .≈ 1) == size(T,1) #rows must sum 1
        if method == "Eigenvector"
            λ_pos=findnearest(real(eigvals(T)) , 1.)
            #Method1: eigenvalues
            λ     = real(eigvals(T))[λ_pos]
            @assert λ≈1
            Reig  = eigvecs(T) #right eigenvalues
            Leig  = inv(Reig)  #left  eigenvalues (stationary distribution)
            v     = vec(real.(Leig[end,:]))
            #scaling v
            edist1=vec( v'T./sum(v'T)  )
            return edist1
        elseif method == "Iteration"
            v=ones(size(T,1))./size(T,1)
            count = 0
            while count < 100_000
                v1 = vec(v'T) # left eigenvalue associated with eigenvalue=1
                dif=norm(v1-v)
                #println(dif)
                v=copy(v1)
                if vec(v'T) ≈ v
                    break
                end
                count += 1
            end
            #scaling v
            edist2=v./sum(v)
            return edist2 
        end
    end

    Q_ss = sdist(Q_mat; method="Iteration")
    return Q_ss
end

Q_ss = stationary_equilibrum(;σ=σ,r=r,w=w)

# σ' * Q_ss
# dot(σ, Q_ss)
plot(a_grid,Q_ss[1:Na])
plot!(a_grid,Q_ss[Na+1:end])



#check with library
# Π_ss = stationary_distributions(MarkovChain(Q_mat))[1]
# plot(a_grid,Π_ss[1:Na])
# plot!(a_grid,Π_ss[Na+1:end])

# ---------------------------------------------------------------------------- #
#                                     FIRMS                                    #
# ---------------------------------------------------------------------------- #

# Firms' parameters
A = 1
N = 1

#production function
α = 0.33 #0.36
δ = 0.05 #0.08
# f(x) = x.^α .- (1-δ) .* x


#f.o.c : 
r_demand(K) = A * α * (N / K) ^ (1 - α) - δ
new_w(r) = A * (1 - α) * (A * α / (r + δ)) ^ (α / (1 - α))


# ---------------------------------------------------------------------------- #
#                              Market Equilibrium                              #
# ---------------------------------------------------------------------------- #

r_grid = collect(range(0.005, 0.06, length = 15))
K      = zeros(length(r_grid))

# V_iter = copy(V)

@time for (i,ri) in enumerate(r_grid)
  println("iter=$i ; IntRate =$ri")
  wi   = new_w(ri)
  ~,σi  = Households(w=wi, r=ri;tol=1e-2,print=false)
  Q_ssi = stationary_equilibrum(;σ=σi,w=wi, r=ri)
  K[i] = σi' * Q_ssi #supply of capital
end

# plot(K)

pyplot()
pl_eq = plot(K,[r_demand.(K) r_grid])

summ = plot(pl_eq,pl_sigma,pl_value)
cd(dirname(@__FILE__))
# savefig(summ,"fig.png")


# ---------------------------------------------------------------------------- #
#                                 Distribution                                 #
# ---------------------------------------------------------------------------- #
r_i = 0.03235
K_d(r) = N* ((A*α)/(r+δ))^(1/(1-α))
Kd = K_d(r_i)
dif=100
while dif > 1e-3
  println("Eq. Int.Rate =$r_i")
  w_i    = new_w(r_i)

  #compute supply of capital
  Vi,σi  = Households(w=w_i, r=r_i;tol=1e-2,print=false)
  Q_ssi  = stationary_equilibrum(;σ=σi,w=w_i, r=r_i)
  Ks     = σi' * Q_ssi 
  r_s    = r_demand(Ks) #tasa consistente S

  #compute demand of capital
  Kd     = K_d(r_i)

  dif = abs(r_i-r_s) 
  println("convergence= $dif")
  r_i = (r_i + r_s)/2
end

r_i
w_i    = new_w(r_i)
~,σ  = Households(w=w_i, r=r_i;tol=1e-2,print=true)
Q_ss = stationary_equilibrum(;σ=σ,w=w_i, r=r_i)

plot(a_grid,Q_ss[1:Na])
plot!(a_grid,Q_ss[Na+1:end])