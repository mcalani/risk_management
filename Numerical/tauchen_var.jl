# Author: Marco Piña
# Last revised: Jan, 2021



# ---------------------------------------------------------------------------- #
#                               Discretize VAR(1)                              #
# ---------------------------------------------------------------------------- #
# y_t = A*y_t-1 + ε_t ,  ε_t ∼ N(0,Σ)

using LinearAlgebra
using Statistics
using Distributions

function Tauchen_VAR(A,Σ,knots;m=3)
     NVAR  = size(Σ,1) # number of variables  
     global nodes = copy(knots)
     #step 1: LDLt decomposition
     # L*Σ*L' = Λ
     Λ=diagm(eigen(Σ).values)
     Q=eigen(Σ).vectors
     #Q*Λ*Q' ≈ Σ #check

     #step2 : Transform the model
     # y'_t = B*y'_t-1 + ε'_t
     # where
     # y' = Q'*y
     # B  = Q'*A*Q 
     # ε' = Q'ε
     # Variance(ε') = Λ
     # Variance(y') = B*Σyp*B' + Q'*Σε*Q = B*Σyp*B' + Λ


     #Variance of yp
     B = Q' * A * Q
     K= kron(B,B)
     Σ_yp=reshape(inv(I-K)*vec(Λ),NVAR,NVAR)

     YP_grid = zeros(nodes,NVAR)
     step    = zeros(NVAR)
     for i in 1:NVAR
          σy = sqrt(Σ_yp[i,i])
          YP_grid[:,i] = collect(range(-m*σy,m*σy,length=nodes))
          step[i] = diff(YP_grid[:,i])[1]        
     end

     #step 3: make grid of all possible states
     states_yp = zeros(nodes^NVAR,NVAR)
     states_yp[:,1] = kron(YP_grid[:,1],ones(nodes)*(NVAR-1))
     for i in 2:NVAR
          states_yp[:,i]=repeat(YP_grid[:,i],nodes^(NVAR-i+1))
     end

     #step 4: Compute transition matrices
     Φ(x;μ=0.,σ=1.)=cdf(Normal(μ,σ),x)

     global T_yp = [] #preallocation for each y-process
     for i in 1:NVAR
          @eval $(Symbol("T$(i)_yp")) = zeros(nodes,nodes)
          @eval T_yp = push!(T_yp,$(Symbol("T$(i)_yp")))
     end

     for i in 1:NVAR #Process
          T       = T_yp[i]      #transition
          yp_grid = YP_grid[:,i] #grid
          d       = step[i]      #step
          ρ       = B[i,:]       #row of coefficients
          σ_εp    = sqrt(Λ[i,i]) #σ of error
          for s in 1:nodes  #current state of vector yp_t
               y_s = states_yp[s,:]
               for sp in 1:nodes  #next state of scalar yp_t+1
                    y_sp = yp_grid[sp]
                    if sp == 1
                         T[s,sp] = Φ( (y_sp - ρ'*y_s  + d/2)/sqrt(σ_εp) )
                    elseif sp==nodes
                         T[s,sp] = 1 - Φ( (y_sp - ρ'*y_s - d/2)/sqrt(σ_εp) )
                    else
                         T[s,sp] = Φ( (y_sp - ρ'*y_s + d/2)/sqrt(σ_εp) ) -
                                   Φ( (y_sp - ρ'*y_s  - d/2)/sqrt(σ_εp) )
                    end
               end
          end
          T_yp[i]=T
     end
     #sum.(T_yp) #check

     T = copy(T_yp[1])
     for i in 2:NVAR
          T = kron(T,T_yp[i])
     end 

     #Step 5: recover states for y such thath yp = Q'y => (Q')^-1 * yp = y
     States_y = zeros(nodes^NVAR,NVAR)
     for j in 1:nodes^NVAR
               States_y[j,:] = inv(Q')*states_yp[j,:]
     end

     return T, States_y
end


# Example: Bianchi2011 Calibration
# A = [ 0.901009 0.495189 ;
#      -0.453210 0.225130]
# Σ = [0.00219105027284 0.00162283372642;
#      0.00162283372642 0.00169597630831]
# nodes = 4  # nodes
# T,States=Tauchen_VAR(A,Σ,nodes)
# T
# exp.(States)
# using Plots
# plot(exp.(States))
# plot(T,legend=:none)
# sum(T) ≈ size(T,1) #True

# #Example 2:
# A = [ 0.745676 0.108413 ;
#       0.029050 0.739263]
# Σ = [0.000287 0.000039;
#      0.000039 0.000107]
# #nodes = 4  # nodes
# #nodes=Any

# T,States=Tauchen_VAR(A,Σ,4)



# using Plots
# plot(exp.(States))
# plot(T,legend=:none)
# sum(T) ≈ size(T,1) #True

# exp.(States[:,1]) .+ .03
# exp.(States[:,1])
