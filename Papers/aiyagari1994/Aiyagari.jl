
## #
#   Aiyagari (1994)
#
## #
using LinearAlgebra, Statistics
using Parameters, Plots, QuantEcon

#Initial parameters
r = 0.01
w = 1.0
σ = 1.0
β = 0.96

#Grid
a_min  = 1e-10
a_max  = 18.0
a_size = 200
a_vals = range(a_min, a_max, length=a_size)

#Stochastic process
z_chain = MarkovChain([0.95 0.05; 0.1 0.9],[0.1;1.0]) # 2 states: 0.1 or 1
z_states       = z_chain.state_values
z_transition = z_chain.p
z_size  = length(z_chain.state_values)
n       = a_size*z_size
s_vals  = gridmake(a_vals, z_chain.state_values) #values
s_i_vals= gridmake(1:a_size, 1:z_size)          #indices
u = σ == 1 ? x -> log(x) : x -> (x^(1-σ)-1)/(1-σ)

## ##
# Household problem:
# max E Σ_{t=0}^{∞} β^{t}u(c_{t})
# st. a_{t+1} + c_{t} ≦ w * z_{t} + (1+r)a_{t},  c_{t} ≥ 0 and a_{t}≥ -B
#
# where w: wage rate ; z: markov process ; B: maximum borrowing
# state : s_{t} ≔ (a_{t}, z_{t})  (asset and shocks)
# action : the choce of a_{t+1}   (next period asset level)
## ##

#initial values
V0 = ones(a_size,z_size);
C  = ones(a_size,z_size);
C[:,1]  = w*z_states[1] .+ (1+r).*a_vals; #bad state
C[:,2]  = w*z_states[2] .+ (1+r).*a_vals; #good state
U = u.(C)


size(z_transition)[1]
function expec!(z_prob,z_states,V)
  n = size(z_prob)[1]
  expec=0
  for i in 1:n
    expec = expec + z_prob[1,i]*V



  end
