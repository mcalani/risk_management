#------------------------------------------------------------------------------#
# QuantEcon Cap.47: The Aiyagari Model: Heterogeneous Agent
#------------------------------------------------------------------------------#
using InstantiateFromURL
using LinearAlgebra, Statistics
using Parameters, Plots, QuantEcon

#------------------------------------------------------------------------------#
# Household problem:
#------------------
# max E Σ_{t=0}^{∞} β^{t}u(c_{t})
# st. a_{t+1} + c_{t} ≦ w * z_{t} + (1+r)a_{t},  c_{t} ≥ 0 and a_{t}≥ -B
#
#where w: wage rate ; z: markov process ; B: maximum borrowing
# state : s_{t} ≔ (a_{t}, z_{t})  (asset and shocks)
# action : the choce of a_{t+1}   (next period asset level)
#------------------------------------------------------------------------------#

#code that maps parameters for a household problem into the
# R and Q matrices needed to generate an instance of DiscreteDP.
Household= @with_kw (r= 0.01,
                     w=1.0,
                     σ=1.0,
                     β=0.96,
                     a_min=1e-10,
                     a_max=18.0,
                     a_size=200,
                     a_vals = range(a_min, a_max, length=a_size),
                     z_chain = MarkovChain([0.95 0.05; 0.1 0.9],[0.1;1.0]), # 2 states: 0.1 or 1
                     z_size= length(z_chain.state_values),
                     n= a_size*z_size,
                     s_vals= gridmake(a_vals, z_chain.state_values), #values
                     s_i_vals=gridmake(1:a_size, 1:z_size),          #indices
                     u= σ == 1 ? x -> log(x) : x -> (x^(1-σ)-1)/(1-σ),
                     R=setup_R!(fill(-Inf,n, a_size),a_vals, s_vals,r,w,u),
                     Q=setup_Q!(zeros(n, a_size,n),s_i_vals, z_chain));

#Transition: Q(s,a,s')
function setup_Q!(Q, s_i_vals, z_chain)
    for next_s_i in 1:size(Q,3)
        for a_i in 1: size(Q,2)
            for s_i in 1:size(Q,1)
                z_i = s_i_vals[s_i,2]
                next_z_i=s_i_vals[next_s_i,2]
                next_a_i = s_i_vals[next_s_i,1]
                if next_a_i==a_i  # a: a'=g(a,s)
                    Q[s_i, a_i, next_s_i]= z_chain.p[z_i, next_z_i]
                end
            end
        end
    end
    return Q
end

#Rewards: R(s,a)
function setup_R!(R, a_vals, s_vals, r, w,u)
    for new_a_i in 1:size(R,2)
        a_new = a_vals[new_a_i]
        for s_i in 1:size(R,1)
            a= s_vals[s_i,1]
            z= s_vals[s_i,2]
            c=w*z+(1+r)*a - a_new
            if c>0
                R[s_i, new_a_i]= u(c)
            end
        end
    end
    return R
end

#create an instance of Household
am= Household(a_max=20.0, r=0.03, w=0.956);

#Use the instance to build a discrete dynamic program
am_ddp= DiscreteDP(am.R, am.Q, am.β);

#solve using policy function iteration
results= solve(am_ddp, PFI);
results.Tv
#simplify names
@unpack z_size, a_size, n, a_vals = am;
z_vals= am.z_chain.state_values;

#Get all optimal actions across the set of a indices
# with z fixed in each column

a_star= reshape([a_vals[results.sigma[s_i]] for s_i in 1:n],
 a_size, z_size);

#The plot shows asset accumulation policies at
# different values of the exogenous state.
 labels = ["z= z_vals[1]", "z=z_vals[2]"];
plot(a_vals, a_star, label= labels, lw=2, alpha=0.6)
plot!(a_vals, a_vals, label="", color= :black, linestyle= :dash)
plot!(xlabel="current assests", ylabel= "next period assets", grid= false)

#------------------------------------------------------------------------------#
# Firm problem:
#--------------
# max A*K_{t}^{α} *N^{1-α} -(r-δ)K-wN
#
# from FOC:
# r=Aα(N/K)^{1-α} -δ
# w(r)= A(1-α)(Aα/(r+δ))^(α/(1-α))
#------------------------------------------------------------------------------#
# Now we want to calculate the equilibrium: to do this we need
# the intersection between interest rate and capital to obtain
# the equilibrium

const A=1;
const N=1;
const α=0.33;
const β=0.96;
const δ=0.05;

# w(r)
function r_to_w(r)
    return A*(1-α)*(A*α/(r+δ))^(α/(1-α))
end

# r(K): inverse of firm's demand of capital
function rd(K)
    return A * α*(N/K) ^ (1-α) - δ
end

function prices_to_capital_stock(am,r)
    #set up problem
    w=r_to_w(r)  #w(r)
    @unpack a_vals, s_vals, u = am
    setup_R!(am.R, a_vals, s_vals, r, w, u)

    aiyagari_ddp=DiscreteDP(am.R, am.Q, am.β)

    #compute the optimal policy
    results = solve(aiyagari_ddp, PFI)

    #compute the stationary Distribution
    stationary_probs=stationary_distributions(results.mc)[:,1][1]

    #return K
    return dot(am.s_vals[:,1], stationary_probs) #dot product
end

#Create an instance of Household
am= Household(β=β, a_max= 20.0);

#Create a grid of r values at wich to compute demand and supply of prices_to_capital_stock
r_vals= range(0.005, 0.04, length=20);

#compute supply of capital
k_vals= prices_to_capital_stock.(Ref(am), r_vals);

# plot against demand for capial by Firms
demand = rd.(k_vals);
labels= ["demand for capital " "supply of capital"];
plot(k_vals, [demand r_vals], label=labels, lw=2, alpha= 0.6)
plot!(xlabel= "capital", ylabel = "interest rate", xlim=(2,14), ylim=(0.0,0.1))
