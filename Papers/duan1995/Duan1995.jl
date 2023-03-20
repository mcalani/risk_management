# ---------------------------------------------------------------------------- #
#                                   Duan 1995                                  #
# ---------------------------------------------------------------------------- #
using LinearAlgebra    #norm function (compute norm of a matrix)
using Plots            #plots
using Random         #seed
using Distributions  
using ProgressMeter  #to show progress



# ---------------------------------------------------------------------------- #
#                       Black-Scholes-Merton (P measure)                       #
# ---------------------------------------------------------------------------- #
function BS(;ω	= 0.00001524,
             α	= 0.1883,
             β  = 0.7162 )
    ρ_P  = α + β
    σ2_P = ω/(1-ρ_P)
    σA   = sqrt(365*σ2_P) #annual vol
            
    t30  = σA*sqrt(30/365)
    t90  = σA*sqrt(90/365)
    t180 = σA*sqrt(180/365)

    B = zeros(9,4)
    cnorm(x) = cdf(Normal(0,1),x)

    for  k in 1:9 # each k is a different strike price
        x         = 0.75 + 0.05*k
        B[k,1]    = x

        #T=30
        d1        = (log(x) + 0.5*t30*t30)/t30
        d2        = d1 - t30
        B[k,2]   = 10000*(x*cnorm(d1) - cnorm(d2))   

        #T=90
        d1 = (log(x) + 0.5*t90*t90)/t90
        d2 = d1 - t90
        B[k,3]   = 10000*(x*cnorm(d1) - cnorm(d2))   

        #T=180
        d1 = (log(x) + 0.5*t180*t180)/t180
        d2 = d1 - t180
        B[k,4]   = 10000*(x*cnorm(d1) - cnorm(d2))   
    end
    return B
end


# ---------------------------------------------------------------------------- #
#                                   Duan 1995                                  #
# ---------------------------------------------------------------------------- #
function Duan1995(; ω	= 0.00001524,
                    α	= 0.1883,
                    β   = 0.7162,
                    λ	= 0.007452,
                    n_sims=30_000,
                    seed=123)
    # ---------------------------------------------------------------------------- #
    #                            Monte Carlo (Q measure)                           #
    # ---------------------------------------------------------------------------- #

    # Parámetros adicionales
    ρ_P  = α + β
    σ2_P = ω/(1-ρ_P)
    σA   = sqrt(365*σ2_P) #annual vol
    ρ_Q  = α*(1+λ^2) + β
    σ2_Q = ω/(1-ρ_Q)
    σAQ  = sqrt(365*σ2_Q)

    #Sims (a,b,c == columnds of table 4.1 in Duan(1995))
    #Preallocation (vol, ret, prices)
    h_a = zeros(180,n_sims)
    h_b = zeros(180,n_sims)
    h_c = zeros(180,n_sims)

    r_a = zeros(180,n_sims)
    r_b = zeros(180,n_sims)
    r_c = zeros(180,n_sims)

    S30a	= zeros(n_sims)
    S30b	= zeros(n_sims)
    S30c	= zeros(n_sims)
    S90a	= zeros(n_sims)
    S90b	= zeros(n_sims)
    S90c	= zeros(n_sims)
    S180a	= zeros(n_sims)
    S180b	= zeros(n_sims)
    S180c	= zeros(n_sims)

    C30a   = zeros(n_sims,9)
    C30b   = zeros(n_sims,9)
    C30c   = zeros(n_sims,9)
    C90a   = zeros(n_sims,9)
    C90b   = zeros(n_sims,9)
    C90c   = zeros(n_sims,9)
    C180a  = zeros(n_sims,9)
    C180b  = zeros(n_sims,9)         
    C180c  = zeros(n_sims,9)
    #initial values:
    h_a[1,:] .= (0.8)^2*σ2_Q
    h_b[1,:] .= (1.00)^2*σ2_Q
    h_c[1,:] .= (1.2)^2*σ2_Q

    #set seed
    Random.seed!(seed)
    ξ = rand(Normal(0,1),180,n_sims)

    for t in 2:180
        # volatilities (180days)
        h_a[t,:] = ω .+ α.*h_a[t-1,:].*((ξ[t-1,:].-λ).^2) .+ β.*h_a[t-1,:] 
        h_b[t,:] = ω .+ α.*h_b[t-1,:].*((ξ[t-1,:].-λ).^2) .+ β.*h_b[t-1,:]
        h_c[t,:] = ω .+ α.*h_c[t-1,:].*((ξ[t-1,:].-λ).^2) .+ β.*h_c[t-1,:] 
    end
    # plot(h_c[:,1])
    # plot!(h_b[:,1])
    # plot!(h_a[:,1])

    # returns (180days)
    r_a = -0.5 .*h_a .+ sqrt.(h_a).*ξ
    r_b = -0.5 .*h_b .+ sqrt.(h_b).*ξ
    r_c = -0.5 .*h_c .+ sqrt.(h_c).*ξ

    # Stock Prices (30, 90 and 180 days)
    S30a  = exp.(sum(r_a[1:30,:],dims=1))
    S30b  = exp.(sum(r_b[1:30,:],dims=1))
    S30c  = exp.(sum(r_c[1:30,:],dims=1))  

    S90a  = exp.(sum(r_a[1:90,:],dims=1))
    S90b  = exp.(sum(r_b[1:90,:],dims=1))
    S90c  = exp.(sum(r_c[1:90,:],dims=1))

    S180a = exp.(sum(r_a[1:end,:],dims=1))
    S180b = exp.(sum(r_b[1:end,:],dims=1))
    S180c = exp.(sum(r_c[1:end,:],dims=1))

    # Call options (30, 90 and 180 days)
    for k in 1:9 
        x          = 0.75 + 0.05*k
        C30a[:,k]  .= 10000*max.(x*vec(S30a) .-1,0)
        C30b[:,k]  .= 10000*max.(x*vec(S30b) .-1,0)
        C30c[:,k]  .= 10000*max.(x*vec(S30c) .-1,0)    

        C90a[:,k]  .= 10000*max.(x*vec(S90a) .-1,0)
        C90b[:,k]  .= 10000*max.(x*vec(S90b) .-1,0)
        C90c[:,k]  .= 10000*max.(x*vec(S90c) .-1,0)

        C180a[:,k] .= 10000*max.(x*vec(S180a) .-1,0)
        C180b[:,k] .= 10000*max.(x*vec(S180b) .-1,0)                
        C180c[:,k] .= 10000*max.(x*vec(S180c) .-1,0)
    end 
    return C30a,C30b,C30c, C90a,C90b,C90c ,C180a,C180b,C180c
end


# ---------------------------------------------------------------------------- #
#                                    Output                                    #
# ---------------------------------------------------------------------------- #
#GARCH-M parameters
# const ω	= 0.00001524
# const α	= 0.1883
# const β   = 0.7162
# const λ	= 0.007452

@time B  = BS();
@time C30a,C30b,C30c, C90a,C90b,C90c ,C180a,C180b,C180c = Duan1995(;n_sims=50_000);


T30 = [ mean(C30a,dims=1)[[1;3:7;9]]'
            mean(C30b,dims=1)[[1;3:7;9]]'
            mean(C30c,dims=1)[[1;3:7;9]]']'

T90 = [ mean(C90a,dims=1)[[1;3:7;9]]' 
            mean(C90b,dims=1)[[1;3:7;9]]' 
            mean(C90c,dims=1)[[1;3:7;9]]']'

T180 = [mean(C180a,dims=1)[[1;3:7;9]]' 
            mean(C180b,dims=1)[[1;3:7;9]]'
            mean(C180c,dims=1)[[1;3:7;9]]']'

GH = [ T30; T90; T180]


SX = repeat(B[[1;3:7;9],1],3)

BST = vec(B[[1;3:7;9],2:end])

IND = vec([repeat(["T30"],7) repeat(["T90"],7) repeat(["T180"],7)])

HEADER = ["MAT" "X/S" "BS" "GARCHa" "GARCHb" "GARCHc"]

TABLE41=[HEADER;IND SX BST GH]


# using DataFrames
# TABLE41
