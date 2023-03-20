# ---------------------------------------------------------------------------- #
#                                 RAMPINI 2010                                 #
# ---------------------------------------------------------------------------- #
using LinearAlgebra  #norm function (compute norm of a matrix)
using Plots          #plots
using Roots
using NLsolve
using LaTeXStrings   #to write LaTeX in legends
using Statistics     #kernel density
# using Random         #seed
# using Parameters
using JLD            #to save and load results
# using Distributions

# ---------------------------------------------------------------------------- #
#                                  State Space                                 #
# ---------------------------------------------------------------------------- #

#grid: productivity
Ns     = 5;
S_grid = collect(range(1,Ns,step=1));
π      = 1/Ns;                    # transition prob
Π      = ones(Ns).*π;             # transition prob vector
q      = 1;#ones(Ns);         # q_0 = q_1 = q_2 = 1

#grid: Net worth ω
Nw      = 800
w0_min  = 0.01#1e-3
w0_max  = 2
w0_grid = collect(range(w0_min,w0_max,length=Nw))

χ_1 = 0.08#0.04 #fixed cost

# inv(R)*χ_1

# ---------------------------------------------------------------------------- #
#                                  calibration                                 #
# ---------------------------------------------------------------------------- #
# T               = 3;         #time: t= ∈{0,1,2} (not used yet)
# Nt              = T+1;       #Matrix dimensions 0,...,T+1 (not used yet)
β               = 0.95;      #preferences
θ               = 0.8;       #collateralization rate
α               = 0.33;      #technology
f(k;alpha=α)    =  k  > 0 ? k^alpha : 0;   #technology
∂f∂k(k;alpha=α) = k  > 0 ? alpha*k^(alpha-1) : 999_999;   #technology
R               = 1/β;       #expected return(?)
γ               = q .- inv(R).*(q).*θ ; # minimum down payment: q_t - inv(R)*(E[q_{t+1}])*θ 
# A1(s)           = s./10
ν = 1.
A1(s;ν= ν)    = s.^ν./10 #.- mean(S_grid.^ν./10) .+ 0.3 #volatilidad: mirar ν=0.2 vs ν=1.2
A2              = 1.5
tol             = 1e-9


# A1(S_grid)f

#FUNCIONES AUXILIARES
function λ1(w,s;n_bind_states=Ns)
    k  = w/γ
    if n_bind_states == Ns
        h1        = zeros(Ns);
    elseif n_bind_states == (Ns-1)
        b1        = zeros(Ns)
        b1[1]     = ( (γ + inv(R)*q*θ -4*π*inv(R)*q*θ )*k  - w0  )*inv(π) 
        b1[2:end] .= q*θ*k/R    
        h1        = q*θ*k .- R*b1
    elseif n_bind_states == (Ns-2)
        b1        = zeros(Ns)
        b1[1]     = ( (γ + inv(R)*q*θ -3*π*inv(R)*q*θ )*k  - w0 - π*b2 )*inv(π) 
        b1[2]     = b2
        b1[3:end] .= q*θ*k/R
        h1        = q*θ*k .- R*b1
    elseif n_bind_states == (Ns-3)
        b1        = zeros(Ns)
        b1[1]     = ( (γ + inv(R)*q*θ -2*π*inv(R)*q*θ )*k  - w0 - π*b2 -π*b3)*inv(π) 
        b1[2]     = b2
        b1[3]     = b3
        b1[4:end] .= q*θ*k/R
        h1        = q*θ*k .- R*b1
    elseif n_bind_states == (Ns-4)
        b1        = zeros(Ns)
        b1[1]     = ( (γ + inv(R)*q*θ -π*inv(R)*q*θ )*k  - w0 - π*b2 -π*b3 - π*b4)*inv(π) 
        b1[2]     = b2
        b1[3]     = b3
        b1[4]     = b4
        b1[5]     = q*θ*k/R 
        h1        = q*θ*k .- R*b1
    end

    w1 = A1(s).*f.(k) .+ q.*k.*(1-θ) .+ h1
    k1 = w1./γ  
    R1 = (A1(s).*∂f∂k.(k)  .+ q.*(1-θ))./γ
    R2 = (A2.*∂f∂k.(k1) .+ q.*(1-θ))./γ
    # LHS = inv(R)*Π'*(R1.*R2)
    # RHS = copy(R2)

    LHS = inv(R)*Π'*(R1.*R2) *β^2 
    RHS = copy(R2) *β^2

    return LHS .- RHS
end

#FUNCION λ1 COMPATIBLE CON LIBRERIA NLOPT
function λ_1_aux(F,par;w0=0.25, s=S_grid, n_bind_states=NBS , out=false,k0_constant=2.5676481767282504) 
    #preallocation
    b1        = zeros(Ns)
    if n_bind_states == Ns-1
        #varibles 
        k0_var=par[1]; b1_var1=par[2]
        b1[1]=b1_var1 ; b1[2:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-2
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ;b1[3:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-3
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3] ; b1_var3=par[4]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ; b1[3]=b1_var3 ; b1[4:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-4
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3] ; b1_var3=par[4] ; b1_var4=par[5]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ;  b1[3]=b1_var3 ; b1[4]=b1_var4 ; b1[5] =q*θ*k0_var/R
    elseif n_bind_states == 0
        #varibles 
        k0_var = copy(k0_constant);
        b1_var1=par[1] ; b1_var2=par[2] ; b1_var3=par[3] ; b1_var4=par[4] ; b1_var5=par[5]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ;  b1[3]=b1_var3 ; b1[4]=b1_var4 ; b1[5] =b1_var5

    end

    #COMPUTING EQUATION OF LAMBDA1
    h1 = q*θ*k0_var .- R*b1
    w1 = A1(s).*f(k0_var).+ q.*k0_var.*(1-θ) + h1
    k1 = w1./γ  
    R1 = (A1(s).*∂f∂k.(k0_var) .+ q.*(1-θ))./γ
    R2 = (A2.*∂f∂k.(k1) .+ q.*(1-θ))./γ

    LHS = inv(R)*Π'*(R1.*R2) *β^2 
    RHS = copy(R2) *β^2

    # LHS = inv(R)*Π'*(R1.*R2)*β^2 
    # RHS = R2*β^2

    if out == false
        if n_bind_states == Ns-1
            F[1] = LHS .- RHS[1]                                #FOC
            F[2] = w0 + Π'*b1 - q*k0_var  #Budget constraint
            # F[2] = w0 - (γ*k0_var + inv(R)*Π'*(q*θ*k0_var .- R*b1) ) #Budget constraint
        elseif n_bind_states == Ns-2
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]                                
            F[3] = w0 + Π'*b1 - q*k0_var  #Budget constraint
        elseif n_bind_states == Ns-3 
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]
            F[3] = LHS .- RHS[3]                                                                
            F[4] = w0 + Π'*b1 - q*k0_var  #Budget constraint
        elseif n_bind_states == Ns-4
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]
            F[3] = LHS .- RHS[3]
            F[4] = LHS .- RHS[4]                                                                
            F[5] = w0 + Π'*b1 - q*k0_var  #Budget constraint
        elseif n_bind_states == 0
            # F[1] = w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1))  #equivalent budget constraint                              #FOC        
            F[1] = (LHS .- RHS[1]) - (Π'*b1 - q*k0_var + w0)                           #FOC        
            F[2] = (LHS .- RHS[2]) - (Π'*b1 - q*k0_var + w0) 
            F[3] = (LHS .- RHS[3]) - (Π'*b1 - q*k0_var + w0) 
            F[4] = (LHS .- RHS[4]) - (Π'*b1 - q*k0_var + w0)                                                                
            F[5] = (LHS .- RHS[5]) - (Π'*b1 - q*k0_var + w0) 
            # F[6] = w0 + Π'*b1 - q*k0_var  #Budget constraint
            #OPT2
            # F[1] = (w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1)))^2 /(LHS .- RHS[1])^2 -1                          #FOC        
            # F[2] = (w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1)))^2 /(LHS .- RHS[2])^2 -1
            # F[3] = (w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1)))^2 /(LHS .- RHS[3])^2 -1
            # F[4] = (w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1)))^2 /(LHS .- RHS[4])^2 -1                                                              
            # F[5] = (w0 - (γ*k0_var + Π'*(q*θ*k0_var .- R*b1)))^2 /(LHS .- RHS[5])^2 -1
        end

    elseif out ==true
        return (LHS .- RHS),h1
    end

end

# ---------------------------------------------------------------------------- #
#                                   MAIN LOOP                                  #
# ---------------------------------------------------------------------------- #
function Rampini_capital_collateral(;Nw=Nw,w0_grid=w0_grid,Ns=Ns,S_grid=S_grid,Π=Π,π=π,q=q,β=β,θ=θ,α=α,γ=γ,R=R,A2=A2,tol=tol)


    #preallocation

    #t=0
    B1 = zeros(Ns,Nw)
    H1 = zeros(Ns,Nw)
    Hw1 = zeros(Ns,Nw)
    K0 = zeros(Nw)
    W1 = zeros(Ns,Nw)
    L1 = zeros(Ns,Nw)
    R1 = zeros(Ns,Nw)

    M0 = zeros(Nw)
    M1 = zeros(Ns,Nw)
    M2 = copy(β^2)
    #t=1
    K1 = zeros(Ns,Nw)
    W2 = zeros(Ns,Nw)
    L2 = zeros(Ns,Nw)
    R2 = zeros(Ns,Nw)
    B2 = zeros(Ns,Nw)

    #t=2
    D2 = zeros(Ns,Nw)

    J  =zeros(1)
    NBS=copy(Ns) #este será el contador de lambdas > 0 (partimos con 5)
    println("Number of binding states = $NBS")

    for (j,w) in enumerate(w0_grid)
        w0 = w0_grid[j]
        # println(j)
        # w0=copy(w)

        # ---------------------------------------------------------------------------- #
        #                        CASE I: check λ_1 > 0 ∀ s ∈ S.                        #
        # ---------------------------------------------------------------------------- #
        if sum(λ1(w0,S_grid) .< 0) == 0 # {=1 if λ_1 ≤ 0 ∀ s ∈ S; =0 otherwise}

            #t=0
            H1[:,j] .= 0
            K0[j]   = (w0 .- inv(R)*Π'*H1[:,j])/γ
            B1[:,j] .= q*θ*K0[j]*inv(R)

            w1_aux = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j] .- χ_1
            if sum(w1_aux .> 0)==Ns  #positivo en todos los estados
                W1[:,j] = w1_aux
            else #find k0 and h1
                println("negative net worth")
                Num_neg_w1 = Ns-sum(w1_aux .> 0)              # cuantos estados negativos
                @assert sum(-w1_aux[1:Num_neg_w1]) < sum(w1_aux[Num_neg_w1:Ns])  "Fixed cost too high"  # {=1 if ∃ s' : λ_1(s') = 0, s' ∈ S ; =0 otherwise}

                # h1_aux = zeros(Ns)

                H1[1:Num_neg_w1,j]     = -w1_aux[1:Num_neg_w1] # h1 será tal que w1 no es negativo
                Hw1[Num_neg_w1+1:Ns,j].= -sum(-w1_aux[1:Num_neg_w1])/(Ns-Num_neg_w1) #pedire el mismo monto a los demas estados
                


                W1[:,j] = clamp(A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j] .+ Hw1[:,j] .- χ_1, 0, W1[:,j] )
            end

            L1[:,j] = λ1(w0,S_grid)
            R1[:,j] = (A1(S_grid).*∂f∂k.(K0[j]) .+ q.*(1-θ))./γ

            #t=1
            K1[:,j] = W1[:,j]./γ 
            W2[:,j] = A2.*f.(K1[:,j]).+ q.*K1[:,j].*(1-θ) # .+ H1[:,j]
            R2[:,j] = ((A2.*∂f∂k.(K1[:,j]) .+ q.*(1-θ))./γ)
            B2[:,j] .= q*θ*K1[:,j]*inv(R)

            #t=2
            D2[:,j] = W2[:,j]

            #multipliers
            M1[:,j] = R2[:,j]*M2
            M0[j]   = Π'*(R2[:,j].*M1[:,j])
            L2[:,j] = M1[:,j]/R .- M2
            # L2[:,j] = (R2[:,j]/R .- 1)*β^2

        else
            J = copy(j) #solo como referencia
            K0[j]   = (w0 .- Π'*H1[:,j])/γ
            B1[:,j] .= q*θ*K0[j]*inv(R)
            W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j]
            L1[:,j] = λ1(w0,S_grid)
            NBS = Ns - 1
            println("Number of binding states = $NBS")
            break

        end
    end

    # ---------------------------------------------------------------------------- #
    #                  CASE II: check ∃ s' : λ_1(s') = 0, s' ∈ S.                  #
    # ---------------------------------------------------------------------------- #
    @assert 0 < sum(L1[:,J] .≤ 0) < Ns "Constraint is binding in all states"  # {=1 if ∃ s' : λ_1(s') = 0, s' ∈ S ; =0 otherwise}
        #preallocation
    F=zeros(Ns-NBS+1)
    J2=zeros(1)

    for (j,w) in enumerate(w0_grid[J:end])
        j=J+j-1
        # println(j)
        w0=w0_grid[j]
        #w0=copy(w)

        F=zeros(Ns-NBS+1)
        par_init=[w0/γ  ; repeat([q*θ*w0/γ*inv(R)],Ns-NBS,1)[:]]
        find_k0_b1(F,pars) = λ_1_aux(F,pars;w0=w0,s=S_grid, n_bind_states=NBS)
        res          = nlsolve(find_k0_b1,par_init)

        #t=0
        K0[j]        = res.zero[1]
        B1[1:Ns-NBS,j]  = res.zero[2:(1+Ns-NBS)]
        B1[Ns-NBS+1:Ns,j] = repeat([q*θ*K0[j]*inv(R)],length(Ns-NBS+1:Ns),1)  
        L1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS,out=true)[1]
        H1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS,out=true)[2]
        W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j] .- χ_1
        R1[:,j] = (A1(S_grid).*∂f∂k.(K0[j]) .+ q.*(1-θ))./γ

        #t=1
        K1[:,j] = W1[:,j]./γ 
        W2[:,j] = A2.*f.(K1[:,j]).+ q.*K1[:,j].*(1-θ) #.+ H1[:,j]
        R2[:,j] = ((A2.*∂f∂k.(K1[:,j]) .+ q.*(1-θ))./γ)
        B2[:,j] .= q*θ*K1[:,j]*inv(R)

        #t=2
        D2[:,j] = W2[:,j]

        #multipliers
        M1[:,j] = R2[:,j]*M2
        M0[j]   = Π'*(R2[:,j].*M1[:,j])
        L2[:,j] = M1[:,j]/R .- M2
        # L2[:,j] = (R2[:,j]/R .- 1)*β^2

        if 0 < sum(L1[:,j] .> tol) < Ns #checkar que aun estamos en caso II
            if sum(L1[:,j] .> tol) < NBS #si algun estado deja de estar binding, cambiamos el contador NBS
                NBS = copy(sum(L1[:,j] .> tol))
                println("Number of binding states = $NBS")
            else
                nothing
            end
        else #si pasamos al caso III, todos los lambdas son 0, entonces NBS=0
            J2=copy(j)
            NBS=0
            println("Number of binding states = $NBS")
            break
        end
    end

    # ---------------------------------------------------------------------------- #
    #                        CASE III: check λ_1 = 0 ∀ s ∈ S                       #
    # ---------------------------------------------------------------------------- #
    @assert sum(L1[:,J2] .> tol) == 0 "There is still a state where the constraint is binding"   # {=1 if λ_1 ≤ 0 ∀ s ∈ S; =0 otherwise}

    #find constant k0
    R1_func(k,s) = (A1(s).*∂f∂k(k) .+ q.*(1-θ))./γ
    find_k0_constant(x) = R - Π'*R1_func(x,S_grid)
    k0_bar=find_zero(find_k0_constant,[0,50]) #constant level of capital
    K0[J2:end]        .= k0_bar


    #preallocation
    F=zeros(Ns)


    for (j,w) in enumerate(w0_grid[J2:end])
        j=J2+j-1
        # println(j)
        w0=w0_grid[j]

        #par_init= [k0_bar ; repeat([q*θ*k0_bar*inv(R)],Ns,1)[:] .* [0.8 , 0.7, 0.6, 0.5, 0.4 ]]

        #par_init= [0.5,0.2,0.1,0.2,0.1]
        par_init= repeat([q*θ*k0_bar*inv(R)],Ns,1)[:] .* [0.8 , 0.7, 0.6, 0.5, 0.4 ]
        #par_init= B1[:,j-1]
        find_k0_b1(F,pars) = λ_1_aux(F,pars;w0=w0,s=S_grid, n_bind_states=NBS, k0_constant=k0_bar)
        res          = nlsolve(find_k0_b1,par_init.*.92)
        
        #t=0
        B1[:,j] = copy(res.zero[1:5])
        L1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS, k0_constant=k0_bar,out=true)[1]
        H1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS, k0_constant=k0_bar,out=true)[2]
        W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j] .- χ_1
        R1[:,j] = (A1(S_grid).*∂f∂k.(K0[j]) .+ q.*(1-θ))./γ

        #t=1
        K1[:,j] = W1[:,j]./γ 
        W2[:,j] = A2.*f.(K1[:,j]).+ q.*K1[:,j].*(1-θ) #.+ H1[:,j]
        R2[:,j] = ((A2.*∂f∂k.(K1[:,j]) .+ q.*(1-θ))./γ)
        B2[:,j] .= q*θ*K1[:,j]*inv(R)

        #t=2
        D2[:,j] = W2[:,j]

        #multipliers
        M1[:,j] = R2[:,j]*M2
        M0[j]   = Π'*(R2[:,j].*M1[:,j])
        L2[:,j] = M1[:,j]/R .- M2
        # L2[:,j] = (R2[:,j]/R .- 1)*β^2

    end
    println("solved! :)")
    return B1,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2

end


B1 ,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2 = Rampini_capital_collateral()

# w0_grid[600]
# K0[600]-Π'B1[:,600]
# γ*K0[600] + inv(R).*Π'H1[:,600]

#save for comparision
# cd(dirname(@__FILE__))
# RES=(B1 ,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2)
# save(string(pwd(), "/JLD_files/collateral_capital.jld2"), "capital_results", RES)
# save(string(pwd(), "/JLD_files/collateral_capital_low_vol.jld2"), "capital_results", RES)
# save(string(pwd(), "/JLD_files/collateral_capital_high_vol.jld2"), "capital_results", RES)

# ---------------------------------------------------------------------------- #
#                                     PLOTS                                    #
# ---------------------------------------------------------------------------- #
col  = palette([:skyblue, :blue], 5);

pyplot()
# plotly()
# gr()
# t=0
# t=0
pl_k0=plot(w0_grid,K0, palette=col, legend=:none,title=L"k_{0}",xaxis=(L"w_0"), framestyle = :box)
pl_b1=plot(w0_grid,B1', palette=col,legend=:none,title=L"b_{1}(s)",xaxis=(L"w_0"), framestyle = :box)
pl_l1=plot(w0_grid,L1', palette=col,legend=:none,title=L"\lambda_{1}(s)",xaxis=(L"w_0"), framestyle = :box)
pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s)",xaxis=(L"w_0"), framestyle = :box)
pl_w1=plot(w0_grid,W1', palette=col,legend=:none,title=L"w_{1}(s)",xaxis=(L"w_0"), framestyle = :box)
pl_r1=plot(w0_grid,R1', palette=col,legend=:none,title=L"R_1(s)",xaxis=(L"w_0"), framestyle = :box)

#t=1
pl_k1=plot(w0_grid,K1',palette=col,legend=:none,title=L"k_{1}(s)",xaxis=(L"w_0"), framestyle = :box)
pl_w2=plot(w0_grid,W2',palette=col,legend=:none,title=L"w_2(s)",xaxis=(L"w_0"), framestyle = :box)
pl_r2=plot(w0_grid,R2',palette=col,legend=:none,title=L"R_2(s)",xaxis=(L"w_0"), framestyle = :box)

#t=2
pl_d2=plot(w0_grid,D2',palette=col,legend=:none,title=L"d_2(s)", framestyle = :box)
pl_l2=plot(w0_grid,L2',palette=col,legend=:none,title=L"\lambda_2(s)",xaxis=(L"w_0"), framestyle = :box)

#multipliers
pl_m0=plot(w0_grid,M0,palette=col,legend=:none,title=L"\mu_{0}", framestyle = :box)
pl_m1=plot(w0_grid,M1',palette=col,legend=:none,title=L"\mu_{1}", framestyle = :box)

#cash flow
cash1=θ.*A1(S_grid).*f.(K0')
cash2=θ.*A2.*f.(K1)
pl_cash1=plot(w0_grid,cash1',palette=col,legend=:none,title=L"Cash flow_{t=1}", framestyle = :box)
pl_cash2=plot(w0_grid,cash2',palette=col,legend=:none,title=L"Cash flow_{t=2}", framestyle = :box)

#hedge ratio
pl_hr=plot(w0_grid, (H1./W1)',palette=col,legend=:none,title=L"Hedge Ratio $(h_1/w_1)$",
         framestyle = :box,xaxis=(L"w_0"))

summ=plot(pl_k0,pl_b1,pl_h1,pl_cash1,pl_cash2,pl_k1,pl_w1,pl_w2,pl_hr
    , dpi=1200 , size=(1200,600))





cd(dirname(@__FILE__))
pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s)", framestyle = :box)


cd(dirname(@__FILE__))

summ_k0_h1=plot(pl_k0,pl_h1)
# savefig(summ_k0_h1,"./Figures/summary_k0_h1.pdf")


pyplot()
pl_w1_short=plot(w0_grid[1:50],clamp.(W1[:,1:50]',0,Inf),palette=col,legend=:none,title="Next period net worth, "*L"$w_1(s)$ [zoom]", framestyle = :box)
pl_h1_short=plot(w0_grid[1:50],H1[:,1:50]',palette=col,legend=:none,title="Hedging, "*L"$h_1(s)$ [zoom]", framestyle = :box)
plot(pl_w1_short,pl_h1_short, dpi=1200 , size=(800,600), xaxis="Net worth, "*L"w_0")


pl_w1_all=plot(w0_grid,W1',palette=col,legend=:none,title="Next period net worth, "*L"$w_1(s)$", framestyle = :box)
pl_h1_all=plot(w0_grid,H1',palette=col,legend=:none,title="Hedging, "*L"$h_1(s)$", framestyle = :box)
plot(pl_w1_all,pl_h1_all, dpi=1200 , size=(800,600))


res3=plot(pl_w1_all,pl_h1_all,pl_w1_short,pl_h1_short, dpi=1200 , size=(1200,600))
cd(dirname(@__FILE__))

cd(dirname(@__FILE__))
# savefig(res3,"./Figures/ExFixedWages_K.pdf")





# ---------------------------------------------------------------------------- #
#                                 Distribution                                 #
# ---------------------------------------------------------------------------- #

using Distributions
using Random         #seed
using StatsPlots

ind_w0_h1_A=findfirst(H1[1,:].<tol)
ind_w0_h1_B=findlast(H1[1,:].==0)
ind_w0_h1_C=findfirst(H1[5,:].>tol)


ind_w0_NH=findlast(H1[1,:].==0)
ind_w0_PH=findfirst(H1[5,:].>tol)


#simulate net worth
n_sim= 100_000
w0_sim=clamp.(rand(LogNormal(log(0.5),0.3), n_sim),w0_min,w0_max  ) 


pyplot()
alpha = 0.3

p1=plot(title="Hedging without fixed costs (wages)")
histogram!(w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
 framestyle = :box,label="Net worth distribution")
vspan!([0, w0_grid[ind_w0_NH]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="No hedging")
vspan!([w0_grid[ind_w0_NH],w0_grid[end]], linecolor = :lightgreen,
         fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")
# vspan!([w0_grid[ind_w0_PH],w0_max], linecolor = :green, fillcolor = :green,
#         fillalpha = alpha, label="Fully hedged")


# density(w0_sim)
p2=plot(title= "Hedging with fixed costs (wages)")
histogram!(w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
 framestyle = :box,label="Net worth distribution")

vspan!([w0_grid[1], w0_grid[ind_w0_h1_A]], linecolor = :darkgreen, fillcolor = :darkgreen,
 fillalpha = 0.5, label="hedging states "*L"w_1<0")

vspan!([w0_grid[ind_w0_h1_A], w0_grid[ind_w0_h1_B]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="No hedging")

vspan!([w0_grid[ind_w0_h1_B], w0_grid[end]], linecolor = :lightgreen, fillcolor = :lightgreen,
        fillalpha = alpha, label="hedging")




p=plot(p1,p2,layout=(1,2),size=(1000,600))
cd(dirname(@__FILE__))
savefig(p,"./Figures/dist_wages.pdf")


pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s)", framestyle = :box)
