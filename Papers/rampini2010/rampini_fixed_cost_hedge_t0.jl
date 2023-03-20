# ---------------------------------------------------------------------------- #
#                  EXERCISE: Fixed cost of h1/SIZE                             #
# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
#                                 RAMPINI 2010                                 #
# ---------------------------------------------------------------------------- #
using LinearAlgebra  #norm function (compute norm of a matrix)
using Plots          #plots
using Roots
using NLsolve
using LaTeXStrings   #to write LaTeX in legends
using Statistics     #kernel density
using JLD            #to save and load results
# using Distributions
# using Random         #seed
# using Parameters

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
w0_min  = 1e-3
w0_max  = 2
w0_grid = collect(range(w0_min,w0_max,length=Nw))

# ---------------------------------------------------------------------------- #
#                                  calibration                                 #
# ---------------------------------------------------------------------------- #
# T               = 3;         #time: t= ∈{0,1,2} (not used yet)
# Nt              = T+1;       #Matrix dimensions 0,...,T+1 (not used yet)
β               = 0.95;      #preferences
θ               = 0.8;       #collateralization rate
α               = 0.33;      #technology
R               = 1/β;       #expected return(?)
γ               = q .- inv(R).*(q).*θ ; # minimum down payment: q_t - inv(R)*(E[q_{t+1}])*θ 
A2              = 1.5
tol             = 1e-8
τ               = 0.0 #fixed cost of hedging
# ν               = 1               #variable que incluye volatilidad


f(k;alpha=α)    =  k  > 0 ? k^alpha : 0;   #technology
∂f∂k(k;alpha=α) = alpha*k^(alpha-1);   #technology
A1(s;ν= 1)      = s.^ν./10 #.- mean(S_grid.^ν./10) .+ 0.3 #volatilidad: mirar ν=0.2 vs ν=1.2

#FUNCIONES AUXILIARES
function λ1(w,s;n_bind_states=Ns,θ=θ)
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
function λ_1_aux(F,par;w0=0.25, s=S_grid, n_bind_states=NBS , out=false,k0_constant=2.5676481767282504,θ=θ) 
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
function Rampini_capital_collateral(;Nw=Nw,w0_grid=w0_grid,Ns=Ns,S_grid=S_grid,Π=Π,π=π,q=q,β=β,θ=θ,α=α,γ=γ,R=R,A2=A2,tol=tol,τ=τ,print=0)
    #preallocation


    #t=0
    B1 = zeros(Ns,Nw)
    H1 = zeros(Ns,Nw)
    K0 = zeros(Nw)
    W1 = zeros(Ns,Nw)
    L1 = zeros(Ns,Nw)
    R1 = zeros(Ns,Nw)
    ϵ_k= zeros(Nw)     #overinvestment
    ϵ_b1=zeros(Ns,Nw)  #overborrowing
    ϵ  = zeros(Nw)

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

    print== 1 ? println("Number of binding states = $NBS") : nothing

    for (j,w) in enumerate(w0_grid)
        w0=copy(w)
        # ---------------------------------------------------------------------------- #
        #                        CASE I: check λ_1 > 0 ∀ s ∈ S.                        #
        # ---------------------------------------------------------------------------- #
        if sum(λ1(w0,S_grid,θ=θ) .< 0) == 0 # {=1 if λ_1 ≤ 0 ∀ s ∈ S; =0 otherwise}

            #t=0
            H1[:,j] .= 0
            K0[j]   = (w0 .- Π'*H1[:,j])/γ
            B1[:,j] .= q*θ*K0[j]*inv(R)
            W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j] #Π'W1[:,j]
            L1[:,j] = λ1(w0,S_grid,θ=θ)
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
            L1[:,j] = λ1(w0,S_grid,θ=θ)
            NBS = Ns - 1
            print == 1 ? println("Number of binding states = $NBS") : nothing
            break

        end
    end

    # ---------------------------------------------------------------------------- #
    #                  CASE II: check ∃ s' : λ_1(s') = 0, s' ∈ S.                  #
    # ---------------------------------------------------------------------------- #
    @assert 0 < sum(L1[:,J] .≤ 0) < Ns "Constraint is binding in all states"  # {=1 if ∃ s' : λ_1(s') = 0, s' ∈ S ; =0 otherwise}
        #preallocation
    F  = zeros(Ns-NBS+1)
    J2 = zeros(1)

    for (j,w) in enumerate(w0_grid[J:end])
        #j=1
        j=J+j-1

        # println(j)
        if w0_grid[j]-τ ≥ 0
            # Firma puede pagar el costo τ (no quiebra por pagar el costo)
            # Ahora debe decidir si le conviene o no pagarlo
            # para decidir qué opcion es mas conveniente, miraremos E[w_1|s]
            # la opcion que maximice E[w_1|s], es la opcion que maximiza los dividendos E[d_2|s]

            #option: hedging
            w0_hedge    = w0_grid[j] - τ

            F                  = zeros(Ns-NBS+1)
            par_init           = [w0_hedge/γ  ; repeat([q*θ*w0_hedge/γ*inv(R)],Ns-NBS,1)[:]]
            find_k0_b1(F,pars) = λ_1_aux(F,pars; w0=w0_hedge, s=S_grid, θ=θ, n_bind_states=NBS)
            res                = nlsolve(find_k0_b1,par_init)
            k0_h               = res.zero[1]
            B1[1:Ns-NBS,j]     = res.zero[2:(1+Ns-NBS)]                            #en estos estados h1>0
            B1[Ns-NBS+1:Ns,j]  = repeat([q*θ*k0_h*inv(R)],length(Ns-NBS+1:Ns),1)  #en estos estados h1=0
            L1[:,j]            = λ_1_aux(F,res.zero;w0=w0_hedge,s=S_grid,θ=θ,n_bind_states=NBS,out=true)[1]
            h1_h               = λ_1_aux(F,res.zero;w0=w0_hedge,s=S_grid,θ=θ,n_bind_states=NBS,out=true)[2]
            w1_h               = A1(S_grid).*f(k0_h).+ q.*k0_h.*(1-θ) .+ h1_h
            k1_h = w1_h./γ 
            d2_h = A2.*f.(k1_h).+ q.*k1_h.*(1-θ)

            #option: No hedging
            w0_no_hedge = w0_grid[j]

            h1_nh = 0
            k0_nh  = (w0_no_hedge)/γ
            b1_nh  = q*θ*k0_nh*inv(R)
            w1_nh  = A1(S_grid).*f(k0_nh).+ q.*k0_nh.*(1-θ) .+ h1_nh
            k1_nh = w1_nh./γ 
            d2_nh = A2.*f.(k1_nh).+ q.*k1_nh.*(1-θ)


            ED2H=Π'd2_h
            ED2NH=Π'd2_nh
            # println("w0=$w0_no_hedge ; E[d2|NH]= $ED2NH ; E[d2|H]= $ED2H")
            

            if ED2NH < ED2H #si es mas rentable hacer hedging entonces:
                w0 = copy(w0_hedge)
                H1[:,j] .= h1_h
                K0[j]    = k0_h
                # B1[:,j] .= b1_h
                W1[:,j]  = w1_h        
            else 
                # nothing
                w0 = copy(w0_no_hedge)
                H1[:,j] .= h1_nh
                K0[j]    = k0_nh
                B1[:,j] .= b1_nh
                W1[:,j]  = w1_nh
                L1[:,j] = λ1(w0,S_grid,θ=θ)
            end


        elseif w0_grid[j]-τ < 0 
            # Costo es muy alto entonces firma no paga tau y h=0
            # si h=0, entonces la restriccion es binding
            # println("Firma en quiebra si paga tau")
            w0=w0_grid[j]
            H1[:,j] .= 0
            K0[j]   = (w0 .- Π'*H1[:,j])/γ
            B1[:,j] .= q*θ*K0[j]*inv(R)
            W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j]
            L1[:,j] = λ1(w0,S_grid,θ=θ)
        end

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

        if 0 < sum(L1[:,j] .> tol) < Ns #checkar que aun estamos en caso II
            if sum(L1[:,j] .> tol) < NBS #si algun estado deja de estar binding, cambiamos el contador NBS
                NBS = copy(sum(L1[:,j] .> tol))
                print == 1 ? println("Number of binding states = $NBS") : nothing
            else
                nothing
            end
        else #si pasamos al caso III, todos los lambdas son 0, entonces NBS=0
            J2=copy(j)
            NBS=0
            print == 1 ? println("Number of binding states = $NBS") : nothing
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
    k0_bar=find_zero(find_k0_constant,[0,50]) #+ maximum(ϵ_k)#constant level of capital
    K0[J2:end]        .= k0_bar


    #preallocation
    F=zeros(Ns)


    for (j,w) in enumerate(w0_grid[J2:end])
        j=J2+j-1


        if w0_grid[j]-τ ≥ 0
            # Firma puede pagar el costo τ (no quiebra por pagar el costo)
            # Ahora debe decidir si le conviene o no pagarlo
            # para decidir qué opcion es mas conveniente, miraremos E[w_1|s]
            # la opcion que maximice E[w_1|s], es la opcion que maximiza los dividendos E[d_2|s]

            #option: hedging
            w0_hedge    = w0_grid[j] - τ

            F                  = zeros(Ns-NBS+1)
            par_init= repeat([q*θ*k0_bar*inv(R)],Ns,1)[:] .* [0.8 , 0.7, 0.6, 0.5, 0.4 ]
            find_k0_b1(F,pars) = λ_1_aux(F,pars;w0=w0_hedge,s=S_grid,θ=θ, n_bind_states=NBS, k0_constant=k0_bar)
            res                = nlsolve(find_k0_b1,par_init.*.92)

            B1[:,j] = copy(res.zero[1:5])
            L1[:,j] = λ_1_aux(F,res.zero;w0=w0_hedge,s=S_grid,θ=θ,n_bind_states=NBS, k0_constant=k0_bar,out=true)[1]
            h1_h = λ_1_aux(F,res.zero;w0=w0_hedge,s=S_grid,θ=θ,n_bind_states=NBS, k0_constant=k0_bar,out=true)[2]
            w1_h = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ h1_h
            k1_h = w1_h./γ 
            d2_h = A2.*f.(k1_h).+ q.*k1_h.*(1-θ)

            #option: No hedging
            w0_no_hedge = w0_grid[j]

            h1_nh = 0
            k0_nh  = (w0_no_hedge)/γ
            b1_nh  = q*θ*k0_nh*inv(R)
            w1_nh  = A1(S_grid).*f(k0_nh).+ q.*k0_nh.*(1-θ) .+ h1_nh
            k1_nh = w1_nh./γ 
            d2_nh = A2.*f.(k1_nh).+ q.*k1_nh.*(1-θ)


            ED2H=Π'd2_h
            ED2NH=Π'd2_nh
            # println("w0=$w0_no_hedge ; E[d2|NH]= $ED2NH ; E[d2|H]= $ED2H")
            

            if ED2NH < ED2H #si es mas rentable hacer hedging entonces:
                w0 = copy(w0_hedge)
                H1[:,j] .= h1_h
                # K0[j]    = k0_h
                # B1[:,j] .= b1_h
                W1[:,j]  = w1_h        
            else 
                # nothing
                w0 = copy(w0_no_hedge)
                H1[:,j] .= h1_nh
                # K0[j]    = k0_nh
                B1[:,j] .= b1_nh
                W1[:,j]  = w1_nh
                L1[:,j] = λ1(w0,S_grid,θ=θ)
            end


        elseif w0_grid[j]-τ < 0 
            # Costo es muy alto entonces firma no paga tau y h=0
            # si h=0, entonces la restriccion es binding
            # println("Firma en quiebra si paga tau")
            w0=w0_grid[j]
            H1[:,j] .= 0
            K0[j]   = (w0 .- Π'*H1[:,j])/γ
            B1[:,j] .= q*θ*K0[j]*inv(R)
            W1[:,j] = A1(S_grid).*f(K0[j]).+ q.*K0[j].*(1-θ) .+ H1[:,j]
            L1[:,j] = λ1(w0,S_grid,θ=θ)
        end

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

    end
    print == 1 ? println("solved! :)") : nothing
    return B1, H1, K0, W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2, ϵ_k,ϵ_b1, ϵ

end

τ=0.0025 #0.0025

@time B1 ,H1, K0, W1, L1,
     R1, M0, M1, M2, K1,
      W2, L2, R2, B2, D2,
       ϵ_k,ϵ_b1, ϵ = Rampini_capital_collateral(;τ=τ, print=1)

# w0_grid[600]
# K0[600]-Π'B1[:,600]
# γ*K0[600] + inv(R).*Π'H1[:,600]
# w0_grid[600] - (K0[600]-Π'B1[:,600])

#comparar sum{ d0+βE[d1] + β^2E[d2]}
# pyplot()
# col  = palette([:skyblue, :blue], 5);
# pl_d2=plot(w0_grid,D2',palette=col,legend=:none,title=L"d_2(s) \ (\tau ="*"$τ )", framestyle = :box)

# ---------------------------------------------------------------------------- #
#                              FIGURAS I: dado τ                               #
# ---------------------------------------------------------------------------- #

pyplot()
# col  = palette([:skyblue, :blue], 5);
col  = palette([:pink, :red], 5);

pl_k0=plot(w0_grid,K0,  palette=col, legend=:none,title=L"k_{0} \ (\tau ="*"$τ )", framestyle = :box)
pl_b1=plot(w0_grid,B1', palette=col,legend=:none,title=L"b_{1}(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_l1=plot(w0_grid,L1', palette=col,legend=:none,title=L"\lambda_{1}(s) \ (\tau ="*"$τ )",xaxis=(L"w_0"), framestyle = :box)
pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_w1=plot(w0_grid,W1', palette=col,legend=:none,title=L"w_{1}(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_r1=plot(w0_grid,R1', palette=col,legend=:none,title=L"R_1(s) \ (\tau ="*"$τ )", framestyle = :box)

#t=1
pl_k1=plot(w0_grid,K1',palette=col,legend=:none,title=L"k_{1}(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_w2=plot(w0_grid,W2',palette=col,legend=:none,title=L"w_2(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_r2=plot(w0_grid,R2',palette=col,legend=:none,title=L"R_2(s) \ (\tau ="*"$τ )", framestyle = :box)

#t=2
pl_d2=plot(w0_grid,D2',palette=col,legend=:none,title=L"d_2(s) \ (\tau ="*"$τ )", framestyle = :box)
pl_l2=plot(w0_grid,L2',palette=col,legend=:none,title=L"\lambda_2(s) \ (\tau ="*"$τ )",xaxis=(L"w_0"), framestyle = :box)

#multipliers
pl_m0=plot(w0_grid,M0,palette=col,legend=:none,title=L"\mu_{0} \ (\tau ="*"$τ )", framestyle = :box)
pl_m1=plot(w0_grid,M1',palette=col,legend=:none,title=L"\mu_{1} \ (\tau ="*"$τ )", framestyle = :box)

#cash flow
cash1=A1(S_grid).*f.(K0') 
cash2=A2.*f.(K1)
pl_cash1=plot(w0_grid,cash1',palette=col,legend=:none,title=L"Cash flow_{t=1} \ (\tau ="*"$τ )", framestyle = :box)
pl_cash2=plot(w0_grid,cash2',palette=col,legend=:none,title=L"Cash flow_{t=2} \ (\tau ="*"$τ )", framestyle = :box)

#hedge ratio
pl_hr=plot(w0_grid, (H1./W1)',palette=col,legend=:none,title=L"Hedge Ratio $(h_1/w_1)$"* L" \ ( \tau ="*"$τ )",
         framestyle = :box,xaxis=(L"w_0"))

summ=plot(pl_k0,pl_b1,pl_h1,pl_cash1,pl_cash2,pl_k1,pl_w1,pl_w2,pl_hr,
     dpi=1200 , size=(1200,600))

cd(dirname(@__FILE__))
# savefig(summ,"./Figures/summary_tau_t0.pdf")




# ---------------------------------------------------------------------------- #
#                                EJERCICIO DE FF                               #
# ---------------------------------------------------------------------------- #

Nτ     = 50
τ_min  = 0.0
τ_max  = 0.01
τ_grid = collect(range(τ_min,τ_max,length=Nτ))

#preallocation
ϵKK0 = zeros(Nw,Nτ)
ϵHH1  = zeros(Nw,Nτ)
ϵBB1  = zeros(Ns,Nw,Nτ)

#t=0
BB1 = zeros(Ns,Nw,Nτ)
HH1 = zeros(Ns,Nw,Nτ)
KK0 = zeros(Nw,Nτ)
WW1 = zeros(Ns,Nw,Nτ)
LL1 = zeros(Ns,Nw,Nτ)
RR1 = zeros(Ns,Nw,Nτ)

MM0 = zeros(Nw,Nτ)
MM1 = zeros(Ns,Nw,Nτ)
MM2 = copy(β^2)
#t=1
KK1 = zeros(Ns,Nw,Nτ)
WW2 = zeros(Ns,Nw,Nτ)
LL2 = zeros(Ns,Nw,Nτ)
RR2 = zeros(Ns,Nw,Nτ)
BB2 = zeros(Ns,Nw,Nτ)

#t=2
DD2 = zeros(Ns,Nw,Nτ)


for (j,τ) in enumerate(τ_grid)
    println("tau = $τ")
    # γ = q .- inv(R).*(q).*θ ; # minimum down payment: q_t - inv(R)*(E[q_{t+1}])*θ 
    BB1[:,:,j] ,HH1[:,:,j], KK0[:,j],WW1[:,:,j],
    LL1[:,:,j], RR1[:,:,j], MM0[:,j], MM1[:,:,j],
     MM2, KK1[:,:,j], WW2[:,:,j], LL2[:,:,j],
      RR2[:,:,j], BB2[:,:,j], DD2[:,:,j],
       ϵKK0[:,j],ϵBB1[:,:,j],ϵHH1[:,j] = Rampini_capital_collateral(;τ=τ,print=0)
end

# ---------------------------------------------------------------------------- #
#                                 Distribution                                 #
# ---------------------------------------------------------------------------- #

using Distributions
using Random         #seed
using StatsPlots


#indexes
τ_0_ind  =findmin(τ_grid)[2]
τ_pos_ind=findmin(abs.(τ_grid.-0.002))[2]


#positive τ
τ_pos = round(τ_grid[τ_pos_ind],sigdigits=1)#,digits=3
ind_w0_NH_pos=findfirst(HH1[1,:,τ_pos_ind].>tol) #no hay hedging
ind_w0_PH_pos=findfirst(HH1[5,:,τ_pos_ind].>tol) #hay hedging

#zero τ
τ_0 =τ_grid[τ_0_ind]
ind_w0_NH=findfirst(HH1[1,:,1].>tol) #no hay hedging
ind_w0_PH=findfirst(HH1[5,:,1].>tol) #hay hedging


#tau as % of capital (in the threshold)
tau_mean_k= round(τ_pos / mean(KK0[1:ind_w0_NH_pos,τ_pos_ind]) *100, digits=2)
tau_lim_k = round(τ_pos / KK0[ind_w0_NH_pos,τ_pos_ind] *100, digits=2)
# τ_pos / KK0[1,τ_pos_ind] *100

#tau as % of cashflow (in the threshold)
cash1_pos=Π'*A1(S_grid).*f.(KK0[:,τ_pos_ind]') 

tau_mean_cf=round(τ_pos / mean(cash1_pos) *100, digits=2)
tau_lim_cf =round(τ_pos / cash1_pos[ind_w0_NH_pos] *100, digits=2)
# τ_pos / cash1_pos[1,τ_pos_ind] *100


#simulate net worth
n_sim= 100_000
w0_sim=clamp.(rand(LogNormal(log(0.5),0.3), n_sim),w0_min,w0_max  ) 



#plot distribution
pyplot()
alpha = 0.3

# Distribucion con costo cero
p0=plot(title="Hedging with a fixed transaction cost, "*L" \tau="*"$τ_0")
histogram!(w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
 framestyle = :box,label="Net worth distribution")
vspan!([0, w0_grid[ind_w0_NH]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="Not hedged")
vspan!([w0_grid[ind_w0_NH],w0_max], linecolor = :lightgreen,
         fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")

vline!([w0_grid[ind_w0_NH]],linewidth=2,color=:red, label=L"\tau="*"0% "*L"k_0; "*"0% "*L"CF_1") 

# vspan!([w0_grid[ind_w0_NH],w0_grid[ind_w0_PH]], linecolor = :lightgreen,
#          fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")

# vspan!([w0_grid[ind_w0_PH],w0_max], linecolor = :green, fillcolor = :green,
#         fillalpha = alpha, label="Fully hedged")


p1=plot(title="Hedging with a fixed transaction cost, "*L" \tau="*"$τ_pos") #L" \tau="*"$τ_pos"
histogram!(w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
 framestyle = :box,label="Net worth distribution")
vspan!([0, w0_grid[ind_w0_NH_pos]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="Not hedged")
vspan!([w0_grid[ind_w0_NH_pos],w0_max], linecolor = :lightgreen,
         fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")

# vspan!([w0_grid[ind_w0_NH_pos],w0_grid[ind_w0_PH]], linecolor = :lightgreen,
#          fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")
# vspan!([w0_grid[ind_w0_PH],w0_max], linecolor = :green, fillcolor = :green,
#         fillalpha = alpha, label="Fully hedged")
vline!([w0_grid[ind_w0_NH_pos]],linewidth=2,color=:red, label=L"\tau="*"$tau_lim_k% "*L"k_0; "*"$tau_lim_cf% "*L"CF_1") 


fig_dist_tau=plot(p0,p1,layout=(1,2),size=(1000,600))
cd(dirname(@__FILE__))
# savefig(fig_dist_tau,"./Figures/dist_tau0.pdf")


# savefig(fig_dist_tau,"./Figures/dist_tau0_optK.pdf")

# savefig(fig_dist_tau,"./Figures/dist_tau0_optCF.pdf")





col  = palette([:pink, :red], 5);
pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s) \ (\tau ="*"$τ_0 )", framestyle = :box)



# ---------------------------------------------------------------------------- #
#                            Alternative: %cashflow                            #
# ---------------------------------------------------------------------------- #
#Example: double  xaxis
# gr()
# function twiny(sp::Plots.Subplot)
#     sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
#     plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
#     twinsp = sp.plt.subplots[end]
#     twinsp[:xaxis][:mirror] = true
#     twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
#     Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
#     twinsp
# end
# twiny(plt::Plots.Plot = current()) = twiny(plt[1])
# plot(1:10,rand(10), label = "randData", ylabel = "Y axis",color = :red, legend = :topleft, grid = :off, xlabel = "Numbers Rand")
# p = twiny()
# plot!(p,5:15,log.(5:15), label = "log(x)", legend = :topright, box = :on, grid = :off, xlabel = "Log values")


#indexes
# τ_0_ind=findmin(τ_grid)[2]
# τ_pos_ind=findmin(abs.(τ_grid.-0.002))[2]



ind_sim   = [findmin(abs.(w0_sim[i].-w0_grid))[2] for i in 1:length(w0_sim)]

k0_sim_0    = KK0[ind_sim,τ_0_ind]
cash1_sim_0 = Π'*A1(S_grid).*f.(k0_sim_0) 
k1_sim_0    = KK1[:,ind_sim,τ_0_ind]
cash2_sim_0 = vec(A2.*Π'*f.(k1_sim_0)) 

k0_sim_pos    = KK0[ind_sim,τ_pos_ind]
cash1_sim_pos = Π'*A1(S_grid).*f.(k0_sim_pos) 
k1_sim_pos    = KK1[:,ind_sim,τ_pos_ind]
cash2_sim_pos = vec(A2.*Π'*f.(k1_sim_pos)) 

pyplot()
histogram(k0_sim_0,legend=:none, xaxis="Capital, " *L"k_0",fillalpha = alpha,col=:blue,#background_color_legend = nothing ,
        framestyle = :box,label="Capital investment, "*L" \tau="*"$τ_0")

histogram!(k0_sim_pos,legend=:topleft, xaxis="Capital, " *L"k_0",fillalpha = alpha,col=:red,#background_color_legend = nothing ,
        framestyle = :box,label="Capital investment, "*L"\tau="*"$τ_pos")


histogram(cash1_sim_0,legend=:none, xaxis="Cash flow, " *L"t=1", fillalpha = alpha,col=:blue,#background_color_legend = nothing ,
        framestyle = :box,label="Cash flow, "*L" \tau="*"$τ_0")

histogram!(cash1_sim_pos,legend=:topleft, xaxis="Cash flow, " *L"t=1", fillalpha = alpha,col=:red,#background_color_legend = nothing ,
        framestyle = :box,label="Cash flow, "*L" \tau="*"$τ_pos")


density(cash2_sim_0,legend=:none, xaxis="Net worth, " *L"w_0", fillalpha = alpha,col=:blue,#background_color_legend = nothing ,
        framestyle = :box,label="Net worth distribution")

density!(cash2_sim_pos,legend=:none, xaxis="Net worth, " *L"w_0", fillalpha = alpha,col=:red,#background_color_legend = nothing ,
        framestyle = :box,label="Net worth distribution")


# TAU= 0
plot(w0_grid,KK0[:,τ_0_ind],color=:blue,legend=:none, framestyle = :box)
plot!(w0_grid,KK0[:,τ_pos_ind],color=:red,legend=:none, framestyle = :box)

cash1_0=Π'*A1(S_grid).*f.(KK0[:,τ_0_ind]') 
plot(w0_grid,cash1_0',color=:blue,legend=:none, framestyle = :box)
cash1_pos=Π'*A1(S_grid).*f.(KK0[:,τ_pos_ind]') 
plot!(w0_grid,cash1_pos',color=:red,legend=:none, framestyle = :box)

# TAU = positive
τ_per_k_grid  = range(τ_pos./vec(cash1_pos)[1],τ_pos./vec(cash1_pos)[end],length=length(w0_sim))*100
τ_per_cf_grid = range(τ_pos./vec(KK0[:,τ_pos_ind])[1],τ_pos./vec(KK0[:,τ_pos_ind])[end],length=Nw)*100


collect(τ_per_k_grid)


gr()
alpha = 0.3

# Distribucion con costo cero
plot()#title="Hedging with a fixed transaction cost, "*L" \tau="*"$τ_0"

histogram(τ_per_k_grid,w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
        framestyle = :box,label="Net worth distribution")

histogram(w0_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
        framestyle = :box,label="Net worth distribution")


        p0CF = twiny()
plot!(p0CF,τ_per_k_grid,τ_per_k_grid, label = :none, legend = :none, box = :on, grid = :off, xlabel = "tau as %k")
vspan!([0, w0_grid[ind_w0_NH]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="Not hedged")
vspan!([w0_grid[ind_w0_NH],w0_max], linecolor = :lightgreen,
         fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")

# vspan!([w0_grid[ind_w0_NH],w0_grid[ind_w0_PH]], linecolor = :lightgreen,
#          fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")

# vspan!([w0_grid[ind_w0_PH],w0_max], linecolor = :green, fillcolor = :green,
#         fillalpha = alpha, label="Fully hedged")


p1CF=plot(title="Hedging with a fixed transaction cost, "*L" \tau="*"$τ_pos")
histogram!(cash1_sim,legend=:topright, xaxis="Net worth, " *L"w_0",#background_color_legend = nothing ,
 framestyle = :box,label="Net worth distribution")
vspan!([0, w0_grid[ind_w0_NH_pos]], linecolor = :red1, fillcolor = :red1,
        fillalpha = alpha, label="Not hedged")
vspan!([w0_grid[ind_w0_NH_pos],w0_max], linecolor = :lightgreen,
         fillcolor = :lightgreen,fillalpha = alpha, label="Hedged")















# ---------------------------------------------------------------------------- #
#                                    FIGURES                                   #
# ---------------------------------------------------------------------------- #


# ---------------------------------------------------------------------------- #
#                                Figure: Hedging                               #
# ---------------------------------------------------------------------------- #
pyplot()

OBJ = copy(HH1);yaxis= "Financial slack, "*L"h_1(s)" ; yaxis2= "Financial slack, "*L"E[h_1(s)|s_1]"

OBJ = copy(DD2);  yaxis= "Dividends, "*L"d_2(s)" ; yaxis2= "Dividends, "*L"E[d_2(s)|s_1]"

OBJ = copy(BB1); yaxis= "Debt, "*L"b_1(s)" ; yaxis2= "Debt, "*L"E[b_1(s)|s_1]"

col  = palette([:skyblue, :blue], 5);
pyplot(size=(500,500))
# plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s)",xaxis=L"\omega_0", framestyle = :box)


#firmas grandes: complete hedging (w_0 = 1.0)
w0_ind_big = findmin(abs.(w0_grid.-1))[2]

pl_big=plot(τ_grid,OBJ[:,w0_ind_big,:]',title="Big firms"*L"(w_0=1)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis=yaxis,
         framestyle = :box,legend=:none,palette=col)

#firmas medianas: partial hedging (w_0 = 0.5)
w0_ind_mid = findmin(abs.(w0_grid.-0.5))[2]
pl_mid=plot(τ_grid,OBJ[:,w0_ind_mid,:]',title="Mid firms " *L"(w_0=0.5)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis=yaxis,
         framestyle = :box,legend=:none,palette=col)


#firmas pequeñas: no hedging (w_0 = 0.1)
w0_ind_sml = findmin(abs.(w0_grid.-0.15))[2]
pl_sml=plot(τ_grid,OBJ[:,w0_ind_sml,:]',title="Small firms " *L"(w_0=0.15)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis=yaxis,
         framestyle = :box,label= ["Less productive state" "State 2" "State 3" "State 4" "More productive state"]
         ,legend=:topleft,
         background_color_legend = nothing, palette=col)

fig_obj=plot(pl_sml,pl_mid,pl_big,layout=(1,3),size=(1200,500))


cd(dirname(@__FILE__))
# savefig(fig_obj,"./Figures/tau_d2.pdf")
# savefig(fig_h_theta,"./Figures/heding_tau.png")

pl_big2=plot(τ_grid,(Π'OBJ[:,w0_ind_big,:])',title="Big firms"*L"(w_0=1)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis="Expected "*yaxis2,
         framestyle = :box,legend=:none,color="green")

pl_mid2=plot(τ_grid,(Π'OBJ[:,w0_ind_mid,:])',title="Mid firms " *L"(w_0=0.5)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis="Expected "*yaxis2,
         framestyle = :box,legend=:none,color="green")

pl_sml2=plot(τ_grid,(Π'OBJ[:,w0_ind_sml,:])',title="Small firms " *L"(w_0=0.15)",
        xaxis="Fixed transaction cost, " *L"\tau",yaxis="Expected "*yaxis2,
         framestyle = :box,label= ["Expectation of all states"][]
         ,legend=:topright,
         background_color_legend = nothing,color="green")

fig_obj2=plot(pl_sml2,pl_mid2,pl_big2,layout=(1,3),size=(1200,500))
cd(dirname(@__FILE__))
# savefig(fig_obj2,"./Figures/tau_expec_d2.pdf")


# fig_Eb
# fig_Eh

plot(fig_Eh,fig_Eb, layout=(2,1))










# ---------------------------------------------------------------------------- #
#                                  SURFACE 3D                                  #
# ---------------------------------------------------------------------------- #

OBJ3D = copy(BB1) #BB1; DD2
title = "Dividends"
zaxis = "Dividends"

surf_sml(state,theta) = OBJ3D[state,w0_ind_sml,theta]
x=1:Ns
y=1:Nτ
plotly()
sur_sml=surface( x,y,surf_sml,
         xticks= (collect(1:Ns),string.(S_grid)) ,
         yticks= (Int.(floor.(collect(range(1,Nτ,length=5)))),
         round.(τ_grid[Int.(floor.(collect(range(1,Nτ,length=5))))],digits=2)),
         xaxis="Productivity",
         yaxis="tau",
         zaxis=zaxis,
         title= title*" (a small firm)")

surf_mid(state,theta) = OBJ3D[state,w0_ind_mid,theta]
x=1:Ns
y=1:Nτ
plotly()
sur_mid=surface( x,y,surf_mid,
         xticks= (collect(1:Ns),string.(S_grid)) ,
         yticks= (Int.(floor.(collect(range(1,Nτ,length=5)))),
         round.(τ_grid[Int.(floor.(collect(range(1,Nτ,length=5))))],digits=2)),
         xaxis="Productivity",
         yaxis="tau",
         zaxis=zaxis,
         title= title*" (a mid firm)")

surf_big(state,theta) = OBJ3D[state,w0_ind_big,theta]
x=1:Ns
y=1:Nτ
plotly()
sur_big=surface( x,y,surf_big,
         xticks= (collect(1:Ns),string.(S_grid)) ,
         yticks= (Int.(floor.(collect(range(1,Nτ,length=5)))),
         round.(τ_grid[Int.(floor.(collect(range(1,Nτ,length=5))))],digits=2)),
         xaxis="Productivity",
         yaxis="tau",
         zaxis=zaxis,
         title= title*" (a big firm)")



# ---------------------------------------------------------------------------- #
#                                   Figurita                                   #
# ---------------------------------------------------------------------------- #

#encontrar interseccion
pyplot()

Fixed_cost=0.1

S=5
max_k = 1
min_k = 0.0001
k_grid=min_k:0.001:max_k
ind_inters=findmin(abs.(A1.(S)*f.(k_grid) + q.*k_grid.*(1-θ) .- k_grid))[2]
k2_grid= [k_grid[1:ind_inters-1].*NaN;k_grid[ind_inters:end].-Fixed_cost]
ind_inters2=ind_inters + findmin((A1.(S)*f.(k_grid[ind_inters:end])+ q.*k_grid[ind_inters:end].*(1-θ) .- k2_grid[ind_inters:end]).^2)[2]

fig_fixed_tran_cost=plot(k_grid, A1.(S)*f.(k_grid)+ q.*k_grid.*(1-θ), color=:black,legend=:none, grid=false, box=true,
             xaxis=("Units of "*L"w_0",(0.0 , max_k* 1.5)),title="Components of "*L"w_1(s)")
plot!(k_grid, k_grid, color=:blue,legend=:none)
plot!(k_grid, k2_grid, color=:red,legend=:none)
annotate!(max_k, A1.(S)*f.(max_k)+ q.*max_k.*(1-θ), text(L"A_1(s)f(k_0) + qk_0(1-\theta)", :black, :left, 12))
annotate!(max_k, k_grid[end], text(L"h_1(s)", :black, :left, 12))
annotate!(max_k, k_grid[end]-Fixed_cost, text(L"h_1(s)- \tau", :black, :left, 12))

annotate!(k_grid[ind_inters], k_grid[ind_inters]*1.05, text("A", :black, :right, 12))
annotate!(k_grid[ind_inters2], k2_grid[ind_inters2]*1.05, text("B", :black, :right, 12))
annotate!(k_grid[ind_inters].*.98, mean([k_grid[ind_inters],k2_grid[ind_inters]]), text(L"\tau", :black, :right, 12))

vline!([k_grid[ind_inters]],linestyle=:dash, color=:skyblue)
vline!([k_grid[ind_inters2]],linestyle=:dash, color=:pink)

annotate!(mean([k_grid[ind_inters],k_grid[ind_inters2]]), 0.01*0.1 , text(L"\epsilon^k", :black, :center, 12))


cd(dirname(@__FILE__))
# savefig(fig_fixed_tran_cost,"./Figures/fixed_tran_cost.pdf")

#########################################


