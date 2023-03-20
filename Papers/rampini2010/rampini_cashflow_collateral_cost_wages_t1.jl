# ---------------------------------------------------------------------------- #
#                Cash flow-based collateral with fixed cost χ_t                #
# ---------------------------------------------------------------------------- #

using LinearAlgebra    #norm function (compute norm of a matrix)
using Plots            #plots
using Roots            
using NLsolve          
using LaTeXStrings     #to write LaTeX in legends
using Statistics       #kernel density
using JLD              #to save and load results
# using Distributions  
# using Random         #seed
# using Parameters     

# ---------------------------------------------------------------------------- #
#                                  State Space                                 #
# ---------------------------------------------------------------------------- #


#fixed cost
χ_1 = 1.0;
# χ_2 = 0.0;

#grid: productivity
Ns     = 5;
S_grid = collect(range(1,Ns,step=1));
π      = 1/Ns;                    # transition prob
Π      = ones(Ns).*π;             # transition prob vector
q      = 1;#ones(Ns);             # q_0 = q_1 = q_2 = 1

#grid: Net worth ω
Nw      = 800
w0_min  = 0.73  #calibrado para que el costo fijo se pueda cubrir
w0_max  = 3     #5
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
# γ               = q .- inv(R).*(q).*θ ; # minimum down payment: q_t - inv(R)*(E[q_{t+1}])*θ 
A2              = 1.5
tol             = 1e-5

#FUNCIONES AUXILIARES
# A1(s)           = s./10
ν=1.25
A1(s;ν= ν)   = s.^ν./10 .- mean(S_grid.^ν./10) .+ 0.3 #volatilidad: mirar ν=0.2 vs ν=1.2
f(k;alpha=α)    = k^alpha;   #technology
∂f∂k(k;alpha=α) = alpha*k^(alpha-1);   #technology
γ0(k0)          = q - inv(R)*θ*(Π'*A1(S_grid))*∂f∂k(k0)
γ1(k1)          = q - inv(R)*θ*A2*∂f∂k(k1)


function solve_k0_h1_w1(F,par;W0,s_b,sum_w1)
    Nh  = sum(s_b)
    s_g = (s_b.==0)

    K_0 = par[1]
    H_1 = zeros(Ns)   ; H_1[s_b]    .= par[2:(1+Nh)]             #;H2(s_g) = 0
    W_1 = zeros(Ns)   ; W_1[s_g]    .= par[(2+Nh):(Ns+1)]        #;W2(s_b) = 0

    Chi    = zeros(Ns); Chi[s_b]    .= par[(2+Ns):(1+Nh+Ns)]     #chi(s_g) = 0 (porque w1 es positivo)
    Lambda = zeros(Ns); Lambda[s_g] .= par[(2+Nh+Ns):(2*Ns+1)]   #Lambda(s_b) = 0 (porque estamos en una region de bajo w_0)

    #compute multipliers
    W_1[s_g] = A1(S_grid[s_g]).*f(K_0).*(1-θ) .+ q*K_0 .+ H_1[s_g] .- χ_1
    W_1[s_b] .= 0
    K_1 = zeros(size(W_1)) 
    for i in 1:Ns
        find_k1(k_1) =  W_1[i]  - q*k_1  +  inv(R)*θ*A2*f(k_1)# -  inv(R)*(Π'*H_2)
        # min_k1 = find_k1(1e-15)<0 ?  -1e-15 : 1e-15
        K_1[i] = find_zero(find_k1,[0.1,30])
    end
    R1 = (A1.(S_grid).*∂f∂k(K_0).*(1-θ) .+ q)./ γ0.(K_0)
    R2 = (A2.*∂f∂k.(K_1).*(1-θ) .+ q)./ γ1.(K_1) 

    F[1]                   =  (W0  .- q.*K_0  .+  inv(R)*θ*(Π'* A1.(S_grid)).*f.(K_0) .-  inv(R)*(Π'*H_1)) #.+ (sum_w1.-sum(W_1)) #  ensure eq. budget(s) = 0
    F[2:(1+Nh)]           .=  A1.(S_grid[s_b]).*f.(K_0).*(1-θ) .+ q.*K_0 .- χ_1 .+ H_1[s_b]   .- W_1[s_b] #.+ (sum_w1.-sum(W_1))  #eq. w2(s_b)   = 0
    F[(2+Nh):(Ns+1)]      .=  A1.(S_grid[s_g]).*f.(K_0).*(1-θ) .+ q.*K_0 .- χ_1 .+ H_1[s_g]   .- W_1[s_g] #.+ (sum_w1.-sum(W_1))   #eq .w2(s_g)   = 0
    F[(2+Ns):(1+Nh+Ns)]   .=  inv(R)*Π'*(R1.*R2.*(β^2 .+ Chi)) .- R2[s_b]*β^2  .- Lambda[s_b] .- Chi[s_b]  #.+ (sum_w1.-sum(W_1))#FOC ξ[s_b] > 0, λ[s_g] = 0
    F[(2+Nh+Ns):(2*Ns+1)] .=  inv(R)*Π'*(R1.*R2.*(β^2 .+ Chi)) .- R2[s_g]*β^2  .- Lambda[s_g] .- Chi[s_g]  #.+ (sum_w1.-sum(W_1))#FOC λ[s_g] > 0, ξ[s_b] = 0,

end

function multipliers(w,s;n_bind_states=Ns)
    #compute λ and ξ assuming h_1 = 0
    if n_bind_states == Ns
        h1        = zeros(Ns);
    else
        nothing
    end

    find_k0(k_0) = w - q*k_0 + inv(R)*θ*(Π'*A1(s)).*f(k_0) #budget (since h=0)
    k = find_zero(find_k0,[1e-20,29]) #cash-flow based

    # 2. Cash flow based
    w1 = A1(s).*f.(k).*(1-θ) .+ q.*k .+ h1 .- χ_1

    if sum(w1.<0)==0 # w>0 ∀s ∈ S

        b1 = inv(R) .*θ.*A1(s).*f(k)
        k1 = zeros(size(w1)) 
        for i in 1:Ns
            find_k1(k_1) =  w1[i]  - q*k_1  +  inv(R)*θ*A2*f(k_1) #-  inv(R)*(Π'*h1)
            k1[i] = find_zero(find_k1,[1e-15,30])
        end

        R1 = (A1(s).*∂f∂k(k).*(1-θ) .+ q)./ γ0.(k)
        R2 = (A2*(1-θ)*∂f∂k.(k1) .+ q)./ γ1.(k1) 

        #LAMBDA
        LHS = inv(R)*Π'*(R1.*R2)# *β^2 
        RHS = copy(R2)# *β^2

        return (LHS .- RHS),w1
    elseif sum(w1.<0) != 0 # w(s')<0, s' ∈ S
        lambda=ones(Ns)  # idk this values yet
        lambda[w1.<0].=0 # but these must be 0

        # ξ = zeros(Ns) # idk this values yet
        # ξ[w1.<0] .= 1 # but these must be >0
        return lambda, w1 
    end
end

function λ_1_aux(F,par;w0=0.25, s=S_grid, n_bind_states=NBS , out=false,k0_constant=2.5676481767282504) 
    #preallocation
    b1        = zeros(Ns)
    if n_bind_states == Ns-1
        #varibles 
        k0_var=par[1]; b1_var1=par[2]
        b1[1]=b1_var1 ; b1[2:end] = θ.*A1(S_grid)[2:end].*f(k0_var).*inv(R) #b1[2:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-2
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ; b1[3:end] = θ.*A1(S_grid)[3:end].*f(k0_var).*inv(R)  #b1[3:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-3
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3] ; b1_var3=par[4]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ; b1[3]=b1_var3 ; b1[4:end] = θ.*A1(S_grid)[4:end].*f(k0_var).*inv(R) #b1[4:end].=q*θ*k0_var/R

    elseif n_bind_states == Ns-4
        #varibles 
        k0_var=par[1] ; b1_var1=par[2] ; b1_var2=par[3] ; b1_var3=par[4] ; b1_var4=par[5]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ;  b1[3]=b1_var3 ; b1[4]=b1_var4 ; b1[5] = θ.*A1(S_grid)[5].*f(k0_var).*inv(R) #b1[5] =q*θ*k0_var/R
    elseif n_bind_states == 0
        #varibles 
        k0_var = copy(k0_constant);
        b1_var1=par[1] ; b1_var2=par[2] ; b1_var3=par[3] ; b1_var4=par[4] ; b1_var5=par[5]
        b1[1]=b1_var1 ; b1[2]=b1_var2 ;  b1[3]=b1_var3 ; b1[4]=b1_var4 ; b1[5] =b1_var5

    end

    #COMPUTING EQUATION OF LAMBDA1

    # 1. Capital based
    # h1 = q*θ*k0_var .- R*b1
    # w1 = A1(s).*f(k0_var).+ q.*k0_var.*(1-θ) + h1
    # k1 = w1./γ  
    # R1 = (A1(s).*∂f∂k.(k0_var) .+ q.*(1-θ))./γ
    # R2 = (A2.*∂f∂k.(k1) .+ q.*(1-θ))./γ

    # 2. Cash flow based
    h1 = θ*A1(s).*f(k0_var) .- R*b1
    w1 = A1(s).*f(k0_var).*(1-θ) .+ q.*k0_var + h1 .-χ_1
    k1 = zeros(Ns)  #numerical solution in cash-flow based collateral 
    for i in 1:Ns
        # println(i)
        # find_k1(k_1) = k_1*(γ1.(k_1)) - w1[i]
        # find_k1(k_1) = w1[i] - q*k_1 + (Π'*b1)
        # find_k1(k_1) =  w1[i]  - q*k_1  +  inv(R)*θ*A2*f(k_1)
        find_k1(k_1) =  w1[i]  - q*k_1  +  inv(R)*θ*A2*f(k_1) #-  inv(R)*(Π'*h2)

        k1[i] = find_zero(find_k1,[1e-20,29])
    end
    R1 = (A1(s).*∂f∂k(k0_var).*(1-θ) .+ q)./ γ0.(k0_var)
    R2 = (A2*∂f∂k.(k1)*(1-θ) .+ q)./ γ1.(k1) 

    #lambda
    LHS = inv(R)*Π'*(R1.*R2)# *β^2 
    RHS = copy(R2)# *β^2

    # LHS = inv(R)*Π'*(R1.*R2)*β^2 
    # RHS = R2*β^2

    if out == false
        if n_bind_states == Ns-1
            F[1] = LHS .- RHS[1]                                #FOC
            F[2] = w0 + Π'*b1 - q*k0_var  #Budget constraint
            # F[2] = w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1) #Budget constraint
        elseif n_bind_states == Ns-2
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]   
            F[3] = w0 + Π'*b1 - q*k0_var  #Budget constraint                                         
            # F[3] = w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1)  #Budget constraint
        elseif n_bind_states == Ns-3 
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]
            F[3] = LHS .- RHS[3]
            F[4] = w0 + Π'*b1 - q*k0_var  #Budget constraint                                                                  
            # F[4] = w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1)  #Budget constraint
        elseif n_bind_states == Ns-4
            F[1] = LHS .- RHS[1]                                #FOC        
            F[2] = LHS .- RHS[2]
            F[3] = LHS .- RHS[3]
            F[4] = LHS .- RHS[4]
            F[5] = w0 + Π'*b1 - q*k0_var  #Budget constraint                                                                  
            # F[5] = w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1)  #Budget constraint
        elseif n_bind_states == 0
            # F[1] = w0 + Π'*b1 - q*k0_var  #equivalent budget constraint                              #FOC        
            F[1] = (LHS .- RHS[1]) + (w0 + Π'*b1 - q*k0_var)                           #FOC        
            F[2] = (LHS .- RHS[2]) + (w0 + Π'*b1 - q*k0_var)  
            F[3] = (LHS .- RHS[3]) + (w0 + Π'*b1 - q*k0_var) 
            F[4] = (LHS .- RHS[4]) + (w0 + Π'*b1 - q*k0_var)                                                                 
            F[5] = (LHS .- RHS[5]) + (w0 + Π'*b1 - q*k0_var) 
            
            # F[1] = (LHS .- RHS[1]) + (w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1))                           #FOC        
            # F[2] = (LHS .- RHS[2]) + (w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1))  
            # F[3] = (LHS .- RHS[3]) + (w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1)) 
            # F[4] = (LHS .- RHS[4]) + (w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1))                                                                 
            # F[5] = (LHS .- RHS[5]) + (w0 - q*k0_var + inv(R)*θ*(Π'*A1(s))*f(k0_var) -  inv(R)*(Π'*h1)) 
        end
        # w0 - q*k0_var + inv(R)*θ*(Π'*A1(s)).*f(k0_var) + inv(R)*(Π'*h1) #budget
    # find_k0(k) = k - w0/(γ0(k))
    # find_k0(k_0) =  w0 - q*k_0  +  inv(R)*θ*(Π'*A1(S_grid)).*f(k_0)  

    elseif out ==true
        return (LHS .- RHS),h1
    end

end

function Rampini_fixed_cost_wage()
    # ---------------------------------------------------------------------------- #
    #                                   MAIN LOOP                                  #
    # ---------------------------------------------------------------------------- #
    #preallocation

    #t=0
    B1 = zeros(Ns,Nw)
    H1 = zeros(Ns,Nw)
    K0 = zeros(Nw)
    W1 = zeros(Ns,Nw)
    L1 = zeros(Ns,Nw)
    R1 = zeros(Ns,Nw)
    X1 = zeros(Ns,Nw)

    M0 = zeros(Nw)
    M1 = zeros(Ns,Nw)
    M2 = copy(β^2)
    #t=1
    K1 = zeros(Ns,Nw)
    W2 = zeros(Ns,Nw)
    L2 = zeros(Ns,Nw)
    R2 = zeros(Ns,Nw)
    B2 = zeros(Ns,Nw)
    H2 = zeros(Ns,Nw)

    #t=2
    D2 = zeros(Ns,Nw)

    J  =zeros(1)
    NBS=copy(Ns) #este será el contador de lambdas > 0 (partimos con 5)

    println("Number of binding states = $NBS")
    for (j,w) in enumerate(w0_grid)
        w0=copy(w)
        # w0 = w0_grid[j]  
        λ,w1 = multipliers(w0,S_grid)
        # println(j)
        find_k0(k_0) =  w0 - q*k_0  +  inv(R)*θ*(Π'*A1(S_grid)).*f(k_0) #-  inv(R)*(Π'*H1[:,j])
        k0    = find_zero(find_k0,[0.1,15]) #cash-flow based

        if sum(w1.>-tol) == Ns #positivo en todos los casos?
        # ---------------------------------------------------------------------------- #
        #     CASE I: check λ_1 > 0 ∀ s ∈ S. (if w1<0 it must be true that λ_1 <0)     #
        # ---------------------------------------------------------------------------- #
            if sum(λ .> -tol) == Ns #positivo en todos los casos? entonces h1 =0
                #t=0
                H1[:,j] .= 0

                K0[j]    = find_zero(find_k0,[1e-20,500]) #cash-flow based
                B1[:,j] .= θ*A1(S_grid).*f(K0[j])*inv(R)
                W1[:,j]  = A1(S_grid).*f(K0[j]).*(1-θ)  .+ q.*K0[j].+ H1[:,j] .-χ_1

                R1[:,j]  = (A1(S_grid).*∂f∂k.(K0[j]).*(1-θ) .+ q)./ γ0(K0[j]) 
                #t=1
                for i in 1:Ns
                    find_k1(k_1) =  W1[i,j]  - q*k_1  +  inv(R)*θ*A2*f(k_1) #-  inv(R)*(Π'*H2[:,j])
                    K1[i,j] = find_zero(find_k1,[1e-15,500])
                end
                W2[:,j] = A2.*f.(K1[:,j]).*(1-θ) .+ q.*K1[:,j] #.+ H2[:,j]

                R2[:,j] = ((A2.*∂f∂k.(K1[:,j]).*(1-θ) .+ q)./ γ1.(K1[:,j])) 
                B2[:,j] .= θ*A1.(S_grid).*f.(K1[:,j])*inv(R) 

                #t=2
                D2[:,j] = W2[:,j]

                #multipliers
                L1[:,j]  = multipliers(w0,S_grid)[1]
                X1[:,j]  = zeros(Ns)
                M1[:,j]  = R2[:,j]*M2
                M0[j]    = Π'*(R2[:,j].*M1[:,j])
                L2[:,j]  = M1[:,j]/R .- M2
            elseif sum(λ .> tol) != Ns #negativo en algun caso? entonces pasamos al siguiente loop
                J = copy(j) #solo como referencia
                K0[j]      = find_zero(find_k0,[1e-10,25]) #cash-flow based
                B1[:,j]   .= θ*A1(S_grid).*f(K0[j])*inv(R)
                W1[:,j]    = A1(S_grid).*f(K0[j]).*(1-θ)  .+ q.*K0[j].+ H1[:,j] .-χ_1
                L1[:,j]  = multipliers(w0,S_grid)[1]
                X1[:,j]  = zeros(Ns)
        
                NBS = Ns - 1
                println("Number of binding states = $NBS")
                break
            end
        elseif sum(w1.>tol) != Ns #negativo en algún caso? entonces recalculamos
            s_b= w1.<0 #bad states
            s_g= w1.≥0 #good states

            @assert abs(sum(w1[s_b])) < abs(sum(w1[s_g])) "Fixed cost is too high"

            F    =zeros(1+Ns*2)
            par  = [k0 ; abs.(w1[s_b]) ; w1[s_g] ; ones(Ns)[s_b].*94 ; λ[s_g] ]
            aux_func(F,par) = solve_k0_h1_w1(F,par;W0=w0,s_b=s_b,sum_w1=sum(w1))
            res= nlsolve(aux_func,par;ftol=0.,xtol=0.) #;ftol=0.,xtol=0.           

            #t=0
            K0[j]       = res.zero[1]
            H1[s_b,j]  .= res.zero[2:(1+sum(s_b))]
            W1[s_g,j]  .= res.zero[(sum(s_b)+2):(Ns+1)]
            W1[s_b,j]  .= A1.(S_grid[s_b])*f.(K0[j]).*(1-θ) .+ q.*K0[j] .-χ_1 .+ H1[s_b,j]  #0

            X1[s_b,j] .= res.zero[(2+Ns):(1+sum(s_b)+Ns)]
            L1[s_g,j] .= res.zero[(2+sum(s_b)+Ns):(2*Ns+1)]

            #t=1
            W2[:,j] = A2.*f.(K1[:,j]).*(1-θ) .+ q.*K1[:,j]  #.+ H2[:,j]
            R2[:,j] = ((A2.*∂f∂k.(K1[:,j]).*(1-θ) .+ q)./ γ1.(K1[:,j])) 
            B2[:,j] .= θ*A1.(S_grid).*f.(K1[:,j])*inv(R) 

            #t=2
            D2[:,j] = W2[:,j]

            #multipliers
            L1[:,j]  = multipliers(w0,S_grid)[1]
            X1[:,j]  = multipliers(w0,S_grid)[2]
            M1[:,j]  = R2[:,j]*M2
            M0[j]    = Π'*(R2[:,j].*M1[:,j])
            L2[:,j]  = M1[:,j]/R .- M2
        else
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
        # par_init=[w0/γ0(K0[j-1])/0.1 ; θ*A1(S_grid)[1:(Ns-NBS)]*f(w0/γ0(K0[j-1])) *inv(R)/0.95]
        par_init=[w0/γ0(K0[j-1]) ; B1[1:(Ns-NBS),j-1]]
        # par_init=[w0/γ0(K0[j])    ; repeat([q*θ*w0/γ0(K0[j])  *inv(R)],Ns-NBS,1)[:]]
        find_k0_b1(F,pars) = λ_1_aux(F,pars;w0=w0,s=S_grid, n_bind_states=NBS)
        res          = nlsolve(find_k0_b1,par_init)

        #t=0
        K0[j]             = res.zero[1]
        B1[1:Ns-NBS,j]    = res.zero[2:(1+Ns-NBS)]
        B1[Ns-NBS+1:Ns,j] = θ*A1(S_grid)[Ns-NBS+1:Ns]*f(K0[j]) *inv(R) #this states are still binding
        L1[:,j]           = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS,out=true)[1]
        H1[:,j]           = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS,out=true)[2]
        W1[:,j]           = A1(S_grid).*f(K0[j]).*(1-θ)  .+ q.*K0[j].+ H1[:,j] .-χ_1
        R1[:,j]           = (A1(S_grid).*∂f∂k.(K0[j]).*(1-θ) .+ q)./γ0(K0[j])


        #t=1
        for i in 1:Ns
            # println(i)
            # find_k1(k_1) = W1[i,j] - q*k_1 +inv(R)*θ*A2*f(k_1) 
            find_k1(k_1) =  W1[i,j]  - q*k_1  +  inv(R)*θ*A2*f(k_1) #-  inv(R)*(Π'*H2[:,j])

            K1[i,j] = find_zero(find_k1,[1e-15,20])
        end
        W2[:,j] = A2.*f.(K1[:,j]).*(1-θ) .+ q.*K1[:,j]  # .+ H2[:,j]
        R2[:,j] = (A2.*∂f∂k.(K1[:,j]).*(1-θ) .+ q)./ γ1.(K1[:,j])

        B2[:,j] .= θ*A1.(S_grid).*f.(K1[:,j])*inv(R) 

        #t=2
        D2[:,j] = W2[:,j]

        #multipliers
        M1[:,j] = R2[:,j]*M2
        M0[j]   = Π'*(R2[:,j].*M1[:,j])
        L2[:,j] = M1[:,j]/R .- M2
        # L2[:,j] = (R2[:,j]/R .- 1)*β^2

        if 0 < sum((L1[:,j]) .> tol) < Ns #checkar que aun estamos en caso II
            if sum((L1[:,j]) .> tol) < NBS #si algun estado deja de estar binding, cambiamos el contador NBS
                NBS = copy(sum((L1[:,j]) .> tol))
                println("Number of binding states = $NBS")
            else
                nothing
            end
        else #si pasamos al caso III, todos los lambdas son 0, entonces NBS=0
            J2=copy(j)
            NBS=0
            break
        end
    end


    # ---------------------------------------------------------------------------- #
    #                        CASE III: check λ_1 = 0 ∀ s ∈ S                       #
    # ---------------------------------------------------------------------------- #
    @assert sum(L1[:,J2] .> tol) == 0 "There is still a state where the constraint is binding"   # {=1 if λ_1 ≤ 0 ∀ s ∈ S; =0 otherwise}

    #find constant k0
    R1_func(k,s) = (A1(s).*∂f∂k(k).*(1-θ) .+ q)./γ0(k)
    find_k0_constant(x) = R - Π'*R1_func(x,S_grid)
    k0_bar=find_zero(find_k0_constant,[0.1,500]) #constant level of capital
    K0[J2:end]        .= k0_bar


    #preallocation
    F=zeros(Ns)

    for (j,w) in enumerate(w0_grid[J2:end])
        j=J2+j-1
        # println(j)
        w0=w0_grid[j]

        par_init= B1[:,j-1]
        find_k0_b1(F,pars) = λ_1_aux(F,pars;w0=w0,s=S_grid, n_bind_states=NBS, k0_constant=k0_bar)
        res          = nlsolve(find_k0_b1,par_init)
        
        #t=0
        B1[:,j]           = copy(res.zero[1:5])
        L1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS, k0_constant=k0_bar,out=true)[1]
        H1[:,j] = λ_1_aux(F,res.zero;w0=w0,s=S_grid,n_bind_states=NBS, k0_constant=k0_bar,out=true)[2]
        W1[:,j]           = A1(S_grid).*f(K0[j]).*(1-θ)  .+ q.*K0[j].+ H1[:,j].-χ_1
        R1[:,j]           = (A1(S_grid).*∂f∂k.(K0[j]).*(1-θ) .+ q)./γ0(K0[j])

        #t=1
        for i in 1:Ns
            # println(i)
            # find_k1(k_1) = W1[i,j] - q*k_1 + inv(R)*θ*A2*f(k_1) 
            find_k1(k_1) =  W1[i,j]  - q*k_1  +  inv(R)*θ*A2*f(k_1) #-  inv(R)*(Π'*H2[:,j])

            K1[i,j] = find_zero(find_k1,[1e-15,29])
        end
        W2[:,j] = A2.*f.(K1[:,j]).*(1-θ) .+ q.*K1[:,j] # .+ H2s[:,j]
        R2[:,j] = ((A2.*∂f∂k.(K1[:,j]).*(1-θ) .+ q)./ γ1.(K1[:,j])) 

        B2[:,j] .= θ*A1.(S_grid).*f.(K1[:,j])*inv(R) 


        #t=2
        D2[:,j] = W2[:,j]

        #multipliers
        M1[:,j] = R2[:,j]*M2
        M0[j]   = Π'*(R2[:,j].*M1[:,j])
        L2[:,j] = M1[:,j]/R .- M2
        # L2[:,j] = (R2[:,j]/R .- 1)*β^2
    end
    #     println("solved! :)")
     return B1,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2

end


B1 ,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2 = Rampini_fixed_cost_wage()

#save for comparision
# cd(dirname(@__FILE__))
# RES=(B1 ,H1, K0,W1, L1, R1, M0, M1, M2, K1, W2, L2, R2, B2, D2)
# save(string(pwd(), "/JLD_files/collateral_cashflow.jld2"), "cashflow_results", RES)
# save(string(pwd(), "/JLD_files/collateral_cashflow_low_vol.jld2"), "cashflow_results", RES)
# save(string(pwd(), "/JLD_files/collateral_cashflow_high_vol.jld2"), "cashflow_results", RES)

# ---------------------------------------------------------------------------- #
#                                     PLOTS                                    #
# ---------------------------------------------------------------------------- #
col  = palette([:pink, :red], 5);

# pyplot()
# plotly()
# gr()
#t=0
#t=0
pyplot()
pl_k0=plot(w0_grid,K0, palette=col, legend=:none,title=L"k_{0}", framestyle = :box)
pl_b1=plot(w0_grid,B1', palette=col,legend=:none,title=L"b_{1}(s)", framestyle = :box)
pl_l1=plot(w0_grid,L1', palette=col,legend=:none,title=L"\lambda_{1}(s)", framestyle = :box)
pl_h1=plot(w0_grid,H1', palette=col,legend=:none,title=L"h_{1}(s)", framestyle = :box)
pl_w1=plot(w0_grid,W1', palette=col,legend=:none,title=L"w_{1}(s)", framestyle = :box)
pl_r1=plot(w0_grid,R1', palette=col,legend=:none,title=L"R_1(s)", framestyle = :box)

#t=1
pl_k1=plot(w0_grid,K1',palette=col,legend=:none,title=L"k_{1}(s)", framestyle = :box)
pl_w2=plot(w0_grid,W2',palette=col,legend=:none,title=L"w_2(s)", framestyle = :box)
pl_r2=plot(w0_grid,R2',palette=col,legend=:none,title=L"R_2(s)", framestyle = :box)
pl_h2=plot(w0_grid,H2', palette=col,legend=:none,title=L"h_{2}(s)", framestyle = :box)

#t=2
pl_d2=plot(w0_grid,D2',palette=col,legend=:none,title=L"d_2(s)", framestyle = :box)
pl_l2=plot(w0_grid,L2',palette=col,legend=:none,title=L"\lambda_2(s)", framestyle = :box)

#multipliers
pl_m0=plot(w0_grid,M0,palette=col,legend=:none,title=L"\mu_{0}", framestyle = :box)
pl_m1=plot(w0_grid,M1',palette=col,legend=:none,title=L"\mu_{1}", framestyle = :box)

#cash flow
cash1=θ.*A1(S_grid).*f.(K0')
cash2=θ.*A2.*f.(K1)
pl_cash1=plot(w0_grid,cash1',palette=col,legend=:none,title=L"Cash flow_{t=1}", framestyle = :box)
pl_cash2=plot(w0_grid,cash2',palette=col,legend=:none,title=L"Cash flow_{t=2}", framestyle = :box)

#hedge ratio
pl_hr=plot(w0_grid, (H1./W1)',palette=col,legend=:none,title=L"Hedge \ Ratio \ (h_1/w_1)", framestyle = :box,xaxis=(L"Net \ worth (w_0)"))

summ=plot(pl_k0,pl_b1,pl_h1,pl_cash1,pl_cash2,pl_k1,pl_w1,pl_w2,pl_hr)

# cd(dirname(@__FILE__))
# savefig(summ,"./Figures/summary_CC.pdf")


pyplot()
pl_w1_short=plot(w0_grid[1:50],W1[:,1:50]',palette=col,legend=:none,title="Next period net worth, "*L"$w_1(s)$ [zoom]", framestyle = :box)
pl_h1_short=plot(w0_grid[1:50],H1[:,1:50]',palette=col,legend=:none,title="Hedging, "*L"$h_1(s)$ [zoom]", framestyle = :box)
cd(dirname(@__FILE__))
plot(pl_w1_short,pl_h1_short, dpi=1200 , size=(800,600), xaxis="Net worth, "*L"w_0")


pl_w1_all=plot(w0_grid,W1',palette=col,legend=:none,title="Next period net worth, "*L"$w_1(s)$", framestyle = :box)
pl_h1_all=plot(w0_grid,H1',palette=col,legend=:none,title="Hedging, "*L"$h_1(s)$", framestyle = :box)
cd(dirname(@__FILE__))
plot(pl_w1_all,pl_h1_all, dpi=1200 , size=(800,600))


res3=plot(pl_w1_all,pl_h1_all,pl_w1_short,pl_h1_short, dpi=1200 , size=(1200,600))
cd(dirname(@__FILE__))
# savefig(res3,"./Figures/ExFixedWages.pdf")
