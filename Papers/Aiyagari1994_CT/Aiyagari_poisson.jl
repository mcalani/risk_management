using Plots
using Statistics
using LinearAlgebra
using SparseArrays
using Parameters
# Value function:
# V(t,x) = max_{α} ∫ ℯ^{-ρ(s-t)} * u(α,x) ds

# HJB:
# ρV(x) = ∂V/∂t + max_{α}{ u(α,x) + ∑_{n=1,N}(μ_{t,n}*∂V/∂x_n) }

# ---------------------------------------------------------------------------- #
#                                  Parameters                                  #
# ---------------------------------------------------------------------------- #

function parameters(;   α = 1/3,
                        β = 1.05^(-1/4),
                        δ = 0.025,
                        γ = 2,
                        μ = 0.5,             
                        ϕ = 0,  
                        ρ = 1-β,
                        Γ = Inf,

                        # Grid of assets
                        Na    = 200,
                        a_max = 500,        
                        a_min = copy(ϕ),
                        a_grid = collect(range(a_min,a_max,length=Na)),
                        λ_e = 0.05, # An employed individual becomes unemployed with probability
                        λ_u = 0.8,  # An unemployed individual becomes employed with probability
                        )

        # ---------------------------------------------------------------------------- #
        #                               Transition Matrix                              #
        # ---------------------------------------------------------------------------- #

        T1 = [1-λ_e    λ_u;
                λ_e   1-λ_u ]
        # Convert into continuous time dg_t=T*g_t

        T = [-λ_e λ_u;
                λ_e -λ_u]

        Λ=[-λ_e*I(Na)  λ_e*I(Na); 
        λ_u*I(Na)  -λ_u*I(Na)] #will be used later
        # The steady state values of g satisfy 0=T*g. That is g is an eigenvector
        # associated with a zero eigenvalue normalized to sum to one.

        eigen_0=eigen(T).vectors[:,2]

        g = eigen_0/sum(eigen_0)

        e_ss = g[1]
        u_ss = g[2]

        # Tax rate to balance budget: (1-e)*mu*w=e*tau*w.
        τ = (1-e_ss)/e_ss * μ 


        # Utility and consumption functions
        u(c)  = (c.^(1-γ) .-1)./(1-γ)
        ∂u∂c(c) = c.^(-γ)
        c(dv) = dv.^(-1/γ) # from ∂u(c)/∂c = ∂V(a)/∂a

        # ---------------------------------------------------------------------------- #
        #                                  Derivatives                                 #
        # ---------------------------------------------------------------------------- #

        # Forward  : dv/da = ( v(a+Δa) - v(a) )/ (Δa)
        # Backward : dv/da = ( v(a) - v(a-Δa) )/ (Δa)
        # Central  : dv/da = ( v(a+Δa) - v(a-Δa) )/ (2Δa)

        # Forward and Backward derivatives of V such that D*V ≈ ∂V/∂a
        Δa = a_grid[2]-a_grid[1];

        D_f = [(i==Na && j ==Na) ? -1 : i==j ? -1/Δa : i==j-1 ? 1/Δa : 0 for i = 1:Na, j=1:Na]
        D_b = [(i==1 && j == 1) ? 1 : i==j ? 1/Δa : i==j+1 ? -1/Δa : 0 for i = 1:Na, j=1:Na]

        # ---------------------------------------------------------------------------- #
        #                                Solve the model                               #
        # ---------------------------------------------------------------------------- #
        #Guess on value function (we solve for employed and unemployed separately)

        #Guess on prices
        r_h = 1/β -1
        r_l = 0

        r = mean([r_h r_l])
        w = (1-α)*(((r+δ)/α )^(1/(α-1)))^(α);

        V_e = u(w*(1-τ) .+ r.*a_grid) ./ ρ
        V_u = u(w*μ .+ r.*a_grid) ./ ρ

        #preallocation
        P = spzeros(2*Na,2*Na)
        a_dot_e = zeros(Na)
        a_dot_u = zeros(Na)
        c_e = zeros(Na)
        c_u = zeros(Na)
        G = zeros(2*Na)
        K = 0.

        tol_r = 1e-8
        tol_v = 1e-8
        return (α =α, β=β ,δ=δ ,γ=γ ,μ=μ ,ϕ=ϕ ,ρ=ρ ,Γ=Γ,λ_e=λ_e ,λ_u=λ_u ,Na=Na ,a_max=a_max,a_min=a_min,a_grid=a_grid,
                r_h=r_h,r_l=r_l,Δa=Δa,D_f=D_f,D_b=D_b,r=r,w=w,V_e=V_e, V_u=V_u,Λ=Λ, T=T,e_ss=e_ss,u_ss=u_ss,τ=τ,
                u=u, c=c, ∂u∂c=∂u∂c,
                P=P, a_dot_e=a_dot_e,a_dot_u=a_dot_u, c_e=c_e, c_u=c_u, G=G, K=K,
                tol_v=tol_v,tol_r=tol_r)
end


function Aiyagari_poisson(ec)
        @unpack α,β,δ,γ,μ,ϕ,ρ,Γ,λ_e,λ_u,Na,a_max,a_min,a_grid = ec
        @unpack r_h,r_l,Δa,D_f,D_b,r,w,V_e, V_u,Λ, T,e_ss,u_ss,τ = ec
        @unpack u, c, ∂u∂c = ec
        @unpack tol_v,tol_r,P = ec

        conv_r=100
        for i in 1:1000
                r = mean([r_h r_l])
                w = (1-α)*(((r+δ)/α )^(1/(α-1)))^(α)

                #begin while loop
                conv=100;
                # V_e = zeros(Na)
                # V_u = zeros(Na)
                for j in 1:1000

                        # ---------------------------------------------------------------------------- #
                        #                                   Employed                                   #
                        # ---------------------------------------------------------------------------- #
                        ∂V_ef = D_f*V_e ; ∂V_ef[end] = ∂u∂c(w*(1-τ) .+ r.*a_grid[end])#FOC 
                        ∂V_eb = D_b*V_e ; ∂V_eb[1]   = ∂u∂c(w*(1-τ) .+ r.*a_grid[1])  #FOC

                        #derivatives
                        c_ef     = c(∂V_ef)
                        c_eb     = c(∂V_eb)

                        a_dot_ef = w*(1-τ) .+ r.*a_grid .- c_ef
                        a_dot_eb = w*(1-τ) .+ r.*a_grid .- c_eb

                        #Indeces: direction of the state ̇a
                        I_ef = vec(a_dot_ef .>0)
                        I_eb = vec(a_dot_eb .<0)

                        #Optimal decisions
                        c_e = I_ef.*c_ef + I_eb.*c_eb + (1 .-I_ef-I_eb).*(w*(1-τ) .+ r.*a_grid) 
                        copyto!(ec.c_e,I_ef.*c_ef + I_eb.*c_eb)
                        #note: the last term only applies when I_ef = I_eb = 0 (steady state)

                        a_dot_e = spdiagm(0 =>I_ef.*a_dot_ef) + spdiagm(0 => I_eb.*a_dot_eb)
                        copyto!(ec.a_dot_e, I_ef.*a_dot_ef + I_eb.*a_dot_eb)
                        
                        # ---------------------------------------------------------------------------- #
                        #                                  Unemployed                                  #
                        # ---------------------------------------------------------------------------- #
                        ∂V_uf = D_f*V_u ; ∂V_uf[end] = ∂u∂c(w*μ .+ r.*a_grid[end])#FOC
                        ∂V_ub = D_b*V_u ; ∂V_ub[1]   = ∂u∂c(w*μ .+ r.*a_grid[1]) 

                        #derivatives
                        c_uf     = c(∂V_uf)
                        c_ub     = c(∂V_ub)

                        a_dot_uf = w*μ .+ r.*a_grid .- c_uf
                        a_dot_ub = w*μ .+ r.*a_grid .- c_ub

                        #Indeces: direction of the state ̇a
                        I_uf = vec(a_dot_uf .>0)
                        I_ub = vec(a_dot_ub .<0)
                        #optimal decisions
                        c_u = I_uf.*c_uf + I_ub.*c_ub + (1 .-I_uf-I_ub).*(w*μ .+ r.*a_grid)
                        copyto!(ec.c_u,I_uf.*c_uf + I_ub.*c_ub)

                        a_dot_u = spdiagm(0 =>I_uf.*a_dot_uf) + spdiagm(0 => I_ub.*a_dot_ub)
                        copyto!(ec.a_dot_u, I_uf.*a_dot_uf + I_ub.*a_dot_ub)

                        # ---------------------------------------------------------------------------- #
                        #                                  Aggregate                                   #
                        # ---------------------------------------------------------------------------- #

                        A_dot    =[ a_dot_e          spzeros(Na,Na) ;
                                    spzeros(Na,Na)    a_dot_u       ]

                        D    =[spdiagm(0 =>Float64.(I_ef))*sparse(Float64.(D_f)) + spdiagm(0 =>Float64.(I_eb))*sparse(Float64.(D_b))  spzeros(Na,Na)    ;
                                spzeros(Na,Na)     spdiagm(0 =>Float64.(I_uf))*sparse(Float64.(D_f)) .+ spdiagm(0 =>Float64.(I_ub))*sparse(Float64.(D_b)) ]

                        #We can write 
                        # ρV_{n+1} = u(c_n) + A_dot*D*V_{n+1} - Λ*V_{n+1}
                        # ρV_{n+1} = u(c_n) + P_n V_{n+1}, where P_n = A_dot*D + Λ 
                        #so we iterate on V_{n+1} = [(ρ+1/Γ)I - P_n]^-1 [u(c_n) - v_n/Γ] = X\y 
                        P = A_dot*D+ Λ

                        A = (1/Γ+ρ)*I(Na*2)-P
                        b = [u.(c_e) ; u.(c_u)] .+[V_e;V_u] ./Γ

                        V = Array(A)\b

                        V_en = V[1:Na];
                        V_un = V[Na+1:end];

                        conv = maximum(abs.([V_en-V_e V_un-V_u])); #println(conv)
                        V_e = copy(V_en)
                        V_u = copy(V_un)

                        if conv < tol_v
                                break
                        end
                        #display([ V_e V_u])
                        #display([ V_en V_un])
                end

                # ---------------------------------------------------------------------------- #
                #                             Kolmogorov Forward Eq                            #
                # ---------------------------------------------------------------------------- #
                # States evolves according to ̇g = P'g, so we must find 0=P'g
                # and normalize to one
                G = eigen(Array(P')).vectors[:,end]/ sum(eigen(Array(P')).vectors[:,end])
                copyto!(ec.G,G)
                
                #Total amount of capital
                K = real.(G[1:Na]')*a_grid + real.(G[Na+1:end]') * a_grid
                #copyto!(ec.K,K)
                #New price ̂r
                r_prime = α.*(K ./e_ss).^(α-1) .- δ
                #copyto!(ec.r,r_prime)

                frs = r_prime[1] .- r

                if frs > 0
                        r_l = copy(r);
                else
                        r_h = copy(r);
                end


                conv_r = abs(frs) ; println(conv_r)
                if conv_r < tol_r
                        break
                end
        end
end


ec = parameters(Na=500)
Aiyagari_poisson(ec)


p11=plot(ec.a_grid[1:20],[ ec.a_dot_e[1:20] ec.a_dot_u[1:20]],legend=:none,title="Savings"  );
p21=plot(ec.a_grid[1:20],[ ec.c_e[1:20] ec.c_u[1:20]],legend=:none,title="Consumption");

p1=plot(p11,p21)

factor=5
grid=ec.a_grid[1:Int(floor(ec.Na/factor))]
dist_e = ec.G[1:Int(floor(ec.Na/factor))]
dist_u = ec.G[ec.Na+1:(ec.Na +Int(floor(ec.Na/factor)))] 
p2=plot(grid,[ dist_e dist_u],legend=:none,title="Distribution" );

p3=plot(p1,p2)

savefig(p3,"aiyagari.poisson.svg")