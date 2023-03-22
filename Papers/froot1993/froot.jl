using LinearAlgebra   #norm function (compute norm of a matrix)
using Plots           #plots
using Roots
using NLsolve
using LaTeXStrings    #to write LaTeX in legends
using Statistics      #kernel density
using JLD             #to save and load results
using Distributions
# using Random        #seed
# using Parameters


#NOTA: Sensibilidades del modelo
# i)   Parámetro de correlacion α: entre 0.2-0.3 da resultados con más sentido
# ii)  Desviación estándar de los flujos σ: tiene que ser pequeño pq si es muy grande (>0.2 aprox), h* tiende a ser negativo
# iii) Grilla de w0: se tiene que cumplir que I ≥ w0 ( I = w + e ).
#      En ocasiones, la inversión optima es negativa fuera de la grilla w0_grid= 0.0:0.75
#      Esta grilla se puede calibrar (medianamente) moviendo α y σ.
# iv) También se pueden mover los parámetros de la función de profit,γ y ν, para modificar la concavidad de P(w)

# ---------------------------------------------------------------------------- #
#        Stochastic processes: w = w0(h+(1-h)ϵ) , θ = α(ϵ-̄ϵ)+1, ϵ∼N(1,σ^2)     #
# ---------------------------------------------------------------------------- #

#ϵ process
const σ       = 1/5
const ϵ_bar   = 1
const N_ϵ     = 5
const ϵ_grid  = collect(range(1-2*σ, 1+2*σ, length=5))
const π       = 1/N_ϵ

#θ process
const α       = 0.3 #ENTRE 2 Y 3 PARA QUE DE ALGO RAZONABLE.
const θ_grid  = α*(ϵ_grid .- ϵ_bar) .+ 1

# ---------------------------------------------------------------------------- #
#                                  Parameters                                  #
# ---------------------------------------------------------------------------- #
#wealth grid: la definimos al final
# w0_min   = 0.05
# w0_max   = 2.00
# w0_grid  = collect(w0_min:0.01:w0_max)
# Nw       = length(w0_grid)

#net present value of investment expenditures F(I)
const γ           = 1/3;      #technology
f(I;gamma=γ)= I  > 0 ?  (1/gamma)*I^gamma : 0;   #technology
F(I;θ=θ_grid)        = θ_grid.*f.(I) .- I

#cost of external financing C(e)
const ν           = (1+γ)
C(e;ν=ν)    = e  > 0 ? (e).^ν : 0

# x_grid  = collect(0.01:0.001:10)
# plot(x_grid,C.(x_grid),legend=:bottomright, xaxis=L"x", yaxis="C(x)",label="C(e)")
# plot!(x_grid,F.(x_grid), xaxis=L"x", yaxis="P(w)",label="F(I)")


#profit function
function P(w ; e=0, θ=θ_grid)
    I      = w + e
    Profit = F.(I;θ=θ) .- C.(e)
    return Profit
end

# pyplot()# plotly()
# plot(w0_grid ,P.(w0_grid),legend=:bottomleft, xaxis=L"w", yaxis=L"P(w)",title="Profit function",label="e=0")
# plot!(w0_grid,P.(w0_grid;e=0.01), xaxis=L"w" , yaxis=L"P(w)",title="Profit function",label="e=0.01")
# plot!(w0_grid,P.(w0_grid;e=0.1), xaxis=L"w"  , yaxis=L"P(w)",title="Profit function",label="e=0.1")


# ---------------------------------------------------------------------------- #
#                                  Derivatives                                 #
# ---------------------------------------------------------------------------- #
# first derivatives
f_I(I)  =  I  > 0 ? I.^(γ-1) : 0
C_e(e)  =  e  > 0 ? ν*e.^(ν-1) : 0

# second derivatives
f_II(I) =  I  > 0 ? (γ-1)*I.^(γ-2) : -0
C_ee(e) =  e  > 0 ? ν*(ν-1)*e.^(ν-2) : -0

# derivative dI*/dw
Iast_w(Iast,w,ϵ) = -C_ee(Iast-w) ./ ( (α*(ϵ .- ϵ_bar) .+1).*f_II(Iast) .- C_ee(Iast-w) )

#OPT1: dan igual
# P_ww(Iast,w,ϵ) = (α*(ϵ .- ϵ_bar) .+1).*f_II(Iast).*Iast_w(Iast,w,ϵ).^2 .- C_ee(Iast-w).*( Iast_w(Iast,w,ϵ).- 1).^2

#OPT2: dan igual
P_ww(I,w,ϵ) = (α*(ϵ .- ϵ_bar) .+1).*f_II.(I).*Iast_w.(I,w,ϵ)




#funcion auxiliar: condicion de primer orden en t=1
function FOC_t_1(F , par ; w0=0.01, hi=0, ϵ_grid=ϵ_grid)
    #determine number of equations
    N     = length(ϵ_grid)

    #determine parameters
    I_ast = par[1:N]
    h_ast = hi #par[N]

        #F.O.C t=1: decide optimal I
    w1      = w0.*(h_ast .+ (1-h_ast).*ϵ_grid)
    F[1:N]  = ( (α.*(ϵ_grid .- ϵ_bar) .+1).*f_I.(I_ast) .- 1 .-  C_e.(I_ast.-w1) ).^2
end


# ---------------------------------------------------------------------------- #
#                        Resolviendo el modelo FSS 1993                        #
# ---------------------------------------------------------------------------- #

function FOC_t_0(par;w0=0.75, optim=1)
    hi = par[1]
    # hi = copy(h_ast)
    # w0 = 0.75
    
    # ---------------------------------------------------------------------------- #
    #                         I) Find optimal investment I*                        #
    # ---------------------------------------------------------------------------- #
    par_init         = [0.8 , 0.8, 0.8, 0.8, 0.8]
    G                = zeros(N_ϵ)
    find_I(G,pars)   = FOC_t_1(G , pars ; w0=w0, hi=hi, ϵ_grid=ϵ_grid)
    res              = nlsolve(find_I,par_init;show_trace=false,xtol=1e-10,ftol=1e-10)
    I_ast            = res.zero

    # ---------------------------------------------------------------------------- #
    #                          II) Find optimal hedging h*                         #
    # ---------------------------------------------------------------------------- #
    w1    = w0.*(hi .+ (1-hi).*ϵ_grid)

    #CHECK sum(I_ast .> w1) ==  N_ϵ

    num   = mean( f_I.(I_ast) .* P_ww.(I_ast,w1,ϵ_grid)./ ((α*(ϵ_grid .- ϵ_bar) .+1).*f_II.(I_ast))  )
    den   =  w0*mean(P_ww.(I_ast,w1,ϵ_grid))       
    h_ast = 1 + α * num / den

    if optim == 1
        return h_ast - hi
    else
        I_check=sum(I_ast .> w1)
        return h_ast, I_check, I_ast
    end
end


# ******** DETERMINAR LOS PUNTOS DE GRILLA DONDE SE PUEDE SOLUCIONAR RESOLVIENDO PARA PUNTOS EXTREMOS ********
# w0_aux= 0.05
# f_aux(par) = find_h(par;w0=copy(w0_aux), optim=1)
# h_ast      = find_zero(f_aux,[-50,50],verbose=true,xtol=1e-12,ftol=1e-12) #constant level of capital
# h_check, I_check, Iast=find_h(h_ast;w0=w0_aux, optim=0)
# h_ast-h_check< tol
# ************************************************************************************************************

function FSS_1993(w0_grid;h="optim")
    Nw       = length(w0_grid)
    #preallocation
    H_opt = zeros(Nw)
    I_opt = zeros(Nw, N_ϵ)
    w1_grid=zeros(size(I_opt))

    if h=="optim"
        @fastmath @inbounds begin
            for (i,w0) in enumerate(w0_grid)
                # w0 = w0_grid[i]

                f_aux(par)              = FOC_t_0(par;w0=w0, optim=1)
                h_ast                   = find_zero(f_aux,[-50,50],verbose=false,xtol=1e-12,ftol=1e-12) #constant level of capital
                h_check, I_check, I_ast = FOC_t_0(h_ast;w0=w0, optim=0)

                H_opt[i]   = h_ast
                I_opt[i,:] = I_ast

                err = h_ast-h_check
                println("Solved for w0=$w0 (round N° $i/$Nw) ; error = $err")
            end
        end

    else 
        H_opt .= copy(h) 
        for (j,ϵ) in enumerate(ϵ_grid) 
            w1_grid[:,j] = w0_grid.*(h .+ (1 .- h)*ϵ )
        end

        #optimal investment given h = x
        par_init         = [0.8 , 0.8, 0.8, 0.8, 0.8]
        G                = zeros(N_ϵ)
        for (i,w0) in enumerate(w0_grid)
            find_I(G,pars)   = FOC_t_1(G , pars ; w0=w0, hi=h, ϵ_grid=ϵ_grid)
            res              = nlsolve(find_I,par_init;show_trace=false,xtol=1e-10,ftol=1e-10)
            I_res            = res.zero
            I_opt[i,:]        = I_res
            println("Solved for w0=$w0 (round N° $i/$Nw))")
        end
    end

    # recovering optimal values
    F_opt = zeros(size(I_opt))
    for (j,ϵ) in enumerate(ϵ_grid) 
        w1_grid[:,j] = w0_grid.*(H_opt .+ (1 .- H_opt)*ϵ )
        θ            = α*(ϵ .- ϵ_bar) .+ 1
        F_opt[:,j]   = θ*f.(I_opt[:,j]) .- I_opt[:,j]
    end

    e_opt = I_opt .- w1_grid
    C_opt = C.(e_opt)
    P_opt = F_opt .- C_opt
    return H_opt, I_opt, e_opt, F_opt, C_opt, P_opt, w1_grid
end


#wealth grid
const Nw       = 100
const w0_min   = 0.05
const w0_max   = 0.8#.95
const w0_grid  = collect(range(w0_min,w0_max,length=Nw))


@time H_opt, I_opt, e_opt, F_opt, C_opt, P_opt, w1_grid = FSS_1993(w0_grid;h="optim")


##### CHECK I ≥ w #####
@assert sum(I_opt .< w1_grid) == 0 " I ≥ w must always be true !!!! "



# ---------------------------------------------------------------------------- #
#                                    FIGURE                                    #
# ---------------------------------------------------------------------------- #
pyplot()
col  = palette([:skyblue, :blue], 5);

# plotly()
ind=findmin(H_opt.^2)[2] #para mirar valores positivos de h
xaxis="cash flow, "*L"w_0";
yaxis_h = [0 , 1]
yaxis_w = [0 , 1]
yaxis_e = [0 , 0.3]
yaxis_I = [0.25 , 1.2]
yaxis_F = [.75, 2.5]
yaxis_P = [1.2 , 2.5]




pl_h=plot(w0_grid[ind:end],H_opt[ind:end], color=:black   ,legend=:none,yaxis=yaxis_h, title= "Hedging ratio, "*L"h^*"     ,xaxis=xaxis);
pl_w=plot(w0_grid[ind:end],w1_grid[ind:end,:], palette=col,legend=:none,yaxis=yaxis_w, title= "Cash-flow realizations, "*L"w"     ,xaxis=xaxis);
pl_e=plot(w0_grid[ind:end],e_opt[ind:end,:], palette=col,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(e_opt[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
pl_I=plot(w0_grid[ind:end],I_opt[ind:end,:], palette=col,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(I_opt[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
pl_F=plot(w0_grid[ind:end],F_opt[ind:end,:], palette=col,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(F_opt[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
pl_P=plot(w0_grid[ind:end],P_opt[ind:end,:], palette=col,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w_0)"          ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(P_opt[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w0)"          ,xaxis=xaxis);

pl_C=plot(w0_grid[ind:end],C_opt[ind:end,:], palette=col,legend=:none,yaxis=[0 , .15], title= "Costs, "*L"C(e^*)"          ,xaxis=xaxis)
# plot!(w0_grid[ind:end],mean(C_opt[ind:end,:],dims=2),  color=:red,legend=:none,yaxis=[0 , .2], title= "Costs, "*L"C(e^*)"          ,xaxis=xaxis)

summ=plot(pl_h,pl_w,pl_I,pl_e,pl_F,pl_P,layout=(2,3),dpi=600 ,size=(900,600))


# cd(dirname(@__FILE__))
# savefig(summ,"./Figures/sum_froot.pdf")



# ---------------------------------------------------------------------------- #
#                             contrafactual: h = 0                             #
# ---------------------------------------------------------------------------- #

@time H_nh, I_nh, e_nh, F_nh, C_nh, P_nh, w1_nh = FSS_1993(w0_grid;h=0)
@assert sum(I_nh .< w1_nh) == 0 " I ≥ w must always be true !!!! "



pl_h_nh=plot(w0_grid[ind:end],H_opt[ind:end].*0, color=:black,legend=:none,yaxis=yaxis_h, title= "Hedging ratio, "*L"h=0"   ,xaxis=xaxis);
pl_w_nh=plot(w0_grid[ind:end],w1_nh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_w, title= "Cash-flow realizations, "*L"w"     ,xaxis=xaxis);
pl_e_nh=plot(w0_grid[ind:end],e_nh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(e_nh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
pl_I_nh=plot(w0_grid[ind:end],I_nh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(I_nh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
pl_F_nh=plot(w0_grid[ind:end],F_nh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(F_nh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
pl_P_nh=plot(w0_grid[ind:end],P_nh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w_0)"          ,xaxis=xaxis);
    # plot!(w0_grid[ind:end],mean(P_nh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w0)"          ,xaxis=xaxis);
pl_C_nh=plot(w0_grid[ind:end],C_nh[ind:end,:], palette=col,legend=:none,yaxis=[0 , .2], title= "Costs, "*L"C(e^*)"          ,xaxis=xaxis)
    # plot!(w0_grid[ind:end],mean(C_nh[ind:end,:],dims=2),  color=:red,legend=:none,yaxis=[0 , .2], title= "Costs, "*L"C(e^*)"          ,xaxis=xaxis);

summ_nh=plot(pl_h_nh,pl_w_nh,pl_I_nh,pl_e_nh,pl_F_nh,pl_P_nh,layout=(2,3),dpi=600 ,size=(900,600))



# ---------------------------------------------------------------------------- #
#                             contrafactual: h = 1                             #
# ---------------------------------------------------------------------------- #
@time H_fh, I_fh, e_fh, F_fh, C_fh, P_fh, w1_fh = FSS_1993(w0_grid;h=1)
@assert sum(I_fh .< w1_fh) == 0 " I ≥ w must always be true !!!! "


pl_h_fh=plot(w0_grid[ind:end],H_opt[ind:end].*0 .+1, color=:black,legend=:none,yaxis=yaxis_h, title= "Hedging ratio, "*L"h=1"     ,xaxis=xaxis);
pl_w_fh=plot(w0_grid[ind:end],w1_fh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_w, title= "Cash-flow realizations, "*L"w"     ,xaxis=xaxis);
pl_e_fh=plot(w0_grid[ind:end],e_fh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
        # plot!(w0_grid[ind:end],mean(e_fh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_e, title= "External financing, "*L"e^*",xaxis=xaxis);
pl_I_fh=plot(w0_grid[ind:end],I_fh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
        # plot!(w0_grid[ind:end],mean(I_fh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_I, title= "Investment, "*L"I^*"        ,xaxis=xaxis);
pl_F_fh=plot(w0_grid[ind:end],F_fh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
        # plot!(w0_grid[ind:end],mean(F_fh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_F, title= "Net value of I, "*L"F(I^*)" ,xaxis=xaxis);
pl_P_fh=plot(w0_grid[ind:end],P_fh[ind:end,:], palette=col,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w_0)"          ,xaxis=xaxis);
        # plot!(w0_grid[ind:end],mean(P_fh[ind:end,:],dims=2), color=:red,legend=:none,yaxis=yaxis_P, title= "Profit, "*L"P(w0)"          ,xaxis=xaxis);

pl_C_fh=plot(w0_grid[ind:end],C_fh[ind:end,:], palette=col,legend=:none,yaxis=[0 , .2], title= "Costs, "*L"C(e^*)"          ,xaxis=xaxis)

summ_fh=plot(pl_h_fh,pl_w_fh,pl_I_fh,pl_e_fh,pl_F_fh,pl_P_fh,layout=(2,3),dpi=600 ,size=(900,600))


sum(mean(P_fh[ind:end,:],dims=2) - mean(P_opt[ind:end,:],dims=2))

# ---------------------------------------------------------------------------- #
#                                    SUMMARY                                   #
# ---------------------------------------------------------------------------- #

summ
summ_nh
summ_fh

cd(dirname(@__FILE__))
savefig(summ,"froot_hast.png")
savefig(summ_nh,"froot_h0.png")
savefig(summ_fh,"froot_h1.png")