# ---------------------------------------------------------------------------- #
#                                 Arellano 2008                                #
# ---------------------------------------------------------------------------- #
using Parameters, Plots, QuantEcon
using LinearAlgebra, Statistics
using Optim
using Interpolations
using ProgressMeter
using BenchmarkTools
# ---------------------------------------------------------------------------- #
#                              Auxiliary Functions                             #
# ---------------------------------------------------------------------------- #

function interp1(x::Vector{Float64},
                 y::AbstractArray{Float64},
                 query::Float64;method="Linear")
    @assert length(x) == size(y,2) || length(x) == size(y,1) "Vector must have the same length"

    x_min = minimum(x)
    x_max = maximum(x)
    nx    = length(x)
    if method == "Cubic" || method == "Quadratic"
        @assert x==collect(range(x_min,x_max,length=nx)) "Grid points are needed to be uniformly distributed to use Quadratic or Cubic method. Use Linear instead"
    end
    
    vector= zeros(size(y,1))
    scalar= 0.0
    if method=="Linear"
        if ndims(y)==2
            for i in 1:size(y,1)
                y_intp      =interpolate((x,), y[i,:], Gridded(Linear()))
                y_intp_extrp=extrapolate( y_intp ,Interpolations.Line() )
                vector[i]   = y_intp_extrp(query)
            end
            return vector
        elseif ndims(y)==1
            y_intp      =interpolate((x,), y, Gridded(Linear()))
            y_intp_extrp=extrapolate( y_intp ,Interpolations.Line() )
            scalar      = y_intp_extrp(query)
            return scalar
        end
    elseif method == "Quadratic"
        if ndims(y)==2
            for i in 1:size(y,1)
                y_intp        = interpolate(y[i,:], BSpline(Quadratic(Line(OnCell()))))
                y_intp_scaled = Interpolations.scale(y_intp,range(x_min,x_max,length=nx))
                y_intp_extrp  = extrapolate( y_intp_scaled ,Interpolations.Line() )
                vector[i]     = y_intp_extrp(query)
            end
            return vector
        elseif ndims(y)==1
            y_intp        = interpolate(y, BSpline(Quadratic(Line(OnCell()))))
            y_intp_scaled = Interpolations.scale(y_intp,range(x_min,x_max,length=nx))
            y_intp_extrp  = extrapolate( y_intp_scaled ,Interpolations.Line() )
            scalar      = y_intp_extrp(query)
            return scalar
        end
    elseif method=="Cubic"
        if ndims(y)==2
            for i in 1:size(y,1)
                y_intp        = interpolate(y[i,:], BSpline(Cubic(Line(OnGrid()))))
                y_intp_scaled = Interpolations.scale(y_intp,range(x_min,x_max,length=nx))
                y_intp_extrp  = extrapolate( y_intp_scaled ,Interpolations.Line() )
                vector[i]     = y_intp_extrp(query)
            end
            return vector
        elseif ndims(y)==1
            y_intp        = interpolate(y, BSpline(Cubic(Line(OnGrid()))))
            y_intp_scaled = Interpolations.scale(y_intp,range(x_min,x_max,length=nx))
            y_intp_extrp  = extrapolate( y_intp_scaled ,Interpolations.Line() )
            scalar      = y_intp_extrp(query)
            return scalar
        end
    end
end;

function max_by_element!(A,B)
    @assert size(A)==size(B) && ndims(A)==2
    C=zeros(size(A,1),size(A,2));
    for i ∈ 1:size(C,1)
        for j ∈ 1:size(C,2)
            C[i,j] = maximum([A[i,j],B[i,j]])
        end
    end
    return C
end;

#Parameters
const β   = 0.953
const γ   = 2.0
const r   = 0.017
const ρ   = 0.945
const η   = 0.025
const θ   = 0.282 #reentering international credit market probability
const ny  = 10#21
const nB  = 100 #251

#Discretization
Π        = tauchen(ny, ρ, η).p
y_grid   = exp.(tauchen(ny, ρ, η).state_values)
ydefgrid = min.(.969 * mean(y_grid), y_grid)
B_grid   = collect(range(-.4, .4, length = nB))

#Initial values
vf = zeros(nB, ny)
vd = zeros(1, ny)
vc = zeros(nB, ny)
policy = zeros(nB, ny)
q = ones(nB, ny) .* (1 / (1 + r))
defprob = zeros(nB, ny)

function u(c;γ=2)
    if c > 0.0
        γ==1 ? log(c) : (c^(1-γ) -1)/(1-γ)
    else
        -Inf
    end
end

# ---------------------------------------------------------------------------- #
#                               Bellman Operator                               #
# ---------------------------------------------------------------------------- #
function T_operator(vd,vc,vf,q)
    Tvd= zeros(1, ny)
    Tvc= zeros(nB, ny)
    Tvf= zeros(nB, ny)
    defs= zeros(nB, ny)
    Bp= zeros(nB, ny)
    for j in 1:ny
            #Default Value
            Tvd[1,j] = u(ydefgrid[j]) + β*Π[j,:]'*interp1(B_grid,(θ*vf .+ (1-θ)*vd)', 0.)

            for i in 1:nB
                #Paying debt Value
                RHSc(Bp) = - (u(y_grid[j] + B_grid[i] - q[i,j]*Bp) + β*Π[j,:]'*interp1(B_grid,vc',Bp) )
                Opt      = optimize(RHSc, minimum(B_grid) ,y_grid[j] + B_grid[j])
                Tvc[i,j] = -Opt.minimum
                Bp[i,j]    = Opt.minimizer

                #Value Function
                Tvf[i,j]   = maximum([Tvc[i,j] Tvd[1,j]])
            end
    end
    return Tvd, Tvc, Tvf, Bp
end


# #test
#@time Tvd1, Tvc1, Tvf1, Bp1 = T_operator(vd,vc,vf,q)
#Tvd, Tvc, Tvf,Bp, defs = T_operator(vd,vc,vf,q)
# plot(B_grid,Tvc)
# plot(B_grid,Tvf)
# plot(Tvd')

# ---------------------------------------------------------------------------- #
#                                      VFI                                     #
# ---------------------------------------------------------------------------- #

function VFI(vd,vc,vf,q,T_operator;show_every=10)
    iter  = 1
    N_iter= 100
    conv  = 100
    tol   = 1e-5
    Tvd= zeros(1, ny)
    Tvc= zeros(nB, ny)
    Tvf= zeros(nB, ny)
    defs= zeros(nB, ny)
    #δ =zeros(nB, ny)
    Bp= zeros(nB, ny)
    for i in 1:N_iter
        #Bellman operator
        Tvd, Tvc, Tvf, Bp = T_operator(vd,vc,vf,q)        
        #converce criterion
        conv  = norm( Tvf - vf)

        if iter % show_every == 0
            println("Diff = $(round(conv,digits=10)) ; Iter= $iter")
        end
        if conv < tol
            break
            println("Convergence achieved :D !")
        end

        defs = repeat(Tvd,nB) .> Tvc
        # for i in 1:nB, j in 1:ny
        #     δ[i,j] = Π[j,:]'*interp1(B_grid,defs', Bp[i,j])
        # end
        δ = defs*Π'
        q = (1 .-δ)./(1+r)
        #loop
        vd=copy(Tvd)
        vc=copy(Tvc)
        vf=copy(Tvf)

        iter += 1
    end
    return vd,vc,vf, Bp, δ, q
end

#test
@time vd,vc,vf, Bp, δ, q = VFI(vd,vc,vf,q,T_operator;show_every=1)


defs== default_states

#opt1
sum((default_states*Π')[1,:])

#opt2
for i in 1:nB, j in 1:ny
    δ[i,j] = Π[j,:]'*interp1(B_grid,defs', Bp[i,j])
end

sum(δ[1,:])


# ---------------------------------------------------------------------------- #
#                          Solving Arellano's Economy                          #
# ---------------------------------------------------------------------------- #

function Arellano(vd,vc,vf,q)
    iter  = 1
    N_iter= 500
    conv  = 100
    tol   = 1e-5
    δ     = zeros(nB, ny)
    for m in 1:N_iter
        vd,vc,vf, Bp, defs = VFI(vd,vc,vf,q,T_operator)

        for i in 1:nB, j in 1:ny
            δ[i,j] = Π[j,:]'*interp1(B_grid,defs',Bp[i,j])
        end
        q_next = (1 .-δ)./(1+r)
        #converce criterion
        conv  = norm( q - q_next)
        println("q(B',y) conv= $(round(conv,digits=9)) ; Iter= $iter")
        if conv < tol
            break
            println("Convergence achieved :D !")
        end

        #Update price and repeat
        q = copy(q_next)
        iter += 1
    end
    return vd,vc,vf, Bp, q, δ
end

vd,vc,vf, Bp, q, δ = Arellano(vd,vc,vf,q)


