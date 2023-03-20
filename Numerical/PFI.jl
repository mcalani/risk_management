#------------------------------------------------------------------------------#
#  Policy function iteration 1
#------------------------------------------------------------------------------#
using Roots
using Interpolations
using Plots
using Optim

α=0.65;
β=0.95;
grid_max=2;
n=150;
N_iter=3000;
kgrid= collect(1e-2:(grid_max-1e-2)/(n-1):grid_max);
f(x)= x^α;
tol=1e-9;

#analytic solution
ab = α * β;
c1 = (log(1 - ab) + log(ab) * ab / (1 - ab)) / (1 - β);
c2 = α/ (1 - ab);
# optimal analytical values
v_star(k) = c1 .+ c2 .* log.(k);
k_star(k) = ab * k.^ α;
c_star(k) = (1-ab) * k.^α;
ufun(x) = log.(x);

#------------------------------------------------------------------------------#
function policy_iter(grid,c0,u_prime,f_prime)
    c1 = zeros(length(grid)) #next guess
    pol_fun = ( extrapolate(interpolate((grid,),
                          c0, Gridded(Linear()) ),
                          Interpolations.Flat()) ) #equivalent to clamp the argument
    #pol_fun = LinearInterpolation(grid,c0)
    #pol_fun = interpolate( (grid,), c0, Gridded(Linear()))

    #loop over current states of current capital
    for (i,k) in enumerate(grid)
        objective(c)= u_prime(c) - β*u_prime(pol_fun(f(k)-c))*f_prime(f(k)-c)
        c1[i]= fzero(objective,1e-10,f(k)-1e-10)
    end
    return c1
end

uprime(x)= 1.0 ./ x ;
fprime(x)= α*x.^(α-1);
c0 = collect(kgrid);
policy_iter(kgrid,c0,uprime,fprime)

#------------------------------------------------------------------------------#

function PFI()
    c_init = kgrid
    for iter in 1:N_iter
        c_next = policy_iter(kgrid,c_init,uprime,fprime)
        #check Convergence
        if maximum(abs,c_init.-c_next) < tol
            perrors = maximum(abs,c_next.-c_star(kgrid))
            println("PFI:")
            println("Found solution after $iter iterations")
            println("max policy function error = $perrors")
            return c_next
        elseif iter==N_iter
            warn("No solution found after $iter iterations")
            return c_next
        end
        c_init = c_next # update guess
     end
end

#------------------------------------------------------------------------------#

function plotPFI()
    v = PFI()
    plot(kgrid,[v v.-c_star(kgrid)],
    lab=["policy" "error"],
    legend=:bottomright,
    layout = 2)
end

plotPFI()

#------------------------------------------------------------------------------#
#revisiones
#example
x1=collect(1:1:1000)
x2=collect(0.1:0.1:100)
y1= log.(x1)
y2= x2 .^2
y=[y= x1[i]+x2[j] for i ∈ 1:length(x1), j ∈ 1:length(x2)]

#multiple interpolation
V0=interpolate((x1,x2),y, Gridded(Linear()))

minimum(x1)
maximum(x1)
minimum(x2)
maximum(x2)

x=34535345
y=45353453
x=clamp(x, minimum(x1),maximum(x1))
y=clamp(y, minimum(x2),maximum(x2))

V0(x,y)


#simple interpolation

v0=(extrapolate(interpolate(((kgrid),),c0, Gridded(Linear()) ),Interpolations.Flat()))
v1=LinearInterpolation(kgrid,c0)
v2=interpolate((kgrid,),c0, Gridded(Linear()))

v0(1)
v1(1)
v2(2)

v0==v1
v1==v2
v0==v2

typeof(v0)
typeof(v1)
typeof(v2)

findmin(v0)
findmin(v1)
findmin(v2)

minimum(kgrid)
maximum(kgrid)

x=2.1
v0(x)
v1(x)
v2(x)
