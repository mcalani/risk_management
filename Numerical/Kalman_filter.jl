#------------------------------------------------------------------------------#
# QuantEcon Cap.21: A first look at the Kalman filter
#------------------------------------------------------------------------------#
# Rocket lauch example
#------------------------------------------------------------------------------#

using InstantiateFromURL
github_project("QuantEcon/quantecon-notebooks-julia", version="0.5.0")
using LinearAlgebra, Statistics, Plots, Distributions
#------------------------------------------------------------------------------#
# I.Prior information
#------------------------------------------------------------------------------#

Σ= [0.4 0.3
    0.3 0.45];
x̂= [ 0.2, -0.2];

#Plot
x_grid= range(-1.5, 2.9, length= 100);
y_grid= range(-3.1, 1.7, length= 100);

#generate distribuions
dist= MvNormal(x̂,Σ);
two_args_to_pdf(dist)= (x,y) -> pdf(dist, [x,y]); #returns a function to be plotted

#plot
contour(x_grid, y_grid, two_args_to_pdf(dist), fill=false,
color= :lighttest, cbar= false)
#El centro de la elipse es x̂

#------------------------------------------------------------------------------#
# II.New information
#------------------------------------------------------------------------------#
# Asumimos G y R como dados

#Define G and R from the equation y= Gx + N(0,R)
y=[2.3, -1.9];
G=I;
R=0.5 .* Σ;

#------------------------------------------------------------------------------#
# III. Filtering
#------------------------------------------------------------------------------#
# update information via Bayes: p(x|y)= P(y|x)p(x)/p(y)
# p(x)~ N(x̂, Σ)
# p(y|x)~ N(Gx,R)
# p(y) does not depend on x

#Define posterior objects
M= Σ*G' * inv(G * Σ * G' + R);
x̂_F= x̂+ M * (y-G*x̂);
Σ_F= Σ -M * G * Σ;

#Plot
newdist= MvNormal(x̂_F, Symmetric(Σ_F));
contour(x_grid, y_grid, two_args_to_pdf(dist), fill=false,
color= :lighttest,levels=10, cbar= false)
contour!(x_grid, y_grid, two_args_to_pdf(newdist), fill=false,
color= :lighttest, cbar= false)
annotate!(y[1], y[2], "y", color= :black)

# Hemos obtenido las probabilidades de la ubicación actual,
# dadas las priors y la current information
# p(x|y)=N(x̂_F, Σ_F) is called the filtering distribution

#------------------------------------------------------------------------------#
# IV. Forecast
#------------------------------------------------------------------------------#
# we need a model of how the states evolves: suppose we have one
# x_t+1 = Ax_t + w_t+1, donde w_t ~ N(0,Q)
# we want to combine this model with a new predictive distribution
# all we have to do is introduce x_F~N(x̂_F, Σ_F) and work out
# the distribution of Ax_F + w

A=[ 1.2 0.0
    0.0 -0.2];

Q= 0.3*Σ;

#get the predictive distribution

new_x̂ = A* x̂_F ;
new_Σ = A* Σ_F * A' + Q;

predictdist = MvNormal(new_x̂, Symmetric(new_Σ));

#plot density
contour(x_grid, y_grid, two_args_to_pdf(predictdist), fill=false,
color= :lighttest,lw=1, cbar= false)
contour!(x_grid, y_grid, two_args_to_pdf(dist), fill=false,
color= :grays, cbar= false)
contour!(x_grid, y_grid, two_args_to_pdf(newdist), fill=false,
color= :grays, cbar= false)

annotate!(y[1], y[2], "y", color= :black)
annotate!(-0.5, 0.5, "x", color= :black)
