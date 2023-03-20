using LinearAlgebra
using Optim
using NLSolversBase
using Statistics

include("nll_lin.jl")
include("Gradp.jl")
include("Hessp.jl")

## a) simulate linear model

T   = 100                      #% set sample size
b0  = [0.2; 0.2; -0.1; 0.8; 0] #% set parameter values
sd  = 1.0                      #% set error standard deviation
x   = [ones(T,1) randn(T,4)]   #% simulate X matrix
err = randn(T,1)*sd            #% simulate error terms
y = x*b0 + err                 #% calculate y_i s

## b) estimate it by OLS
b=(x'*x)\(x'*y)
k = length(b)
y_hat    = x*b
err      = y - y_hat
sigma    = (err'*err)/(T-k)
var_cov  = inv(x' *x) .*  sigma
bse      = sqrt.(diag(var_cov))

## c) estimate ot by ML
datamat = [y x]                                         #% define data matrix for use in nll_lin
theta0  = [mean(y); zeros(size(b0,1)-1,1); std(y)]       #% this sets the initial parameter vector

LL(x) =  nll_lin(x,datamat)
res   = optimize(LL, theta0,LBFGS(),
                 Optim.Options(show_trace=true,
                 iterations=100_000,f_tol=1e-15))

thetaopt = res.minimizer

parameters = [thetaopt [b;sigma]]
## d) estimate different ML standard errors

LL2(x) =  nll_lin(x,datamat;vec=true)
g=gradp!(LL2,thetaopt)
J = g'g
se_J  = sqrt.(diag(inv(J)))

LL3(x) =  nll_lin(x,datamat;vec=false)
H = hessp!(LL3,thetaopt)
se_H  = sqrt.(diag(inv(H)))

se_SW = sqrt.(diag(inv(H)*J*inv(H)))    #% Sandwich variance covariance


## variance covariance estimators
[[bse;0] se_H se_J se_SW]