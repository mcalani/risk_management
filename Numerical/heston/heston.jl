# ---------------------------------------------------------------------------- #
#                               Heston Portfolio                               #
# ---------------------------------------------------------------------------- #

#Boguslavskaya2015: CIR
#Sandjo2017: O-U process

using Random
using Statistics
using Distributions
using LinearAlgebra
using Plots


#VaR constraint

u(W;K=5,γ=0.5,λ=1) = (1/γ).*W.^γ .- λ*[if (W.<K) 1 else 0 end ][] 
u2(W;γ=0.5) = (1/γ).*W.^γ

u3(W;K=5,γ=0.5) = (1/γ).*(W).^γ .*(1+[if (W.<K) 1 .-W./K else 0 end ][]) 


w_range= 00.1:0.01:10

plot(w_range,[u.(w_range) u2.(w_range) u3.(w_range)])


#deep parameters
μ   = 0.02
σ   = 0.2
r_f = 0.01
γ   = 2
u(w,γ)= w^γ


#grid for wealth: dW_t = {[α_t(μ-r) + r]*W_t}*dt + {(1/2)*α*σ*W_t}*dB_t

Nw      = 800
w_min   = 0.01
w_max   = 2
w_grid  = collect(range(w_min,w_max,length=Nw))

dw = (w_max- w_min)/Nw

#grid for time: uniform
t_min=0
t_max=1
Nt = 800

dt = (t_max-t_min)/Nt





#with u(x) = x^γ
α_opt(μ,r_f,γ) = - (μ - r_f)/σ^2  * γ/(1-γ)