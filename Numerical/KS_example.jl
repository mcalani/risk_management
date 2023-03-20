#-------------------------------------------------------------------------------
# KS functions
#-------------------------------------------------------------------------------
using Interpolations # to use interpolation
using Random, LinearAlgebra
using QuantEcon  # to use `gridmake`, `<:AbstractUtility`
using Optim      # to use minimization routine to maximize RHS of bellman equation
using GLM        # to regress
using JLD2       # to save the result
using ProgressMeter # to show progress of iterations
using Parameters # to use type with keyword arguments

include("parameters.jl")
@unpack γ,β,c_b,c_f,ζ_b,ζ_f,κ_b,κ_f,φ_bz,φ_fz,x,ξ,n_grid = Incumbents();
@unpack F_z,F_z_stat,z,G_δb, G_δf,δ_b,δ_f = Stochastic_processes();


#Transition matrix
struct TransitionMatrix
    P::Matrix{Float64}       # 4x4
    Pz::Matrix{Float64}      # 2x2 aggregate shock
    Peps_gg::Matrix{Float64} # 2x2 idiosyncratic shock conditional on good to good
    Peps_bb::Matrix{Float64} # 2x2 idiosyncratic shock conditional on bad to bad
    Peps_gb::Matrix{Float64} # 2x2 idiosyncratic shock conditional on good to bad
    Peps_bg::Matrix{Float64} # 2x2 idiosyncratic shock conditional on bad to good
end

abstract type UMPSolutionMethod end

@with_kw struct VFI <: UMPSolutionMethod
    Howard_on::Bool = false
    Howard_n_iter::Int = 20
end
