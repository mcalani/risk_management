
function interp1(x::Vector{Float64},
                y::AbstractArray{Float64},
                query::Float64)
    @assert length(x) == size(y,2) || length(x) == size(y,1)

    vector= zeros(size(y,1))
    scalar= 0.0
    if ndims(y)==2
        for i in 1:size(y,1)
            y_intp=extrapolate( interpolate((x,), y[i,:],
                        Gridded(Linear())),Interpolations.Line() )
            vector[i]= y_intp(query)
        end
        return vector
    elseif ndims(y)==1
        y_intp=extrapolate( interpolate((x,), y,
                Gridded(Linear())),Interpolations.Line() )
        scalar = y_intp(query)
        return scalar
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

function min_by_element!(A,B)
        @assert size(A)==size(B) && ndims(A)==2
        C=zeros(size(A,1),size(A,2));
    for i ∈ 1:size(C,1)
        for j ∈ 1:size(C,2)
            C[i,j] = minimum([A[i,j],B[i,j]])
        end
    end
    return C
end;

function mc_sample_path(P; init = 1, sample_size = 1000)
    @assert size(P)[1] == size(P)[2] # square required
    N = size(P)[1] # should be square

    # create vector of discrete RVs for each row
    dists = [Categorical(P[i, :]) for i in 1:N]

    # setup the simulation
    X = fill(0, sample_size) # allocate memory, or zeros(Int64, sample_size)
    X[1] = init # set the initial state

    for t in 2:sample_size
        dist = dists[X[t-1]] # get discrete RV from last state's transition distribution
        X[t] = rand(dist) # draw new value
    end
    return X
end

function Simulate!(n_sims,bp,bpsp;cut=Int(floor(n_sims*0.2)))
    y_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(y_shock))
    bpsp_sim = zeros(length(y_shock))

    @showprogress for i in 2:length(y_shock)
        #Simulate y-shocks
        YT_sim=YT[y_shock[i]]
        YN_sim=YN[y_shock[i]]

        #Competitive Equilibrium Simulations
        bp_sim[i]   = interp1(b_grid, bp[y_shock[i],:]  , bp_sim[i-1]   )

#       #Social Planner Simulations
        bpsp_sim[i] = interp1(b_grid, bpsp[y_shock[i],:], bpsp_sim[i-1] )

    end
    return bp_sim[cut:end], bpsp_sim[cut:end]
end


function findnearest(v::AbstractArray, c::Real)
    findmin(abs.(v .- c ))[2]
end;

function Simulate2!(n_simulation,bp,bpsp;cut=Int(floor(n_simulation*0.2)))
    y_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(y_shock))
    bpsp_sim = zeros(length(y_shock))

    @showprogress for i in 2:length(y_shock)
        #Simulate y-shocks
        YT_sim=YT[y_shock[i]]
        YN_sim=YN[y_shock[i]]
        #Competitive Equilibrium Simulations
        bp_sim[i] = bp[findnearest(bp[y_shock[i],:],bp_sim[i-1])]

        #Social Planner Simulations
        bpsp_sim[i] = bpsp[findnearest(bpsp[y_shock[i],:],bpsp_sim[i-1])]

    end
    return bp_sim[cut:end], bpsp_sim[cut:end]
end