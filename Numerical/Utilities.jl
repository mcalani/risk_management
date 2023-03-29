## funciones auxiliares:
#max_by_element() construye una matriz a partir del maximo elemento entre cada
# par de valores de las dos matrices de entradas
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

#findnearest() encuentra el valor de un determinado vector que es más cercano a
# a determinado valor
function findnearest(v::AbstractArray, c::Real)
    findmin(abs.(v .- c ))[2]
end;

# interp1() interpola cada columna de la matrix Y sobre el vector X. luego evalua
# el punto "consulta" sobre cada interpolacion. Por lo tanto, el output es
# un vector de valores interpolados
query=0.5
x=copy(b_grid)
y=bp[s_shock[i],:]
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
            size(y)[1] == length(x) ? y = y' : nothing
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

#function unpack_all!()
function unpack_all!(parameters)
    tuple = parameters()
    names = keys(parameters())
    for i in 1:length(tuple)
        par_name = names[i]
        par_val  = tuple[names[i]]
        @eval $par_name = $par_val
    end
end

#Simulate AR(1)
function simulate_AR(μ, ρ, σ, n)
    dist = Normal(0,σ)
    y_sim = [0.0 for i = 1:n]
    noise = rand(dist, n)

    for i in 1:(n-1)
        y_sim[i+1] = (1-ρ)*μ + ρ*y_sim[i] + noise[i] 
    end
    return y_sim
end

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

function Simulate1!(n_sims,bp,bpsp,bpre;cut=Int(floor(n_sims*0.2)))

    s_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(s_shock))
    bpsp_sim = zeros(length(s_shock))
    bpre_sim = zeros(length(s_shock))

    @showprogress for i in 2:length(s_shock)
        #Simulate y-shocks
        YT_sim = YT[s_shock[i]]
        K_sim  = K[s_shock[i]]

        #Competitive Equilibrium Simulations
        bp_sim[i]   = interp1(b_grid, bp[s_shock[i],:]  , bp_sim[i-1]   )

        #Social Planner Simulations
        bpsp_sim[i] = interp1(b_grid, bpsp[s_shock[i],:], bpsp_sim[i-1] )

        #Reserves Accumulation
        bpre_sim[i] = interp1(b_grid, bpre[s_shock[i],:], bpre_sim[i-1] )
    end
    return bp_sim[cut:end], bpsp_sim[cut:end], bpre_sim[cut:end]
end


function Simulate3!(n_sims,bp,bpsp,bpre,
                    cpde,cpsp,cpre,
                    ppde,ppsp,ppre,apre ;cut=Int(floor(n_sims*0.2)))
    
    s_shock = mc_sample_path(T;sample_size = n_sims+cut)#Simulations
    #preallocations
    bp_sim   = zeros(length(s_shock))
    bpsp_sim = zeros(length(s_shock))
    bpre_sim = zeros(length(s_shock))

    @showprogress for i in 2:length(s_shock)
        #Simulate y-shocks
        YT_sim = YT[s_shock[i]]
        K_sim  = K[s_shock[i]]

        #Competitive Equilibrium Simulations
        bp_sim[i]   = interp1(b_grid, bp[s_shock[i],:]  , bp_sim[i-1]   )./(ppde - YT_sim )

        #Social Planner Simulations
        bpsp_sim[i] = interp1(b_grid, bpsp[s_shock[i],:], bpsp_sim[i-1] )

        #Reserves Accumulation
        bpre_sim[i] = interp1(b_grid, bpre[s_shock[i],:], bpre_sim[i-1] )
        end
    return bp_sim[cut:end], bpsp_sim[cut:end], bpre_sim[cut:end]
end

function Policy!(Pol)
    pol=reshape(Pol,length(κ_grid),length(y_grid),length(b_grid))
    pol_func(k,y,b)=extrapolate( interpolate((κ_grid,y_grid,b_grid), pol,
            Gridded(Linear())),Interpolations.Line() )[k,y,b]
    return pol_func
end

function Simulate2!(n_sims,
                    BPde,BPsp,BPre,
                    CPde,CPsp,CPre,
                    PPde,PPsp,PPre,APre; Out="Absolutes", cut=5000,
                    q95de=0.7316, q95sp=0.7042 , q95re=0.9373,
                    q75de=0.7179 , q75sp=0.6940, q75re= 0.8512)
    #simulate bonds, consumption, price
    @assert n_sims > cut
    global z_shock = simulate(MarkovChain(T),n_sims+cut) #Simulations
    #preallocations    
    for i in [:dbp_sim,  :dbpsp_sim,  :dbpre_sim,
              :dcp_sim, :dcpsp_sim, :dcpre_sim,
              :dpp_sim, :dppsp_sim, :dppre_sim,:da_sim,

            :bp_sim, :bpsp_sim, :bpre_sim,
            :cp_sim, :cpsp_sim, :cpre_sim,
            :pp_sim, :ppsp_sim, :ppre_sim,:a_sim,
            
            :gdp_bp, :gdp_bpsp, :gdp_bpre,
            :gdp_cp, :gdp_cpsp, :gdp_cpre,
            :gdp_pp, :gdp_ppsp, :gdp_ppre, :gdp_a]

            @eval $i = zeros(length(z_shock))
    end

    for i in [:CB_bp_sim, :CB_bpsp_sim, :CB_bpre_sim,
              :CB_cp_sim, :CB_cpsp_sim, :CB_cpre_sim,
              :CB_pp_sim, :CB_ppsp_sim, :CB_ppre_sim,
              :CB_a_sim,
              :CB_gdp_bp, :CB_gdp_bpsp, :CB_gdp_bpre,
              :CB_gdp_cp, :CB_gdp_cpsp, :CB_gdp_cpre,
              :CB_gdp_pp, :CB_gdp_ppsp, :CB_gdp_ppre,
              :CB_gdp_a,
              :CB_dbp, :CB_dbpsp, :CB_dbpre,
              :CB_dcp, :CB_dcpsp, :CB_dcpre,
              :CB_dpp, :CB_dppsp, :CB_dppre,
              :CB_da,
              :CB2_bp_sim, :CB2_bpsp_sim, :CB2_bpre_sim,
              :CB2_cp_sim, :CB2_cpsp_sim, :CB2_cpre_sim,
              :CB2_pp_sim, :CB2_ppsp_sim, :CB2_ppre_sim,
              :CB2_a_sim,
              :CB2_gdp_bp, :CB2_gdp_bpsp, :CB2_gdp_bpre,
              :CB2_gdp_cp, :CB2_gdp_cpsp, :CB2_gdp_cpre,
              :CB2_gdp_pp, :CB2_gdp_ppsp, :CB2_gdp_ppre,
              :CB2_gdp_a,
              :CB2_dbp, :CB2_dbpsp, :CB2_dbpre,
              :CB2_dcp, :CB2_dcpsp, :CB2_dcpre,
              :CB2_dpp, :CB2_dppsp, :CB2_dppre,
              :CB2_da]
        @eval $i = []
    end

    b   = mean(b_grid); bsp = mean(b_grid); bre = mean(b_grid); at  = mean(A_grid)

    BPde_func = Policy!(BPde); BPsp_func = Policy!(BPsp); BPre_func = Policy!(BPre)
    CPde_func = Policy!(CPde); CPsp_func = Policy!(CPsp); CPre_func = Policy!(CPre);
    PPde_func = Policy!(PPde); PPsp_func = Policy!(PPsp); PPre_func = Policy!(PPre)
    APre_func = Policy!(APre)

    @showprogress for i in 2:length(z_shock)
    #Simulate y-shocks
        yt = YT[z_shock[i]]
        κ  = K[z_shock[i]]

        ###Decentralized Equilibrium Simulations
        #absolutes
        bp_sim[i]    = BPde_func(κ,yt,b)
        cp_sim[i]    = CPde_func(κ,yt,b)               
        pp_sim[i]    = PPde_func(κ,yt,b)
        #differentials
        dcp_sim[i]   = (cp_sim[i] - cp_sim[i-1])/(pp_sim[i] +yt).*100
        dpp_sim[i]   = (pp_sim[i] - pp_sim[i-1])/(pp_sim[i] +yt).*100
        dbp_sim[i]   = (bp_sim[i] - b)/(pp_sim[i] +yt)          .*100
        #relative to gdp
        gdp_bp[i]    = bp_sim[i]/(pp_sim[i] +yt).*100
        gdp_cp[i]    = cp_sim[i]/(pp_sim[i] +yt).*100             
        gdp_pp[i]    = pp_sim[i]/(pp_sim[i] +yt).*100
        
        if (Out=="Credit Boom"|| Out=="All")  && bp_sim[i] > q95de #credit boom 95%
            global CB_bp_sim = push!(CB_bp_sim,bp_sim[i])
            global CB_cp_sim  = push!(CB_cp_sim,cp_sim[i])
            global CB_pp_sim = push!(CB_pp_sim,pp_sim[i])

            global CB_gdp_bp = push!(CB_gdp_bp,gdp_bp[i])
            global CB_gdp_cp = push!(CB_gdp_cp,gdp_cp[i])
            global CB_gdp_pp = push!(CB_gdp_pp,gdp_pp[i])

            global CB_dbp = push!(CB_dbp,dbp_sim[i])
            global CB_dcp = push!(CB_dcp,dcp_sim[i])
            global CB_dpp = push!(CB_dpp,dpp_sim[i])

            if i > 3 && (bp_sim[i] > q75de &&bp_sim[i-1] > q75de && bp_sim[i-2]> q75de)#credit boom 75% + 3 periods
                global CB2_bp_sim = push!(CB2_bp_sim,bp_sim[i])
                global CB2_cp_sim  = push!(CB2_cp_sim,cp_sim[i])
                global CB2_pp_sim = push!(CB2_pp_sim,pp_sim[i])
    
                global CB2_gdp_bp = push!(CB2_gdp_bp,gdp_bp[i])
                global CB2_gdp_cp = push!(CB2_gdp_cp,gdp_cp[i])
                global CB2_gdp_pp = push!(CB2_gdp_pp,gdp_pp[i])
    
                global CB2_dbp = push!(CB2_dbp,dbp_sim[i])
                global CB2_dcp = push!(CB2_dcp,dcp_sim[i])
                global CB2_dpp = push!(CB2_dpp,dpp_sim[i])
            end
        end

        ###Social Planner Simulations
        #absolutes
        bpsp_sim[i]  = BPsp_func(κ,yt,bsp) 
        cpsp_sim[i]  = CPsp_func(κ,yt,bsp)             
        ppsp_sim[i]  = PPsp_func(κ,yt,bsp)
        #differentials
        dcpsp_sim[i] = (cpsp_sim[i] - cpsp_sim[i-1])/(ppsp_sim[i] +yt).*100
        dppsp_sim[i] = (ppsp_sim[i] - ppsp_sim[i-1])/(ppsp_sim[i] +yt).*100
        dbpsp_sim[i] = (bpsp_sim[i]-bsp)            /(ppsp_sim[i] +yt).*100
        #relative to gdp
        gdp_bpsp[i]  = bpsp_sim[i] /(ppsp_sim[i] +yt).*100
        gdp_cpsp[i]  = cpsp_sim[i] /(ppsp_sim[i] +yt).*100             
        gdp_ppsp[i]  = ppsp_sim[i] /(ppsp_sim[i] +yt).*100 

        if (Out=="Credit Boom"|| Out=="All")  && bpsp_sim[i] > q95sp #credit boom 95%
            global CB_bpsp_sim = push!(CB_bpsp_sim,bpsp_sim[i])
            global CB_cpsp_sim = push!(CB_cpsp_sim,cpsp_sim[i])
            global CB_ppsp_sim = push!(CB_ppsp_sim,ppsp_sim[i])

            global CB_gdp_bpsp = push!(CB_gdp_bpsp,gdp_bpsp[i])
            global CB_gdp_cpsp = push!(CB_gdp_cpsp,gdp_cpsp[i])
            global CB_gdp_ppsp = push!(CB_gdp_ppsp,gdp_ppsp[i])

            global CB_dbpsp = push!(CB_dbpsp,dbpsp_sim[i])
            global CB_dcpsp = push!(CB_dcpsp,dcpsp_sim[i])
            global CB_dppsp = push!(CB_dppsp,dppsp_sim[i])
            
            if i > 3 && (bpsp_sim[i]>q75sp &&bpsp_sim[i-1]>q75sp&&bpsp_sim[i-2]> q75sp)#credit boom 75% + 3 periods
                global CB2_bpsp_sim = push!(CB2_bpsp_sim,bpsp_sim[i])
                global CB2_cpsp_sim = push!(CB2_cpsp_sim,cpsp_sim[i])
                global CB2_ppsp_sim = push!(CB2_ppsp_sim,ppsp_sim[i])
    
                global CB2_gdp_bpsp = push!(CB2_gdp_bpsp,gdp_bpsp[i])
                global CB2_gdp_cpsp = push!(CB2_gdp_cpsp,gdp_cpsp[i])
                global CB2_gdp_ppsp = push!(CB2_gdp_ppsp,gdp_ppsp[i])
    
                global CB2_dbpsp = push!(CB2_dbpsp,dbpsp_sim[i])
                global CB2_dcpsp = push!(CB2_dcpsp,dcpsp_sim[i])
                global CB2_dppsp = push!(CB2_dppsp,dppsp_sim[i])
            end
        end

        ### RR.II Economy Simulations
        #absolutes
        cpre_sim[i]  = CPre_func(κ,yt, bre - at)  #cpre(k,yt,bre,at)
        ppre_sim[i]  = PPre_func(κ,yt, bre - at)  #ppre(k,yt,bre,at)
        a_sim[i]     = APre_func(κ,yt, bre - at)  #apre(k,yt,bre,at)
        bpre_sim[i]  = BPre_func(κ,yt, bre - at)  #bsp' = bre' - at'

        #differentials
        dcpre_sim[i] = (cpre_sim[i] - cpre_sim[i-1])/(ppre_sim[i] +yt).*100
        dppre_sim[i] = (ppre_sim[i] - ppre_sim[i-1])/(ppre_sim[i] +yt).*100
        dbpre_sim[i] = (bpre_sim[i]-bre)            /(ppre_sim[i] +yt).*100
        da_sim[i]    = (a_sim[i] - at)              /(ppre_sim[i] +yt).*100
        
        #relative to gdp
        gdp_bpre[i]  = bpre_sim[i] /(ppre_sim[i] +yt).*100
        gdp_cpre[i]  = cpre_sim[i] /(ppre_sim[i] +yt).*100        
        gdp_ppre[i]  = ppre_sim[i] /(ppre_sim[i] +yt).*100
        gdp_a[i]     = a_sim[i]    /(ppre_sim[i] +yt).*100

        if (Out=="Credit Boom"|| Out=="All")  && bpre_sim[i] > q95re   #credit boom 95%
            global CB_bpre_sim = push!(CB_bpre_sim,bpre_sim[i])
            global CB_cpre_sim = push!(CB_cpre_sim,cpre_sim[i])
            global CB_ppre_sim = push!(CB_ppre_sim,ppre_sim[i])
            global CB_a_sim    = push!(CB_a_sim,   a_sim[i])

            global CB_gdp_bpre = push!(CB_gdp_bpre,gdp_bpre[i])
            global CB_gdp_cpre = push!(CB_gdp_cpre,gdp_cpre[i])
            global CB_gdp_ppre = push!(CB_gdp_ppre,gdp_ppre[i])
            global CB_gdp_a    = push!(CB_gdp_a,gdp_a[i])

            global CB_dbpre = push!(CB_dbpre,dbpre_sim[i])
            global CB_dcpre = push!(CB_dcpre,dcpre_sim[i])
            global CB_dppre = push!(CB_dppre,dppre_sim[i])
            global CB_da  = push!(CB_da,da_sim[i])

            if i > 3 && (bpre_sim[i] > q75re  &&bpre_sim[i-1]>q75re && bpre_sim[i-2]> q75re) #credit boom 75% + 3 periods
                global CB2_bpre_sim = push!(CB2_bpre_sim,bpre_sim[i])
                global CB2_cpre_sim = push!(CB2_cpre_sim,cpre_sim[i])
                global CB2_ppre_sim = push!(CB2_ppre_sim,ppre_sim[i])
                global CB2_a_sim    = push!(CB2_a_sim,   a_sim[i])
    
                global CB2_gdp_bpre = push!(CB2_gdp_bpre,gdp_bpre[i])
                global CB2_gdp_cpre = push!(CB2_gdp_cpre,gdp_cpre[i])
                global CB2_gdp_ppre = push!(CB2_gdp_ppre,gdp_ppre[i])
                global CB2_gdp_a    = push!(CB2_gdp_a,gdp_a[i])
    
                global CB2_dbpre = push!(CB2_dbpre,dbpre_sim[i])
                global CB2_dcpre = push!(CB2_dcpre,dcpre_sim[i])
                global CB2_dppre = push!(CB2_dppre,dppre_sim[i])
                global CB2_da  = push!(CB2_da,da_sim[i])
            end
        end

        #update state variables
        at  = a_sim[i]
        b   = bp_sim[i]
        bsp = bpsp_sim[i] # = bre - at 
        bre = bpre_sim[i]
        
    end
    #effect of shocks on consumption and prices
    if Out == "Diff"
        return [cc[cut:end],  ccSP[cut:end],  ccRE[cut:end],
                dcp_sim[cut:end], dcpsp_sim[cut:end], dcpre_sim[cut:end],
                dpp_sim[cut:end], dppsp_sim[cut:end], dppre_sim[cut:end],
                da_sim[cut:end]]

    elseif Out == "Absolutes" 
        return [bp_sim[cut:end], bpsp_sim[cut:end], bpre_sim[cut:end],
                cp_sim[cut:end], cpsp_sim[cut:end], cpre_sim[cut:end],
                pp_sim[cut:end], ppsp_sim[cut:end], ppre_sim[cut:end],
                a_sim[cut:end]]

    elseif Out == "Credit Boom" 
        return [CB_bp_sim, CB_bpsp_sim, CB_bpre_sim,
                CB_cp_sim,  CB_cpsp_sim, CB_cpre_sim,
                CB_pp_sim,  CB_ppsp_sim, CB_ppre_sim,
                CB_a_sim,
                
                CB_gdp_bp, CB_gdp_bpsp, CB_gdp_bpre,
                CB_gdp_cp, CB_gdp_cpsp, CB_gdp_cpre,
                CB_gdp_pp, CB_gdp_ppsp,CB_gdp_ppre,
                CB_gdp_a,
                
                CB_dbp, CB_dbpsp, CB_dbpre,
                CB_dcp, CB_dcpsp, CB_dcpre,
                CB_dpp, CB_dppsp, CB_dppre,
                CB_da,
                
                CB2_bp_sim, CB2_bpsp_sim, CB2_bpre_sim,
                CB2_cp_sim,  CB2_cpsp_sim, CB2_cpre_sim,
                CB2_pp_sim,  CB2_ppsp_sim, CB2_ppre_sim,
                CB2_a_sim,
                
                CB2_gdp_bp, CB2_gdp_bpsp, CB2_gdp_bpre,
                CB2_gdp_cp, CB2_gdp_cpsp, CB2_gdp_cpre,
                CB2_gdp_pp, CB2_gdp_ppsp, CB2_gdp_ppre,
                CB2_gdp_a,
                
                CB2_dbp, CB2_dbpsp, CB2_dbpre,
                CB2_dcp, CB2_dcpsp, CB2_dcpre,
                CB2_dpp, CB2_dppsp, CB2_dppre,
                CB2_da]

    elseif Out == "GDP"
        return [gdp_bp[cut:end], gdp_bpsp[cut:end], gdp_bpre[cut:end],
                gdp_cp[cut:end], gdp_cpsp[cut:end], gdp_cpre[cut:end],
                gdp_pp[cut:end], gdp_ppsp[cut:end], gdp_ppre[cut:end],
                gdp_a[cut:end]]

    elseif Out == "All"
        return [dbp_sim[cut:end], dbpsp_sim[cut:end], dbpre_sim[cut:end],
                dcp_sim[cut:end], dcpsp_sim[cut:end], dcpre_sim[cut:end],
                dpp_sim[cut:end], dppsp_sim[cut:end], dppre_sim[cut:end],
                da_sim[cut:end],

                bp_sim[cut:end], bpsp_sim[cut:end], bpre_sim[cut:end],
                cp_sim[cut:end], cpsp_sim[cut:end], cpre_sim[cut:end],
                pp_sim[cut:end], ppsp_sim[cut:end], ppre_sim[cut:end],
                a_sim[cut:end],
                
                gdp_bp[cut:end], gdp_bpsp[cut:end], gdp_bpre[cut:end],
                gdp_cp[cut:end], gdp_cpsp[cut:end], gdp_cpre[cut:end],
                gdp_pp[cut:end], gdp_ppsp[cut:end], gdp_ppre[cut:end],
                gdp_a[cut:end],
                
                CB_bp_sim, CB_bpsp_sim, CB_bpre_sim, #credit boom 95%
                CB_cp_sim,  CB_cpsp_sim, CB_cpre_sim,
                CB_pp_sim,  CB_ppsp_sim, CB_ppre_sim,
                CB_a_sim,
                
                CB_gdp_bp, CB_gdp_bpsp, CB_gdp_bpre,
                CB_gdp_cp, CB_gdp_cpsp, CB_gdp_cpre,
                CB_gdp_pp, CB_gdp_ppsp,CB_gdp_ppre,
                CB_gdp_a,
                
                CB_dbp, CB_dbpsp, CB_dbpre,
                CB_dcp, CB_dcpsp, CB_dcpre,
                CB_dpp, CB_dppsp, CB_dppre,
                CB_da,
                
                CB2_bp_sim, CB2_bpsp_sim, CB2_bpre_sim, #credit boom 75% + 3 periods
                CB2_cp_sim,  CB2_cpsp_sim, CB2_cpre_sim,
                CB2_pp_sim,  CB2_ppsp_sim, CB2_ppre_sim,
                CB2_a_sim,
                
                CB2_gdp_bp, CB2_gdp_bpsp, CB2_gdp_bpre,
                CB2_gdp_cp, CB2_gdp_cpsp, CB2_gdp_cpre,
                CB2_gdp_pp, CB2_gdp_ppsp, CB2_gdp_ppre,
                CB2_gdp_a,
                
                CB2_dbp, CB2_dbpsp, CB2_dbpre,
                CB2_dcp, CB2_dcpsp, CB2_dcpre,
                CB2_dpp, CB2_dppsp, CB2_dppre,
                CB2_da]
    end
end

# v'T = v ; where v is the stationary distribution
function sdist(T ; method="Eigenvector")
    rows =[ sum(T[i,:])  for i in 1:size(T,1)]
    @assert size(T,1) == size(T,2)      #squared matrix
    @assert sum(rows .≈ 1) == size(T,1) #rows must sum 1
    if method == "Eigenvector"
        λ_pos=findnearest(real(eigvals(T)) , 1.)
        #Method1: eigenvalues
        λ     = real(eigvals(T))[λ_pos]
        @assert λ≈1
        Reig  = eigvecs(T) #right eigenvalues
        Leig  = inv(Reig)  #left  eigenvalues (stationary distribution)
        v     = vec(real.(Leig[end,:]))
        #scaling v
        edist1=vec( v'T./sum(v'T)  )
        return edist1
    elseif method == "Iteration"
        v=ones(size(T,1))./size(T,1)
        count = 0
        while count < 100_000
            v1 = vec(v'T) # left eigenvalue associated with eigenvalue=1
            dif=norm(v1-v)
            println(dif)
            v=copy(v1)
            if vec(v'T) ≈ v
                break
            end
            count += 1
        end
        #scaling v
        edist2=v./sum(v)
        return edist2 
    end
end