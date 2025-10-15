module DCPF_H2

using JuMP, Gurobi
using CSV, DataFrames
import MathOptInterface as MOI
import ..Data: NODES, ZONES, GENERATORS, GEN_COST, E_DEMAND, zone_of,
               AC_LINES, DC_LINES, H2_DEMAND, H2_CONVERSION_COST, H2_EFFICIENCY, H2_CONNECTIONS, H2_TRANS_LIMIT,
               ONSHORE_NODES, OFFSHORE_NODES, REP_HOURS, HOUR_WEIGHT,
               H2_CAPEX_PER_MW, H2_CAP_LIMIT, H2_ANNUITY_FACTOR

export dcpf_H2

function dcpf_H2()
    m = Model(Gurobi.Optimizer)
    set_silent(m)

    # Electricity generation
    @variable(m, 0 <= g[h in REP_HOURS, n in NODES] <= GENERATORS[(h,n)])  
    @variable(m, θ[h in REP_HOURS, n in NODES])                            
    # Hydrogen system
    @variable(m, 0 <= h2_cap[n in NODES] <= H2_CAP_LIMIT[n])              
    @variable(m, 0 <= h2_prod[h in REP_HOURS, n in NODES])                
    @variable(m, -H2_TRANS_LIMIT[(i,j)] <= h2_trans[h in REP_HOURS, (i,j) in H2_CONNECTIONS] <= H2_TRANS_LIMIT[(i,j)])
    h2_cap_bind = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for h in REP_HOURS, n in NODES
        h2_cap_bind[(h,n)] = @constraint(m, h2_prod[h,n] <= h2_cap[n])
    end

    h2_balance = Dict{Tuple{Symbol, Symbol}, ConstraintRef}()
    for h in REP_HOURS, n in NODES
        h2_in  = sum(h2_trans[h,(i,n)] for (i,j) in H2_CONNECTIONS if j == n; init=0.0)
        h2_out = sum(h2_trans[h,(n,j)] for (i,j) in H2_CONNECTIONS if i == n; init=0.0)
        h2_balance[(h,n)] = @constraint(m, h2_prod[h,n] + h2_in - h2_out == H2_DEMAND[(h,n)])
    end

    for h in REP_HOURS, (i,j) in H2_CONNECTIONS
        @constraint(m, h2_trans[h,(i,j)] <= H2_TRANS_LIMIT[(i,j)])
    end

    # Fix one angle reference
    adj = Dict(n => Symbol[] for n in NODES)
    for (i,j,_,_) in AC_LINES
        push!(adj[i], j); push!(adj[j], i)
    end
    seen = Set{Symbol}()
    for n in NODES
        if !(n in seen)
            comp = Symbol[]; stack = [n]; push!(seen, n)
            while !isempty(stack)
                u = pop!(stack); push!(comp, u)
                for v2 in adj[u]
                    if !(v2 in seen); push!(seen, v2); push!(stack, v2); end
                end
            end
            @constraint(m, [h in REP_HOURS], θ[h, comp[1]] == 0)
        end
    end

    # AC transmission 
    @expression(m, f_ac[h in REP_HOURS, ℓ in 1:length(AC_LINES)],
        (1 / AC_LINES[ℓ][4]) * (θ[h, AC_LINES[ℓ][1]] - θ[h, AC_LINES[ℓ][2]])
    )
    @variable(m, t_ac[h in REP_HOURS, ℓ in 1:length(AC_LINES)] >= 0)  

    for (ℓ, (_, _, F, _)) in enumerate(AC_LINES)
        @constraint(m, [h in REP_HOURS], -F <= f_ac[h,ℓ] <= F)
        @constraint(m, [h in REP_HOURS], t_ac[h,ℓ] >=  f_ac[h,ℓ])
        @constraint(m, [h in REP_HOURS], t_ac[h,ℓ] >= -f_ac[h,ℓ])
    end

    # DC transmission
    @variable(m, f_dc[h in REP_HOURS, ℓ in 1:length(DC_LINES)])

    for (ℓ, (_, _, F)) in enumerate(DC_LINES)
        @constraint(m, [h in REP_HOURS], -F <= f_dc[h,ℓ] <= F)
    end

    # Net DC inflow 
    dc_out = Dict(ℓ => DC_LINES[ℓ][1] for ℓ in eachindex(DC_LINES))
    dc_in  = Dict(ℓ => DC_LINES[ℓ][2] for ℓ in eachindex(DC_LINES))
    @expression(m, dc_net[h in REP_HOURS, n in NODES],
        sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if dc_in[ℓ] == n) -
        sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if dc_out[ℓ] == n)
    )

    # Nodal balance
    @constraint(m, node_balance[h in REP_HOURS, n in NODES],
        sum(f_ac[h,ℓ] for ℓ in eachindex(AC_LINES) if AC_LINES[ℓ][1] == n) -
        sum(f_ac[h,ℓ] for ℓ in eachindex(AC_LINES) if AC_LINES[ℓ][2] == n) +
        g[h,n] + dc_net[h,n] == E_DEMAND[(h,n)] + h2_prod[h,n] / H2_EFFICIENCY[n]
    )

    ### Objective 
    @objective(m, Min,
        sum(HOUR_WEIGHT[h] * (
            sum(GEN_COST[n] * g[h,n] for n in NODES) +
            sum(H2_CONVERSION_COST[n] * h2_prod[h,n] for n in NODES)
        ) for h in REP_HOURS)
        +
        sum(H2_ANNUITY_FACTOR * H2_CAPEX_PER_MW[n] * h2_cap[n] for n in NODES)
    )

    optimize!(m)
    obj = objective_value(m)

    # Weighted H2 distribution
    total_onshore_h2  = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in ONSHORE_NODES)
    total_offshore_h2 = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in OFFSHORE_NODES)

    return m, g, f_ac, f_dc,
           h2_prod, h2_trans, h2_cap, h2_balance,
           node_balance, total_onshore_h2, total_offshore_h2, obj,
           h2_cap_bind, Dict(n => sum(-dual(h2_cap_bind[(h,n)]) for h in REP_HOURS) for n in NODES)
end

end # module
