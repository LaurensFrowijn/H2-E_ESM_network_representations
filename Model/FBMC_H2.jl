module FBMC_H2

using JuMP, Gurobi
using CSV, DataFrames
import MathOptInterface as MOI
import ..Data: ZONES, NODES, zone_of, GENERATORS, GEN_COST, E_DEMAND,
               AC_LINES, DC_LINES, TCONNECT, H2_DEMAND, H2_CONVERSION_COST, H2_EFFICIENCY, H2_CONNECTIONS,
               H2_TRANS_LIMIT, ONSHORE_NODES, OFFSHORE_NODES, REP_HOURS, HOUR_WEIGHT,
               H2_CAPEX_PER_MW, H2_CAP_LIMIT, H2_ANNUITY_FACTOR

export fbmc_H2

# FBMC with Exact Projection (Aravena et al., 2021), see paper for more information
function fbmc_H2()
    m = Model(Gurobi.Optimizer)
    set_silent(m)

    # Electricity system
    @variable(m, 0 <= v[h in REP_HOURS, n in NODES] <= 1)

    # Hydrogen system 
    @variable(m, 0 <= h2_cap[n in NODES] <= H2_CAP_LIMIT[n])
    @variable(m, 0 <= h2_prod[h in REP_HOURS, n in NODES])

    @variable(m, -H2_TRANS_LIMIT[(i,j)] <= h2_trans[h in REP_HOURS, (i,j) in H2_CONNECTIONS] <= H2_TRANS_LIMIT[(i,j)])

    h2_cap_bind = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for h in REP_HOURS, n in NODES
        h2_cap_bind[(h,n)] = @constraint(m, h2_prod[h,n] <= h2_cap[n])
    end

    # Hydrogen nodal balances
    h2_balance = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for h in REP_HOURS, n in NODES
        h2_in  = sum(h2_trans[h,(i,n)] for (i,j) in H2_CONNECTIONS if j == n; init=0.0)
        h2_out = sum(h2_trans[h,(n,j)] for (i,j) in H2_CONNECTIONS if i == n; init=0.0)
        h2_balance[(h,n)] = @constraint(m, h2_prod[h,n] + h2_in - h2_out == H2_DEMAND[(h,n)])
    end

    # HVDC flows 
    @variable(m, f_dc[h in REP_HOURS, ℓ in 1:length(DC_LINES)])
    for (ℓ, (_, _, F)) in enumerate(DC_LINES)
        @constraint(m, [h in REP_HOURS], -F <= f_dc[h,ℓ] <= F)
    end

    # Zonal net positions 
    @variable(m, net_pos[h in REP_HOURS, z in ZONES])
    zone_balance = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for h in REP_HOURS, z in ZONES
        dc_in  = sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][2]] == z; init=0.0)
        dc_out = sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][1]] == z; init=0.0)
        zone_balance[(h,z)] = @constraint(m,
            sum(GENERATORS[(h,n)] * v[h,n] for n in NODES if zone_of[n] == z; init=0.0)
            - net_pos[h,z] + dc_in - dc_out
            ==
            sum(E_DEMAND[(h,n)] for n in NODES if zone_of[n] == z; init=0.0)
          + sum(h2_prod[h,n] / H2_EFFICIENCY[n] for n in NODES if zone_of[n] == z; init=0.0)
        )
    end
    @constraint(m, [h in REP_HOURS], sum(net_pos[h,z] for z in ZONES) == 0)

    ### Exact projection (Aravena et al. 2021) ###

    # Auxillary variable nodal generation fraction
    @variable(m, 0 <= vbar[h in REP_HOURS, n in NODES] <= 1)

    @variable(m, f_dc_hat[h in REP_HOURS, ℓ in 1:length(DC_LINES)])
    for (ℓ, (_, _, F)) in enumerate(DC_LINES)
        @constraint(m, [h in REP_HOURS], -F <= f_dc_hat[h,ℓ] <= F)
    end

    # Voltage angles 
    @variable(m, θ_hat[h in REP_HOURS, n in NODES])

    # Fix one angle reference
    adj = Dict(n => Symbol[] for n in NODES)
    for (i,j,_,_) in AC_LINES
        push!(adj[i], j); push!(adj[j], i)
    end
    seen = Set{Symbol}()
    for n in NODES
        if !(n in seen)
            comp = Symbol[]
            stack = [n]; push!(seen, n)
            while !isempty(stack)
                u = pop!(stack); push!(comp, u)
                for v2 in adj[u]
                    if !(v2 in seen); push!(seen, v2); push!(stack, v2); end
                end
            end
            @constraint(m, [h in REP_HOURS], θ_hat[h, comp[1]] == 0)
        end
    end

    # Auxillary variable AC flows
    @variable(m, f_ac_hat[h in REP_HOURS, ℓ in 1:length(AC_LINES)])

    for (ℓ, line) in enumerate(AC_LINES)
        i, j, F, X = line
        B = 1.0 / X
        @constraint(m, [h in REP_HOURS], f_ac_hat[h,ℓ] == B * (θ_hat[h, i] - θ_hat[h, j]))
        @constraint(m, [h in REP_HOURS], -F <= f_ac_hat[h,ℓ] <= F)
    end

    # Nodal balance
    @variable(m, 0 <= h2_prod_hat[h in REP_HOURS, n in NODES])
    @constraint(m, [h in REP_HOURS, n in NODES], h2_prod_hat[h,n] <= h2_cap[n])
    @variable(m, -H2_TRANS_LIMIT[(i,j)] <= h2_trans_hat[h in REP_HOURS, (i,j) in H2_CONNECTIONS] <= H2_TRANS_LIMIT[(i,j)])
    h2_balance_hat = Dict{Tuple{Symbol,Symbol},ConstraintRef}()
    for h in REP_HOURS, n in NODES
        h2_in_hat  = sum(h2_trans_hat[h,(i,n)] for (i,j) in H2_CONNECTIONS if j == n; init=0.0)
        h2_out_hat = sum(h2_trans_hat[h,(n,j)] for (i,j) in H2_CONNECTIONS if i == n; init=0.0)
        h2_balance_hat[(h,n)] = @constraint(m, h2_prod_hat[h,n] + h2_in_hat - h2_out_hat == H2_DEMAND[(h,n)])
    end
    for h in REP_HOURS, n in NODES
        ac_out = sum(f_ac_hat[h,ℓ] for (ℓ,(i,_,_,_)) in enumerate(AC_LINES) if i == n; init=0.0) 
        ac_in  = sum(f_ac_hat[h,ℓ] for (ℓ,(_,j,_,_)) in enumerate(AC_LINES) if j == n; init=0.0)
        dc_inj = (  sum(f_dc_hat[h,ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if j == n; init=0.0)
                  - sum(f_dc_hat[h,ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if i == n; init=0.0) )
        inj = GENERATORS[(h,n)] * vbar[h,n] - (E_DEMAND[(h,n)] + h2_prod_hat[h,n] / H2_EFFICIENCY[n]) + dc_inj
        @constraint(m, ac_out - ac_in == inj)
    end

    # Zonal balance with auxillary variable DC 
    for h in REP_HOURS, z in ZONES
        dc_in_z  = sum(f_dc_hat[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][2]] == z; init=0.0)
        dc_out_z = sum(f_dc_hat[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][1]] == z; init=0.0)
        lhs =  sum(GENERATORS[(h,n)] * vbar[h,n] for n in NODES if zone_of[n] == z; init=0.0)
        lhs -= (sum(E_DEMAND[(h,n)] for n in NODES if zone_of[n] == z; init=0.0) +
                sum(h2_prod_hat[h,n] / H2_EFFICIENCY[n] for n in NODES if zone_of[n] == z; init=0.0))
        lhs += dc_in_z - dc_out_z
        @constraint(m, lhs == net_pos[h,z])
    end

    ### Objective
    @objective(m, Min,
        sum(HOUR_WEIGHT[h] * (
            sum(GEN_COST[n] * GENERATORS[(h,n)] * v[h,n] for n in NODES) +
            sum(H2_CONVERSION_COST[n] * h2_prod[h,n] for n in NODES)
        ) for h in REP_HOURS)
        + sum(H2_ANNUITY_FACTOR * H2_CAPEX_PER_MW[n] * h2_cap[n] for n in NODES)
    )

    optimize!(m)
    obj = objective_value(m)

    # Weighted H2 distribution
    total_onshore_h2  = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in ONSHORE_NODES)
    total_offshore_h2 = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in OFFSHORE_NODES)

    return m, v, f_ac_hat, f_dc, net_pos,
           h2_prod, h2_trans, h2_cap, h2_balance, zone_balance,
           total_onshore_h2, total_offshore_h2, obj,
           h2_cap_bind, Dict(n => sum(-dual(h2_cap_bind[(h,n)]) for h in REP_HOURS) for n in NODES)
end

end # module
