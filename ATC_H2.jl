module ATC_H2

using JuMP, Gurobi
using CSV, DataFrames
import MathOptInterface as MOI

import ..Data: ZONES, NODES, zone_of, GENERATORS, GEN_COST, E_DEMAND,
               AC_LINES, DC_LINES, TCONNECT,
               H2_DEMAND, H2_CONVERSION_COST, H2_EFFICIENCY,
               H2_CONNECTIONS, H2_TRANS_LIMIT,
               ONSHORE_NODES, OFFSHORE_NODES, REP_HOURS, HOUR_WEIGHT,
               H2_CAP_LIMIT, H2_CAPEX_PER_MW, H2_ANNUITY_FACTOR

export atc_H2

# Onshore AC grid zones
const Onshore_AC_zone_set = unique(vcat([t[1] for t in TCONNECT], [t[2] for t in TCONNECT]))

# ATC with Exact Projection (Aravena et al., 2021), see paper for more information
# ATC bounds via Exact Projection
function _bounds_atc()
    m = Model(Gurobi.Optimizer)
    set_silent(m)

    # ATC bounds and box-vertex enumeration
    @variable(m, atc_plus[t in TCONNECT])
    @variable(m, atc_minus[t in TCONNECT])
    @constraint(m, [t in TCONNECT], -atc_minus[t] <= atc_plus[t])

    signs = collect(Iterators.product(fill([-1, 1], length(TCONNECT))...))
    verts = length(signs)

    # Electricity system
    @variable(m, 0 <= vbar[v = 1:verts, h in REP_HOURS, n in NODES] <= 1)     

    @variable(m, f_dc[v = 1:verts, h in REP_HOURS, ℓ = 1:length(DC_LINES)])
    @variable(m, net_pos[v = 1:verts, h in REP_HOURS, z in Onshore_AC_zone_set])               
    @variable(m, e[v = 1:verts, h in REP_HOURS, t in TCONNECT])               

    # Voltage angles and flows with nodal limits
    @variable(m, θ_hat[v = 1:verts, h in REP_HOURS, n in NODES])
 
    @variable(m, f_ac_hat[v = 1:verts, h in REP_HOURS, ℓ = 1:length(AC_LINES)])
    for (ℓ, (i, j, F, X)) in enumerate(AC_LINES)
        B = 1.0 / X
        @constraint(m, [v in 1:verts, h in REP_HOURS], f_ac_hat[v, h, ℓ] == B * (θ_hat[v, h, i] - θ_hat[v, h, j]))
        @constraint(m, [v in 1:verts, h in REP_HOURS], -F <= f_ac_hat[v, h, ℓ] <= F)
    end

    # Hydrogen system 
    @variable(m, 0 <= h2_prod_hat[v = 1:verts, h in REP_HOURS, n in NODES])
    @variable(m, -H2_TRANS_LIMIT[(i,j)] <= h2_trans_hat[v = 1:verts, h in REP_HOURS, (i,j) in H2_CONNECTIONS] <= H2_TRANS_LIMIT[(i,j)])
    @constraint(m, [v in 1:verts, h in REP_HOURS, n in NODES], h2_prod_hat[v,h,n] <= H2_CAP_LIMIT[n])

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
                    if !(v2 in seen)
                        push!(seen, v2); push!(stack, v2)
                    end
                end
            end
            @constraint(m, [v in 1:verts, h in REP_HOURS], θ_hat[v, h, comp[1]] == 0)
        end
    end

    # Constraints for each vertex of the polytope
    for (vtx, s) in enumerate(signs), h in REP_HOURS
        for (k, t) in enumerate(TCONNECT)
            @constraint(m, e[vtx, h, t] == (s[k] == 1 ? atc_plus[t] : -atc_minus[t]))
        end

        for (ℓ, (_, _, F)) in enumerate(DC_LINES)
            @constraint(m, -F <= f_dc[vtx, h, ℓ] <= F)
        end

        for z in Onshore_AC_zone_set
            @constraint(m, net_pos[vtx, h, z] ==
                sum(e[vtx, h, t] for t in TCONNECT if t[1] == z; init=0.0) -
                sum(e[vtx, h, t] for t in TCONNECT if t[2] == z; init=0.0))
        end
        @constraint(m, sum(net_pos[vtx, h, z] for z in Onshore_AC_zone_set) == 0)
        
        # H2 balance
        @constraint(m, [n in NODES],
            h2_prod_hat[vtx,h,n] +
            sum(h2_trans_hat[vtx,h,(i,n)] for (i,j) in H2_CONNECTIONS if j == n; init=0.0) -
            sum(h2_trans_hat[vtx,h,(n,j)] for (i,j) in H2_CONNECTIONS if i == n; init=0.0)
            == H2_DEMAND[(h,n)]
        )

        # Nodal balance
        for n in NODES
            ac_out = sum(f_ac_hat[vtx, h, ℓ] for (ℓ,(i,_,_,_)) in enumerate(AC_LINES) if i == n; init=0.0)
            ac_in  = sum(f_ac_hat[vtx, h, ℓ] for (ℓ,(_,j,_,_)) in enumerate(AC_LINES) if j == n; init=0.0)
            dc_inj = (+ sum(f_dc[vtx, h, ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if j == n; init=0.0)
                      - sum(f_dc[vtx, h, ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if i == n; init=0.0))
            inj = GENERATORS[(h, n)] * vbar[vtx, h, n] -
                  (E_DEMAND[(h, n)] + h2_prod_hat[vtx,h,n] / H2_EFFICIENCY[n]) +
                  dc_inj
            @constraint(m, ac_out - ac_in == inj)
        end
        # Zonal balance with auxillary variable DC
        for z in Onshore_AC_zone_set
            dc_in_z  = sum(f_dc[vtx, h, ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if zone_of[j] == z; init=0.0)
            dc_out_z = sum(f_dc[vtx, h, ℓ] for (ℓ,(i,j,_)) in enumerate(DC_LINES) if zone_of[i] == z; init=0.0)
            lhs =  sum(GENERATORS[(h,n)] * vbar[vtx, h, n] for n in NODES if zone_of[n] == z; init=0.0)
            lhs -= (sum(E_DEMAND[(h,n)] for n in NODES if zone_of[n] == z; init=0.0) +
                    sum(h2_prod_hat[vtx,h,n] / H2_EFFICIENCY[n] for n in NODES if zone_of[n] == z; init=0.0))
            lhs += dc_in_z - dc_out_z
            @constraint(m, lhs == net_pos[vtx, h, z])
        end

    end
    # Maximize ATC size (Exact Projection objective)
    @objective(m, Max, sum(atc_plus[t] + atc_minus[t] for t in TCONNECT))
    optimize!(m)

    ATCplus  = value.(atc_plus)
    ATCminus = value.(atc_minus)
    atc_vector = [(abs(ATCplus[t]), -abs(ATCminus[t])) for t in TCONNECT]
    return atc_vector
end

function _atc_market(atc::Vector{Tuple{Float64,Float64}})
    atc_plus  = Dict(t => atc[i][1] for (i, t) in enumerate(TCONNECT))
    atc_minus = Dict(t => atc[i][2] for (i, t) in enumerate(TCONNECT))

    m = Model(Gurobi.Optimizer)
    set_silent(m)

    @variable(m, 0 <= v[h in REP_HOURS, n in NODES] <= 1)

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

    # DC 
    @variable(m, f_dc[h in REP_HOURS, ℓ in 1:length(DC_LINES)])
    for (ℓ, (_, _, F)) in enumerate(DC_LINES)
        @constraint(m, [h in REP_HOURS], -F <= f_dc[h,ℓ] <= F)
    end

    # Zonal exchanges & net positions
    @variable(m, e[h in REP_HOURS, t in TCONNECT])         
    @variable(m, net_pos[h in REP_HOURS, z in ZONES])      

    # For onshore AC grid
    @constraint(m, [h in REP_HOURS, z in Onshore_AC_zone_set],
        net_pos[h,z] ==
            sum(e[h,t] for t in TCONNECT if t[1] == z; init=0.0) -
            sum(e[h,t] for t in TCONNECT if t[2] == z; init=0.0))
    @constraint(m, [h in REP_HOURS, z in setdiff(ZONES, Onshore_AC_zone_set)], net_pos[h,z] == 0)

    # Limits on exchanges
    @constraint(m, [h in REP_HOURS, t in TCONNECT], atc_minus[t] <= e[h,t] <= atc_plus[t])

    # Net position must sum to zero across all zones
    @constraint(m, [h in REP_HOURS], sum(net_pos[h,z] for z in ZONES) == 0)

    # Zonal balance
    zone_balance = Dict{Tuple{Symbol, Symbol}, ConstraintRef}()
    for h in REP_HOURS, z in ZONES
        dc_in  = sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][2]] == z; init=0.0)
        dc_out = sum(f_dc[h,ℓ] for ℓ in eachindex(DC_LINES) if zone_of[DC_LINES[ℓ][1]] == z; init=0.0)
        zone_balance[(h,z)] = @constraint(m,
            sum(GENERATORS[(h,n)] * v[h,n] for n in NODES if zone_of[n] == z; init=0.0)
            - net_pos[h,z] + dc_in - dc_out ==
            sum(E_DEMAND[(h,n)] for n in NODES if zone_of[n] == z; init=0.0) +
            sum(h2_prod[h,n] / H2_EFFICIENCY[n] for n in NODES if zone_of[n] == z; init=0.0)
        )
    end

    ### Objective
    @objective(m, Min,
        sum(HOUR_WEIGHT[h] * (
            sum(GEN_COST[n] * GENERATORS[(h,n)] * v[h,n] for n in NODES; init=0.0) +
            sum(H2_CONVERSION_COST[n] * h2_prod[h,n] for n in NODES; init=0.0)
        ) for h in REP_HOURS; init=0.0)
        + sum(H2_ANNUITY_FACTOR * H2_CAPEX_PER_MW[n] * h2_cap[n] for n in NODES; init=0.0)
    )

    optimize!(m)
    obj = objective_value(m)

    # Weighted H2 distribution
    total_onshore_h2  = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in ONSHORE_NODES; init=0.0)
    total_offshore_h2 = sum(HOUR_WEIGHT[h] * value(h2_prod[h,n]) for h in REP_HOURS for n in OFFSHORE_NODES; init=0.0)

    return m, v, f_dc, net_pos, e,
           h2_prod, h2_trans, h2_cap, h2_balance, zone_balance,
           atc_plus, atc_minus, total_onshore_h2, total_offshore_h2, obj,
           h2_cap_bind, Dict(n => sum(-dual(h2_cap_bind[(h,n)]) for h in REP_HOURS; init=0.0) for n in NODES)
end

function atc_H2()
    atc_bounds = _bounds_atc()
    atc_run = _atc_market(atc_bounds)
    return atc_run, atc_bounds
end

end # module
