module Data
import CSV, DataFrames
using CSV, DataFrames

export ZONES, NODES, zone_of, GENERATORS, GEN_COST, E_DEMAND,
       AC_LINES, DC_LINES, TCONNECT, H2_DEMAND, H2_CONVERSION_COST, H2_EFFICIENCY,
       H2_CONNECTIONS, H2_TRANS_LIMIT, ONSHORE_NODES, OFFSHORE_NODES, REP_HOURS, HOUR_WEIGHT,
       H2_CAPEX_PER_MW, H2_CAP_LIMIT, H2_ANNUITY_FACTOR

### Represenative Hours
const REP_HOURS = [:h1, :h2, :h3, :h4]  
const HOUR_WEIGHT = Dict(:h1 => 2232, :h2 => 2136, :h3 => 2712, :h4 => 1680) # weights for each Hour

### Load demand data from CSV file
using CSV, DataFrames

### Load demand and generator data from CSV files using path relative to this file
const basepath = joinpath(@__DIR__, "profiles")

const electricity_demand_df = CSV.read(joinpath(basepath, "electricity_demand.csv"), DataFrame)
const hydrogen_demand_df    = CSV.read(joinpath(basepath, "hydrogen_demand.csv"), DataFrame)
const generation_df         = CSV.read(joinpath(basepath, "generation.csv"), DataFrame)

### Zones and nodes
# Zones in the system
const ZONES = [:A, :B, :C, :D, :E, :F, :G]

# Nodes
const NODES = [:n1, :n2, :n3, :n4, :n5, :n6, :n7, :n8]

# Onshore and offshore nodes
const ONSHORE_NODES = [:n1, :n2, :n3, :n4, :n8] 
const OFFSHORE_NODES = [:n5, :n6, :n7] 

# Node-to-zone mapping
const zone_of = Dict(
    :n1 => :A,
    :n2 => :A,
    :n3 => :B,
    :n4 => :C,
    :n5 => :D,
    :n6 => :E,
    :n7 => :F,
    :n8 => :G   
)

### Generators for electricty generation
# Generators defined by (node, capacity in MW)
const GENERATORS = Dict(
    (Symbol(row.Hour), Symbol(col)) => generation_df[row_idx, col]
    for (row_idx, row) in enumerate(eachrow(generation_df)),
        col in names(generation_df)[2:end]  
)

# Specific generation cost per generator in €/MWh
const GEN_COST = Dict(
    :n1 => 50, 
    :n2 => 150,
    :n3 => 200,
    :n4 => 50,
    :n5 => 0,
    :n6 => 0,
    :n7 => 0,
    :n8 => 100
)

### Electricity systems

# Demand per node in MW
const E_DEMAND = Dict(
    (Symbol(row.Hour), Symbol(col)) => electricity_demand_df[row_idx, col]
    for (row_idx, row) in enumerate(eachrow(electricity_demand_df)),
        col in names(electricity_demand_df)[2:end]  
)

# AC lines
REACTANCE = 0.1  # reactance value in p.u.
THERMAL_LIMIT = 500 # thermal limit in MW
const AC_LINES = [
    (:n1, :n2, THERMAL_LIMIT, REACTANCE),   # A–A
    (:n1, :n4, THERMAL_LIMIT, REACTANCE),   # A–C
    (:n4, :n3, THERMAL_LIMIT, REACTANCE),  # C–B
    (:n3, :n2, THERMAL_LIMIT, REACTANCE)   # B–A
]

# DC lines
const DC_LINES = [
    (:n5, :n6, 1000),   # D-E
    (:n5, :n4, 1000),   # D–C
    (:n6, :n3, 1000),   # E–C
    (:n7, :n6, 1000),   # F–E
    (:n7, :n5, 1000),   # F–D
    (:n8, :n7, 1000)    # G-F
]

# Inter-zonal AC connections (for onshore AC grid)
const TCONNECT = [(:A, :C), (:C, :B), (:B, :A)]

### Hydrogen system

# Hydrogen CAPEX per MW per node
offshore_multipier = 1.25
const H2_CAPEX_PER_MW = Dict(
    :n1 => 0.8e6, 
    :n2 => 0.8e6, 
    :n3 => 0.8e6, 
    :n4 => 0.8e6,  
    :n5 => 0.8e6 * offshore_multipier,  # Offshore node
    :n6 => 0.8e6 * offshore_multipier,  # Offshore node
    :n7 => 0.8e6 * offshore_multipier,  # Offshore node
    :n8 => 0.8e6
)

# Hydrogen capacity limit per node in MW
const H2_CAP_LIMIT = Dict(
    :n1 => 1000, 
    :n2 => 1000, 
    :n3 => 1000, 
    :n4 => 1000,
    :n5 => 1000, 
    :n6 => 1000, 
    :n7 => 1000, 
    :n8 => 1000
)

# Annuity factor for hydrogen CAPEX
function annuity_factor(r, T)
    return r / (1 - (1 + r)^(-T))
end

DISCOUNT_RATE = 0.05
H2_LIFETIME_YEARS = 20
H2_ANNUITY_FACTOR = annuity_factor(DISCOUNT_RATE, H2_LIFETIME_YEARS)

# Hydrogen demand per node in MWh
const H2_DEMAND = Dict(
    (Symbol(row.Hour), Symbol(col)) => hydrogen_demand_df[row_idx, col]
    for (row_idx, row) in enumerate(eachrow(hydrogen_demand_df)),
        col in names(hydrogen_demand_df)[2:end]  # skip de "Hour"-kolom
)

# Specific hydrogen conversion cost per node in €/MWh H2
const H2_CONVERSION_COST = Dict(
    :n1 => 5,
    :n2 => 5,
    :n3 => 5,
    :n4 => 5,
    :n5 => 5,
    :n6 => 5,
    :n7 => 5,
    :n8 => 5
)

# Hydrogen conversion efficiency
const H2_EFFICIENCY = Dict(
    :n1 => 0.77,
    :n2 => 0.77,
    :n3 => 0.77,
    :n4 => 0.77,
    :n5 => 0.77,
    :n6 => 0.77,
    :n7 => 0.77,
    :n8 => 0.77
)

# Possible hydrogen transmission connections 
const H2_CONNECTIONS = [
    (:n1, :n2),  # A → A
    (:n1, :n4),  # B → A
    (:n2, :n3),  # A → C
    (:n3, :n4),  # B → C
    (:n3, :n6),  # B → E        
    (:n4, :n5),  # C → D
    (:n6, :n5),  # E → D
    (:n7, :n6),  # F → E
    (:n7, :n5),  # F → D
    (:n8, :n7)   # G → F
]

# Transmission limits for hydrogen connections in MWh
const H2_TRANS_LIMIT = Dict(
    (:n1, :n2) => 10000,  # A → A
    (:n1, :n4) => 10000,  # B → A
    (:n2, :n3) => 10000,  # A → C
    (:n3, :n4) => 10000,  # B → C
    (:n3, :n6) => 10000,  # B → C
    (:n4, :n5) => 10000,  # C → D
    (:n6, :n5) => 10000,  # E → D
    (:n7, :n6) => 10000,  # F → E
    (:n7, :n5) => 10000,  # F → D
    (:n8, :n7) => 10000   # G → F
)

end # module