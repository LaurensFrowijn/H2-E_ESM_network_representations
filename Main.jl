include("Data.jl")
include("FBMC_H2.jl")
include("ATC_H2.jl")
include("DCPF_H2.jl")



using JuMP, Gurobi
using .Data, .FBMC_H2, .ATC_H2, .DCPF_H2

function main()
    # === ATC ===
    atc_run, atc_vec = atc_H2()

    total_onshore_h2_atc   = atc_run[13]
    total_offshore_h2_atc  = atc_run[14]
    obj_atc                = atc_run[15]

    # === DCPF ===
    dcpf_run = dcpf_H2()
    total_onshore_h2_dcpf  = dcpf_run[10]
    total_offshore_h2_dcpf = dcpf_run[11]
    obj_dcpf               = dcpf_run[12]

    # === FBMC ===
    fb_run = fbmc_H2()
    total_onshore_h2_fb    = fb_run[11]
    total_offshore_h2_fb   = fb_run[12]
    obj_fb                 = fb_run[13]

    println("\n=== Objective & Hydrogen Production Summary ===")
    println("\nATC:")
    println("  Objective = $(round(obj_atc / 1e6, digits=3)) MEUR")
    println("  Total Hydrogen production Onshore = $(round(total_onshore_h2_atc / 1000, digits=0)) GWh")
    println("  Total Hydrogen production Offshore = $(round(total_offshore_h2_atc / 1000, digits=0)) GWh")

    println("\nDCPF:")
    println("  Objective = $(round(obj_dcpf / 1e6, digits=3)) MEUR")
    println("  Total Hydrogen production Onshore = $(round(total_onshore_h2_dcpf / 1000, digits=0)) GWh")
    println("  Total Hydrogen production Offshore = $(round(total_offshore_h2_dcpf / 1000, digits=0)) GWh")

    println("\nFlow-Based:")
    println("  Objective = $(round(obj_fb / 1e6, digits=3)) MEUR")
    println("  Total Hydrogen production Onshore = $(round(total_onshore_h2_fb / 1000, digits=0)) GWh")
    println("  Total Hydrogen production Offshore = $(round(total_offshore_h2_fb / 1000, digits=0)) GWh")

end

main()
