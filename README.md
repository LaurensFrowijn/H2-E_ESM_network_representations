# Hydrogen-Electricity energy system model with three network representations: ATC, DC PF, and FBMC

This repository contains the Julia models accompanying the paper “On- or offshore electrolytic hydrogen production? The role of electricity network representation” (Frowijn et al., 2025).
The models implement three electricity–hydrogen market formulations; Available Transfer Capacity (ATC), DC Power Flow (DCPF), and Flow-Based Market Coupling (FBMC), and compute minimized cost for electricity dispatch, hydrogen production/transport, and investment decisions on a stylized multi-node system. For scientific context, assumptions, and interpretation, please see the paper; the repository provides runnable code.
When citing, cite both the paper and this model repository.

Reproducibility notes:
The solver is set to Gurobi in each model file. To use another solver, replace Model(Gurobi.Optimizer) accordingly and ensure feature parity. The default run prints a console summary; extend model/main.jl if you need CSV exports or additional KPIs.

Requirements:
- Julia 1.9 or newer
- Gurobi installed with a valid license
- Julia packages: JuMP, Gurobi, CSV, DataFrames (MathOptInterface is imported internally by the model files)

Installation (packages):
Use the Julia package manager to add the required packages: JuMP, Gurobi, CSV, DataFrames. Make sure Gurobi.jl is linked to your local Gurobi installation.

Model:
Input data (profiles):
Place electricity_demand.csv, generation.csv, and hydrogen_demand.csv in the profiles/ directory. Each file should contain hourly values per node in MWh. The model/Data.jl script loads these profiles and assembles all other inputs (nodes/zones, AC/DC lines, costs, efficiencies, capacity limits, and representative hours).

Running:
Run model/Main.jl in Julia. This executes the ATC, DCPF, and FBMC models and prints, for each model, the objective value (MEUR) and the weighted totals of onshore and offshore hydrogen production (GWh). Data.jl is located in model/ and is used by all three models.
To run different scenarios or with a different dataset, change the input data in Data.jl

Repository structure:
Repository structure

model/
- Data.jl — Model data (nodes/zones, lines, costs, limits, representative hours, etc.)
- ATC.jl — ATC model
- DCPF.jl — DC power flow model
- FBMC.jl — Flow-Based model
- main.jl — Runs ATC, DCPF, FBMC

-profiles/
-- electricity_demand.csv — MWh per hour per node
-- generation.csv — MWh per hour per node (capacities)
-- hydrogen_demand.csv — MWh (LHV) per hour per node

Contact / issues:
Open a GitHub issue for questions or bug reports. For conceptual questions or interpretation, refer to and cite the paper alongside this model repository.
