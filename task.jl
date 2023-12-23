## Step 0: Activate environment
using Pkg
Pkg.instantiate()
using PowerModels, Ipopt, JuMP,Cbc, Juniper#, Gurobi
using Plots


# Define solver
optimizer = Juniper.Optimizer
nl_solver= optimizer_with_attributes(Ipopt.Optimizer,  "print_level" => 0)
#mip_solver = optimizer_with_attributes(Gurobi.Optimizer,"logLevel" => 0)
#juniper = optimizer_with_attributes(optimizer, "nl_solver"=>nl_solver,  "mip_solver"=>mip_solver)
juniper = optimizer_with_attributes(Ipopt.Optimizer)
m = Model(juniper)


# ipopt = optimizer_with_attributes(Juniper.Optimizer)

##### Step 1: Import the grid data and initialize the JuMP model
# case_file = "PowerModelsACDC.jl-master/test/data/tnep/case9_test.m"
case_file = "PowerModelsACDC.jl-master/test/data/tnep/pglib_opf_case5_pjm.m"
# case_file = "PowerModelsACDC.jl-master/test/data/tnep/PSCC/case9.m"
#case_file = "PowerModelsACDC.jl-master\\test\\data\\tnep\\case14_test.m"


# For convenience, use the parser of Powermodels to convert the MATPOWER format file to a Julia dictionary
data = PowerModels.parse_file(case_file)


# Initialize the JuMP model (an empty JuMP model) with defined solver
#m = Model(optimizer)
#m = Model(optimizer_with_attributes(optimizer))
##### Step 2: create the JuMP model & pass data to model
include("info.jl") # Define functions define_sets! and process_parameters!
define_sets!(m, data) # Pass the sets to the JuMP model
process_parameters!(m, data) # Pass the parameters to the JuMP model




##### Step 3: Build the model
include("builing.jl") # Define build_ac_opf! function
build_ac_opf!(m) # Pass the model to the build_ac_opf! function

##### Step 4: Solve the model
optimize!(m) # Solve the model
println(objective_value(m)) # Print the objective value of the model

##### Compare the two objective functions
result_pm = PowerModels.solve_opf(case_file, ACPPowerModel, juniper) # Solve using PowerModels and retrieve the solutions
print(Dict("objective"=>objective_value(m),"objective_pm"=>result_pm["objective"])) # Compare the objective values


pg = value.(m.ext[:variables][:pg])

#####