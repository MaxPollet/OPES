## Step 0: Activate environment
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate()
using PowerModels, Ipopt, JuMP,Cbc, Juniper, Gurobi, PolyhedralRelaxations
using Plots


# Define solver
optimizer = Juniper.Optimizer
nl_solver= optimizer_with_attributes(Ipopt.Optimizer,  "print_level" => 0)
mip_solver = optimizer_with_attributes(Cbc.Optimizer)
juniper = optimizer_with_attributes(Gurobi.Optimizer)
m = Model(juniper)



# ipopt = optimizer_with_attributes(Juniper.Optimizer)

##### Step 1: Import the grid data and initialize the JuMP model
case_file = "data/case9.m"
# case_file = "data/case14.m"


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
# include("LPAC_no_MI.jl") # Working LPAC with no integer
build_ac_opf!(m) # Pass the model to the build_ac_opf! function

##### Step 4: Solve the model
optimize!(m) # Solve the model
println(objective_value(m)) # Print the objective value of the model

##### Compare the two objective functions
# result_pm = PowerModels.solve_opf(case_file, ACPPowerModel, juniper) # Solve using PowerModels and retrieve the solutions
# print(Dict("objectilnve"=>objective_value(m),"objective_pm"=>result_pm["objective"])) # Compare the objective values


# pg = value.(m.ext[:variables][:pg])

#####

####### see what the power is on the ac branch
# println("Power through ac branch")
# for (n,j,i) in m.ext[:sets][:B_ac]
#     print((n,j,i))
#     print("=>")
#     println(value.(m.ext[:variables][:pb][(n,j,i)]))
    
# end



####### see which dc branch is switched on or off

# println("see which dc branch is switched on or off")
# for (n,j,i) in m.ext[:sets][:B_dc_to]
#     c = parse(Float64, n)
#     if (c >= 10)
#         print(n)
#         print("=>")
#         println(value.(m.ext[:variables][:ed][n]))
#     else
#         print(n)
#         print(" =>")
#         println(value.(m.ext[:variables][:ed][n]))
#     end
# end


####### see which dc branch is switched on 

# println("see which dc branch is switched on ")
# for (n,j,i) in m.ext[:sets][:B_dc_to]
#     c = parse(Float64, n)
#     if (c >= 10)
    
#         if value.(m.ext[:variables][:ed][n]) > 0.1
#             print((n,j,i) )
#             print("=>")
#             println(value.(m.ext[:variables][:ed][n,]))
#             # println(value.(m.ext[:variables][:ed][n]))
#         end
#     else
        
#         if value.(m.ext[:variables][:ed][n]) > 0.1
#             print((n,j,i) )
#             print(" =>")
#             println(value.(m.ext[:variables][:ed][n,]))
#             # println(value.(m.ext[:variables][:ed][n]))
#         end
#     end
# end


# println("see the dc branch power if non zero ")
for (n,j,i) in m.ext[:sets][:B_dc_to]
    c = parse(Float64, n)
    if (c >= 10)
    
        if value.(m.ext[:variables][:ed][n]) > 0.1
            print((n,j,i) )
            print("=>")
            println(value.(m.ext[:variables][:pb_dc][(n,j,i)]))
            # println(value.(m.ext[:variables][:ed][n]))
        end
    else
        
        if value.(m.ext[:variables][:ed][n]) > 0.1
            print((n,j,i) )
            print(" =>")
            println(value.(m.ext[:variables][:pb_dc][(n,j,i)]))
            # println(value.(m.ext[:variables][:ed][n]))
        end
    end
end


# println("see the power from converter into dc net if non zero ")
# for n in m.ext[:sets][:N_tf]
#     c = parse(Float64, n)
#     if (c >= 10)
    
#         if value.(m.ext[:variables][:ec][n]) > 0.1
#             print(n )
#             print("=>")
#             println(value.(m.ext[:variables][:pb_dc_node][n]))
#             # println(value.(m.ext[:variables][:ed][n]))
#         end
#     else
        
#         if value.(m.ext[:variables][:ec][n]) > 0.1
#             print(n)
#             print(" =>")
#             println(value.(m.ext[:variables][:pb_dc_node][n]))
#             # println(value.(m.ext[:variables][:ed][n]))
#         end
#     end
# end





####### see which dc converter is switched on or off
# println("see which dc converter is switched on or off ")
# for (n,j,i) in m.ext[:sets][:B_tf_dc_ac_fr]
#     c = parse(Float64, n)
#     if (c >= 10)
#         print(n)
#         print("=>")
#         println(value.(m.ext[:variables][:ec][n]))
#     else
#         print(n)
#         print(" =>")
#         println(value.(m.ext[:variables][:ec][n]))
#     end
# end




####### see which dc converter is switched on 
println("see which dc converter is switched on ")
for (n,j,i) in m.ext[:sets][:B_tf_dc_ac_fr]
    c = parse(Float64, n)
    if (c >= 10)
    
        if value.(m.ext[:variables][:ec][n]) > 0.1
            print(n)
            print("=>")
            println(value.(m.ext[:variables][:ec][n]))
        end
    else
        
        if value.(m.ext[:variables][:ec][n]) > 0.1
            print(n)
            print(" =>")
            println(value.(m.ext[:variables][:ec][n]))
        end
    end
end



# println("see what the power through the transfor to the dc side is: ")
# for n in m.ext[:sets][:N_tf]
#         print(n)
#         print("=>")
#         println(value.(m.ext[:variables][:pb_tf_ac_dc][(n,n,n)]))
# end



# println("Total DC losses ")
# sumy = 0
# for n in m.ext[:sets][:N]
#         sumy = sumy + value(m.ext[:variables][:pb_tf_ac_dc][(n,n,n)])

# end

# print("=>")
# println(sumy)


# println("Total AC branch losses ")
# sumy = 0
# for (n,i,j) in m.ext[:sets][:B_ac]
#         sumy = sumy + value(m.ext[:variables][:pb][(n,i,j)])

# end

# print("=>")
# println(sumy)


