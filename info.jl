

# Create a function to define sets and pass it to the function
function define_sets!(m::Model, data::Dict)
    # Create a dictionary for the sets
    m.ext[:sets] = Dict()
    # Sets of elements
    # Set of AC nodes
    N = m.ext[:sets][:N] = [bus_id for (bus_id,bus) in data["bus"]]
    # Set of slack nodes
    N_sl = m.ext[:sets][:N_sl] = [bus_id for (bus_id,bus) in data["bus"] if bus["bus_type"]==3]
    # Set of AC branches
    B = m.ext[:sets][:B] = [br_id for (br_id,branch) in data["branch"]]
    # Set of shunt elements
    S = m.ext[:sets][:S] = [s_id for (s_id,shunt) in data["shunt"]]
    # Set of generators
    G = m.ext[:sets][:G] = [gen_id for (gen_id,gen) in data["gen"]]
    # Set of loads
    L = m.ext[:sets][:L] = [load_id for (load_id,load) in data["load"]]



    

    # Set of AC topology from side (i->j) and to side (j->i)
    B_ac_fr = m.ext[:sets][:B_ac_fr] = [(br_id, string(br["f_bus"]), string(br["t_bus"])) for (br_id,br) in data["branch"]] 
    B_ac_to = m.ext[:sets][:B_ac_to] = [(br_id, string(br["t_bus"]), string(br["f_bus"])) for (br_id,br) in data["branch"]]
    # Build a union set of both sets above
    B_ac = m.ext[:sets][:B_ac] = [B_ac_fr; B_ac_to]
    # Branch connectivity to buses, i.e. which branches are connected to a certain node, used in nodal power balance equations
    B_arcs = m.ext[:sets][:B_arcs] = Dict((i, []) for i in N) # create a
    for (l,i,j) in B_ac
        push!(B_arcs[i], (l,i,j))
    end


    # Set of DC nodes
    N_dc = m.ext[:sets][:N_dc] = [bus_id for (bus_id,bus) in data["busdc_ne"]]
    # Set of DC branches
    B_dc = m.ext[:sets][:B_dc] = [br_id for (br_id,branch) in data["branchdc_ne"]]

    # Set of DC topology from side (i->j) and to side (j->i)
    B_dc_fr = m.ext[:sets][:B_dc_fr] = [(br_id, string(br["fbusdc"]), string(br["tbusdc"])) for (br_id,br) in data["branchdc_ne"]] 
    B_dc_to = m.ext[:sets][:B_dc_to] = [(br_id, string(br["tbusdc"]), string(br["fbusdc"])) for (br_id,br) in data["branchdc_ne"]]
    # Build a union set of both sets above
    B_dc_to_fr = m.ext[:sets][:B_dc_to_fr] = [B_dc_fr; B_dc_to]
    # Branch connectivity to buses, i.e. which branches are connected to a certain node, used in nodal power balance equations
    B_dc_arcs = m.ext[:sets][:B_dc_arcs] = Dict((i, []) for i in N_dc) # create a
    for (l,i,j) in B_dc_to_fr
        push!(B_dc_arcs[i], (l,i,j))
    end


    B_dc_aux_fr = m.ext[:sets][:B_dc_aux_fr] = [(br_id, string(br["fbusdc"])) for (br_id,br) in data["branchdc_ne"]] 
    B_dc_aux_to = m.ext[:sets][:B_dc_aux_to] = [(br_id, string(br["tbusdc"])) for (br_id,br) in data["branchdc_ne"]]
    # Build a union set of both sets above
    B_dc_aux_to_fr = m.ext[:sets][:B_dc_aux_to_fr] = [B_dc_aux_fr; B_dc_aux_to]



    ###Transformer
    N_tf = m.ext[:sets][:N_tf] = [bus_id for (bus_id,bus) in data["convdc_ne"]]
    B_tf_dc_ac_fr = m.ext[:sets][:B_tf_dc_ac_fr] = [(tf_id, string(tf["busdc_i"]), string(tf["busac_i"])) for (tf_id,tf) in data["convdc_ne"]]
    B_tf_dc_ac_to = m.ext[:sets][:B_tf_dc_ac_to] = [(tf_id, string(tf["busac_i"]), string(tf["busdc_i"])) for (tf_id,tf) in data["convdc_ne"]]
    B_tf_dc_ac = m.ext[:sets][:B_tf_dc_ac] = [(tf_id, string(tf["busdc_i"]), string(tf["busac_i"])) for (tf_id,tf) in data["convdc_ne"]]
    B_tf_ac_dc = m.ext[:sets][:B_tf_ac_dc] = [(tf_id, string(tf["busac_i"]), string(tf["busdc_i"])) for (tf_id,tf) in data["convdc_ne"]]


    B_tf_dc_ac_to_fr = m.ext[:sets][:B_tf_dc_ac_to_fr] = [B_tf_dc_ac_fr; B_tf_dc_ac_to]
    
    N_dc_ac_connected = m.ext[:sets][:N_dc_ac_connected] = [string(tf["busac_i"]) for (tf_id,tf) in data["convdc_ne"]]
    N_no_dc_connected = m.ext[:sets][:N_no_dc_connected] =Vector()
    for j in N
            
        if (j in N_dc_ac_connected) == false
            push!(N_no_dc_connected,j)
        end
    end

    N_all_dc_connected = m.ext[:sets][:N_all_dc_connected] = Dict((i, []) for i in N_dc_ac_connected)

    for (b,j,i) in B_tf_dc_ac_to
        if (j in N_dc_ac_connected)
            push!(N_all_dc_connected[j],(i,j))
        end
    end

    ####Shunt component




    ### phase reactor 




    # Shunt connectivity to buses, i.e. which branches are connected to a certain node, used in nodal power balance equations
    S_ac = m.ext[:sets][:S_ac] = Dict((i, []) for i in N)
    for (s, shunt) in data["shunt"]
        push!(S_ac[string(shunt["shunt_bus"])], s)
    end
    # Generator connectivity, i.e. which generators are connected to a certain node, used in nodal power balance equations
    G_ac = m.ext[:sets][:G_ac] = Dict((i, []) for i in N)
    for (g, gen) in data["gen"]
        push!(G_ac[string(gen["gen_bus"])], g)
    end
    # Load connectivity, i.e. which loads are connected to a certain node, used in nodal power balance equations
    L_ac = m.ext[:sets][:L_ac] = Dict((i, []) for i in N)
    for (l, load) in data["load"]
        push!(L_ac[string(load["load_bus"])], l)
    end



    N_dc_index = m.ext[:sets][:N_dc_index] = Dict((i, []) for i in N)
    for (bus_id, bus) in data["busdc_ne"]
        push!(N_dc_index[string(bus_id)], bus_id)
    end


    return
end

# Create a function to pass the grid data to the JuMP model
function process_parameters!(m::Model, data::Dict)
    # Extract sets
    N = m.ext[:sets][:N]
    B = m.ext[:sets][:B]
    G = m.ext[:sets][:G]
    L = m.ext[:sets][:L]
    S = m.ext[:sets][:S]
    N_dc = m.ext[:sets][:N_dc]
    B_dc = m.ext[:sets][:B_dc]
    N_tf = m.ext[:sets][:N_tf]


    # Create parameter dictionary
    m.ext[:parameters] = Dict()

    baseMVA = m.ext[:parameters][:baseMVA] = data["baseMVA"] # get the base MVA
    # baseMVA = 1 # get the base MVA

    Pacmax = data["convdc_ne"]["1"]["Pacmax"]
    Qacmax = data["convdc_ne"]["1"]["Qacmax"]

    Sa = sqrt(Pacmax^2 + Qacmax^2) 

    # Zbase = (data["bus"]["1"]["base_kv"])^2/100
    Zbase = 1

    # Bus parameters
    vmmin = m.ext[:parameters][:vmmin] = Dict(i => data["bus"][i]["vmin"] for i in N) # minimum voltage magnitude
    vmmax = m.ext[:parameters][:vmmax] = Dict(i => data["bus"][i]["vmax"] for i in N) # maximum voltage magnitude
    vamin = m.ext[:parameters][:vamin] = Dict(i => -pi for i in N) # Arbitrary limit of -pi for minimum bus voltage angle
    vamax = m.ext[:parameters][:vamax] = Dict(i =>  pi for i in N) # Arbitrary limit of  pi for maximum bus voltage angle

    # Bus paramters DC
    vmmin_dc = m.ext[:parameters][:vmmin_dc] = Dict(i => data["busdc_ne"][i]["Vdcmin"] for i in N_dc) # minimum voltage magnitude
    vmmax_dc = m.ext[:parameters][:vmmax_dc] = Dict(i => data["busdc_ne"][i]["Vdcmax"] for i in N_dc) # maximum voltage magnitude


    # Branch parameters
    rb = m.ext[:parameters][:rb] = Dict(b => (data["branch"][b]["br_r"]) for b in B) # branch resistance
    xb = m.ext[:parameters][:xb] = Dict(b => (data["branch"][b]["br_x"])  for b in B) # branch reactance
    gb =  m.ext[:parameters][:gb] = Dict(b => real(1 / (data["branch"][b]["br_r"] + data["branch"][b]["br_x"]im)) for b in B) # branch series conductance
    bb =  m.ext[:parameters][:bb] = Dict(b => imag(1 / (data["branch"][b]["br_r"] + data["branch"][b]["br_x"]im)) for b in B) # branch series admittance
    gfr = m.ext[:parameters][:gb_sh_fr] = Dict(b => data["branch"][b]["g_fr"] for b in B) # branch shunt conductance from side i -> j
    bfr = m.ext[:parameters][:bb_sh_fr] = Dict(b => data["branch"][b]["b_fr"] for b in B) # branch shunt susceptance from side i -> j
    gto = m.ext[:parameters][:gb_sh_to] = Dict(b => data["branch"][b]["g_to"] for b in B) # branch shunt conductance to side j -> i
    bto = m.ext[:parameters][:bb_sh_to] = Dict(b => data["branch"][b]["b_to"] for b in B) # branch shunt susceptance to side  j -> j
    # smax = m.ext[:parameters][:smax] = Dict(b => data["branch"][b]["rate_a"] for b in B) # branch rated power in pu
    smax = m.ext[:parameters][:smax] = Dict(b => (data["branch"][b]["rate_a"])*0.3 for b in B) # branch rated power in pu
    # smax = m.ext[:parameters][:smax] = Dict(b => 0.4 for b in B) # branch rated power in pu
    imax = m.ext[:parameters][:imax] = Dict(b => data["branch"][b]["c_rating_a"] for b in B if haskey(data["branch"][b],"c_rating_a")) # branch rated power in pu
    angmin = m.ext[:parameters][:angmin] = Dict(b => data["branch"][b]["angmin"] for b in B) # minimum voltage angle difference over branch
    angmax = m.ext[:parameters][:angmax] = Dict(b => data["branch"][b]["angmax"] for b in B) # maximum voltage angle difference over branch
    b_shift = m.ext[:parameters][:b_shift] = Dict(b => data["branch"][b]["shift"] for b in B) # maximum voltage angle difference over branch
    b_tap = m.ext[:parameters][:b_tap] = Dict(b => data["branch"][b]["tap"] for b in B)
        
    # Branch parameters DC
    smax_dc = m.ext[:parameters][:smax_dc] = Dict(b => (data["branchdc_ne"][b]["rateA"])/baseMVA for b in B_dc) # branch rated power in pu
    rd = m.ext[:parameters][:rd] = Dict(b => data["branchdc_ne"][b]["r"] for b in B_dc) # branch dc resistance
    branch_cost = m.ext[:parameters][:branch_cost] =  Dict(b => data["branchdc_ne"][b]["cost"] for b in B_dc) # cost branch
    
    # Transformer parameter
    rtf = m.ext[:parameters][:rtf] = Dict(n => data["convdc_ne"][n]["rtf"] for n in N_tf) # branch reactance
    xtf = m.ext[:parameters][:xtf] = Dict(n => (data["convdc_ne"][n]["xtf"]) for n in N_tf) # branch reactance
    gtf =  m.ext[:parameters][:gtf] = Dict(n => real(1 / (data["convdc_ne"][n]["rtf"] + data["convdc_ne"][n]["xtf"]im)) for n in N_tf) # branch series conductance
    btf =  m.ext[:parameters][:btf] = Dict(n => imag(1 / (data["convdc_ne"][n]["rtf"] + data["convdc_ne"][n]["xtf"]im)) for n in N_tf) # branch series admittance


    tc_tf = m.ext[:parameters][:tc_tf] = Dict(n => data["convdc_ne"][n]["tm"] for n in N_tf) # Transformer winding ratio

    vmax_tf =  m.ext[:parameters][:vmax_tf] = Dict(n => data["convdc_ne"][n]["Vmmax"] for n in N_tf) # Maximum AC active Power
    vmin_tf =  m.ext[:parameters][:vmin_tf] = Dict(n => data["convdc_ne"][n]["Vmmin"] for n in N_tf) # Maximum AC active Power

    pmax_tf = m.ext[:parameters][:pmax_tf] = Dict(n => (data["convdc_ne"][n]["Pacmax"])/baseMVA for n in N_tf) # Maximum AC active Power
    pmin_tf = m.ext[:parameters][:pmin_tf] = Dict(n => (data["convdc_ne"][n]["Pacmin"])/baseMVA for n in N_tf) # Minumun AC active Power 
    qmax_tf = m.ext[:parameters][:qmax_tf] = Dict(n => (data["convdc_ne"][n]["Qacmax"])/baseMVA for n in N_tf) # Maximum AC reactive Power
    qmin_tf = m.ext[:parameters][:qmin_tf] = Dict(n => (data["convdc_ne"][n]["Qacmin"])/baseMVA for n in N_tf) # Minumun AC reactive Power 

    # converter 

    a_loss_cv = m.ext[:parameters][:a_loss_cv] = Dict(n => data["convdc_ne"][n]["LossA"]/baseMVA for n in N_tf) # zero order losses
    b_loss_cv = m.ext[:parameters][:b_loss_cv] = Dict(n => data["convdc_ne"][n]["LossB"]/baseMVA for n in N_tf) # first order losses
    c_loss_cv = m.ext[:parameters][:c_loss_cv] = Dict(n => data["convdc_ne"][n]["LossCrec"]/baseMVA for n in N_tf) # second order losses
    
    i_cv_lim = m.ext[:parameters][:i_cv_lim] = Dict(n => data["convdc_ne"][n]["Imax"] for n in N_tf) # zero order losses
    cv_cost = m.ext[:parameters][:cv_cost] = Dict(n => data["convdc_ne"][n]["cost"] for n in N_tf) # cost converter
    # Phase reactor

    rc = m.ext[:parameters][:rc] = Dict(n => data["convdc_ne"][n]["rc"] for n in N_tf) # phase reactance
    xc = m.ext[:parameters][:xc] = Dict(n => data["convdc_ne"][n]["xc"] for n in N_tf) # phase reactance
    gc =  m.ext[:parameters][:gc] = Dict(n => real(1 / (data["convdc_ne"][n]["rc"] + data["convdc_ne"][n]["xc"]im))*Zbase  for n in N_tf) # phase conductance
    bc =  m.ext[:parameters][:bc] = Dict(n => imag(1 / (data["convdc_ne"][n]["rc"] + data["convdc_ne"][n]["xc"]im))*Zbase for n in N_tf) # phase admittance


    # Filter
    
    bf = m.ext[:parameters][:bf] = Dict(n => data["convdc_ne"][n]["bf"] for n in N_tf) # branch resistance



    # Load parameters: Assuming a fixed demand!
    pd = m.ext[:parameters][:pd] = Dict(l => data["load"][l]["pd"] for l in L)  # active power demand in pu
    qd = m.ext[:parameters][:qd] = Dict(l => data["load"][l]["qd"] for l in L)  # reactive power demand in pu
    # pd = m.ext[:parameters][:pd] = Dict(l => 0 for l in L)  # active power demand in pu
    # qd = m.ext[:parameters][:qd] = Dict(l => 2 for l in L)  # reactive power demand in pu

    # il_rated = m.ext[:parameters][:il_rated] = Dict(l => (data["load"][l]["pd"]^2 + data["load"][l]["qd"]^2) / 0.9 for l in L)  # active power demand in pu
    
    # Shunt elements
    gs =  m.ext[:parameters][:gs] = Dict(s => data["shunt"][s]["gs"] for s in S) # branch series conductance
    bs =  m.ext[:parameters][:bs] = Dict(s => data["shunt"][s]["bs"] for s in S) # branch series admittance

    # Generator parameters
    pmax = m.ext[:parameters][:pmax] = Dict(g => data["gen"][g]["pmax"] for g in G)  # maximum active power in pu
    # pmax = m.ext[:parameters][:pmax] = Dict(g => 100 for g in G)  # maximum active power in pu
    pmin = m.ext[:parameters][:pmin] = Dict(g => data["gen"][g]["pmin"] for g in G)  # minimum active power in pu
    qmax = m.ext[:parameters][:qmax] = Dict(g => data["gen"][g]["qmax"] for g in G)  # maximum reactive power in pu
    qmin = m.ext[:parameters][:qmin] = Dict(g => data["gen"][g]["qmin"] for g in G)  # minimum reactive power in pu
    ig_rated = m.ext[:parameters][:ig_rated] = Dict(g => sqrt(data["gen"][g]["pmax"]^2 + data["gen"][g]["qmax"]^2) / 0.9 for g in G)  # maximum active power in pu
    
    max_gen_ncost = m.ext[:parameters][:gen_max_ncost] = maximum([data["gen"][g]["ncost"] for g in G])
    m.ext[:parameters][:gen_ncost] = Dict(g => data["gen"][g]["ncost"] for g in G)
    m.ext[:parameters][:gen_cost] = Dict(g => data["gen"][g]["cost"] for g in G)
    # Uniform the length of the cost vector for all generators
    for (gen_id,gen_cost) in m.ext[:parameters][:gen_cost]
        while (length(m.ext[:parameters][:gen_cost][gen_id]) < max_gen_ncost)
            prepend!(m.ext[:parameters][:gen_cost][gen_id],0)
        end
    end            
    
    return m
end