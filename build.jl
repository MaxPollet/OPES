function build_ac_opf!(m::Model)
    # This function builds the polar form of a nonlinear AC power flow formulation

    # Create m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    N = m.ext[:sets][:N]
    N_sl = m.ext[:sets][:N_sl]
    B = m.ext[:sets][:B]
    B_ac_fr = m.ext[:sets][:B_ac_fr]
    B_ac_to = m.ext[:sets][:B_ac_to]
    G = m.ext[:sets][:G]
    G_ac = m.ext[:sets][:G_ac]
    L = m.ext[:sets][:L]
    L_ac = m.ext[:sets][:L_ac]
    B_ac = m.ext[:sets][:B_ac]
    B_arcs = m.ext[:sets][:B_arcs]
    S = m.ext[:sets][:S]
    S_ac = m.ext[:sets][:S_ac]

    N_dc = m.ext[:sets][:N_dc]
    B_dc = m.ext[:sets][:B_dc]
    B_dc_to = m.ext[:sets][:B_dc_to]
    B_dc_fr = m.ext[:sets][:B_dc_fr]
    B_dc_to_fr = m.ext[:sets][:B_dc_to_fr]
    N_tf = m.ext[:sets][:N_tf]
    B_tf_dc_ac_fr = m.ext[:sets][:B_tf_dc_ac_fr]
    B_tf_dc_ac_to = m.ext[:sets][:B_tf_dc_ac_to]
    B_tf_dc_ac_to_fr = m.ext[:sets][:B_tf_dc_ac_to_fr] 

    # Extract parameters
    vmmin = m.ext[:parameters][:vmmin]
    vmmax = m.ext[:parameters][:vmmax]
    vamin = m.ext[:parameters][:vamin]
    vamax = m.ext[:parameters][:vamax]
    rb =  m.ext[:parameters][:rb]
    xb =  m.ext[:parameters][:xb] 
    gb =  m.ext[:parameters][:gb]
    bb =  m.ext[:parameters][:bb] 
    gs =  m.ext[:parameters][:gs]
    bs =  m.ext[:parameters][:bs] 
    gfr = m.ext[:parameters][:gb_sh_fr]
    bfr = m.ext[:parameters][:bb_sh_fr]
    gto = m.ext[:parameters][:gb_sh_to]
    bto = m.ext[:parameters][:bb_sh_to]
    smax = m.ext[:parameters][:smax]
    angmin = m.ext[:parameters][:angmin]
    angmax = m.ext[:parameters][:angmax]
    b_shift = m.ext[:parameters][:b_shift]
    b_tap = m.ext[:parameters][:b_tap]
    pd = m.ext[:parameters][:pd]
    qd = m.ext[:parameters][:qd]
    pmax = m.ext[:parameters][:pmax]
    pmin = m.ext[:parameters][:pmin]
    qmax = m.ext[:parameters][:qmax]
    qmin = m.ext[:parameters][:qmin]
    gen_cost = m.ext[:parameters][:gen_cost]

    vmmin_dc = m.ext[:parameters][:vmmin_dc]# minimum voltage magnitude
    vmmax_dc = m.ext[:parameters][:vmmax_dc]# maximum voltage magnitude

    #branch DC parameters
    smax_dc = m.ext[:parameters][:smax_dc] # branch rated power in pu
    rd = m.ext[:parameters][:rd] # branch dc resistance
    branch_cost = m.ext[:parameters][:branch_cost] =  Dict(b => data["branchdc_ne"][b]["cost"] for b in B_dc) # cost branch

    # Transformer parameter
    rtf = m.ext[:parameters][:rtf] # branch reactance
    xtf = m.ext[:parameters][:xtf] # branch reactance
    gtf =  m.ext[:parameters][:gb] # branch series conductance
    btf =  m.ext[:parameters][:bb] # branch series admittance


    tc_tf = m.ext[:parameters][:tc_tf]  # Transformer winding ratio

    vmax_tf =  m.ext[:parameters][:vmax_tf] # Maximum AC active Power
    vmin_tf =  m.ext[:parameters][:vmin_tf] # Maximum AC active Power


    pmax_tf =  m.ext[:parameters][:pmax_tf] # Maximum AC active Power
    pmin_tf = m.ext[:parameters][:pmin_tf] # Minumun AC active Power 
    qmax_tf =  m.ext[:parameters][:qmax_tf] # Maximum AC reactive Power
    qmin_tf = m.ext[:parameters][:qmin_tf]  # Minumun AC reactive Power 


    # Phase reactor
    rc = m.ext[:parameters][:rc]  # phase reactance
    xc = m.ext[:parameters][:xc] # phase reactance
    gc =  m.ext[:parameters][:gc] # phase conductance
    bc =  m.ext[:parameters][:bc] # phase admittance


    # Filter  
    bf = m.ext[:parameters][:bf]  # branch resistance

    # converter 
    i_cv_lim = m.ext[:parameters][:i_cv_lim] # current limit converter
    a_loss_cv = m.ext[:parameters][:a_loss_cv] # zero order losses
    b_loss_cv = m.ext[:parameters][:b_loss_cv] # first order losses
    c_loss_cv = m.ext[:parameters][:c_loss_cv]  # second order losses

    cv_cost = m.ext[:parameters][:cv_cost] # cost converter

    ##### Create variables 
    # integer variables
    ed = m.ext[:variables][:ed] = @variable(m, [n = N_dc], binary=true, base_name="DC cable") # DC on off parameter
    ec = m.ext[:variables][:ec] = @variable(m, [n = N_tf], binary=true, base_name="converter")   # Converter on off parameter

    # Bus variables
    vm = m.ext[:variables][:vm] = @variable(m, [i=N], lower_bound = vmmin[i], upper_bound = vmmax[i], base_name = "vm") # voltage magnitude
    va = m.ext[:variables][:va] = @variable(m, [i=N], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va") # voltage angle
    # Bus variables DC
    vm_dc = m.ext[:variables][:vm_dc] = @variable(m, [i=N_dc], lower_bound = vmmin_dc[i], upper_bound = vmmax_dc[i], base_name = "vm_dc") # voltage magnitude
    # voltage Transformer
    vm_tf = m.ext[:variables][:vm_tf] = @variable(m, [i=N_tf], lower_bound = vmmin[i], upper_bound = vmmax[i], base_name = "vm_tf") # voltage magnitude
    va_tf = m.ext[:variables][:va_tf] = @variable(m, [i=N_tf], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_tf") # voltage angle


    #voltage phase reactor
    vm_ph = m.ext[:variables][:vm_ph] = @variable(m, [i=N_tf], lower_bound = vmmin[i], upper_bound = vmmax[i], base_name = "vm_ph") # voltage magnitude
    va_ph = m.ext[:variables][:va_ph] = @variable(m, [i=N_tf], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_ph") # voltage angle


    # # voltage contverter

    vm_cv = m.ext[:variables][:vm_cv] = @variable(m, [i=N_tf], lower_bound = vmmin[i], upper_bound = vmmax[i], base_name = "vm_cv") # voltage magnitude
    va_cv = m.ext[:variables][:va_cv] = @variable(m, [i=N_tf], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_cv") # voltage angle

    # Generator variables
    pg = m.ext[:variables][:pg] = @variable(m, [g=G], lower_bound = pmin[g], upper_bound = pmax[g], base_name = "pg") # active and reactive
    qg = m.ext[:variables][:qg] = @variable(m, [g=G], lower_bound = qmin[g], upper_bound = qmax[g], base_name = "qg") # voltage angle

    # Branch variables
    pb = m.ext[:variables][:pb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "pb") # from side active power flow (i->j)
    qb = m.ext[:variables][:qb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "qb") # from side reactive power flow (i->j)
    # Branch variables DC
    pb_dc = m.ext[:variables][:pb_dc] = @variable(m, [(b,i,j) in B_dc_to_fr], lower_bound = -smax_dc[b], upper_bound = smax_dc[b], base_name = "pb_dc") # from side active power flow (i->j)

    # # Transformer
    pb_tf = m.ext[:variables][:pb_tf] = @variable(m, [(b,i,j) = B_tf_dc_ac_to_fr], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_tf") # from side active power flow (i->j)
    qb_tf = m.ext[:variables][:qb_tf] = @variable(m, [(b,i,j) = B_tf_dc_ac_to_fr], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_tf") # from side reactive power flow (i->j)
     


    # phase 
    pb_ph = m.ext[:variables][:pb_ph] = @variable(m, [(b,i,j) = B_tf_dc_ac_to_fr], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_ph") # from side active power flow (i->j)
    qb_ph = m.ext[:variables][:qb_ph] = @variable(m, [(b,i,j) = B_tf_dc_ac_to_fr], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_ph") # from side reactive power flow (i->j)
 
    # converter current 

    i_cv = m.ext[:variables][:i_cv] = @variable(m, [n = N_tf], lower_bound = -i_cv_lim[n], upper_bound = i_cv_lim[n], base_name = "i_cv") # current of converter


    # DC node power
    pb_dc_node = m.ext[:variables][:pb_dc_node] = @variable(m, [n = N_dc], base_name = "pb_dc_node") # current of converter

    ##### Objective
    # max_gen_ncost = m.ext[:parameters][:gen_max_ncost]
    # if max_gen_ncost == 1
    #     m.ext[:objective] = @objective(m, Min,
    #             sum(gen_cost[g][1]
    #                     for g in G)
    #     )
    # elseif max_gen_ncost == 2
    #     m.ext[:objective] = @objective(m, Min,
    #             sum(gen_cost[g][1]*pg[g] + gen_cost[g][2]
    #                     for g in G)
    #     )
    # elseif max_gen_ncost == 3
    #     m.ext[:objective] = @NLobjective(m, Min,
    #             sum(gen_cost[g][1]*pg[g]^2 + gen_cost[g][2]*pg[g] + gen_cost[g][3]
    #                     for g in G)
    #     )
    # elseif max_gen_ncost == 4
    #     m.ext[:objective] = @NLobjective(m, Min,
    #             sum(gen_cost[g][1]*pg[g]^3 + gen_cost[g][2]*pg[g]^2 + gen_cost[g][3]*pg[g] + gen_cost[g][4]
    #                     for g in G)
    #     )
    # end
    
    m.ext[:objective] = @NLobjective(m, Min, sum(cv_cost[n]*ec[n] for n in N_tf) + sum(branch_cost[b]*ed[b] for (b,i,j) in B_dc_to))


    #on off contstrains
    
    m.ext[:constraints][:pd_dc_upper] = @constraint(m, [(b,i,j) in B_dc_to_fr],
        pb_dc[(b,i,j)] <=  ed[b]*smax_dc[b]
    )

    m.ext[:constraints][:pd_dc_lower] = @constraint(m, [(b,i,j) in B_dc_to_fr],
        pb_dc[(b,i,j)] >=  ed[b]*-smax_dc[b]
    )

    m.ext[:constraints][:pd_tf_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_tf[(b,i,j)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pd_tf_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_tf[(b,i,j)] >=  ec[b]*pmin_tf[b]
    )

    m.ext[:constraints][:qd_tf_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_tf[(b,i,j)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qd_tf_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_tf[(b,i,j)] >=  ec[b]*qmin_tf[b]
    )

    m.ext[:constraints][:i_cv_uper] = @constraint(m, [n = N_tf],
        i_cv[n] <=  ec[n]*i_cv_lim[n]
    )
    
    m.ext[:constraints][:i_cv_lower] = @constraint(m, [n = N_tf],
        i_cv[n] >=  ec[n]*-i_cv_lim[n]
    )

    m.ext[:constraints][:pb_ph_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_ph[(b,i,j)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pb_ph_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        pb_ph[(b,i,j)] >=  ec[b]*pmin_tf[b]
    )

    m.ext[:constraints][:qb_ph_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        qb_ph[(b,i,j)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qb_ph_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac_to_fr],
        qb_ph[(b,i,j)] >=  ec[b]*qmin_tf[b]
    )

    

    m.ext[:constraints][:pb_dc_node_lower] = @constraint(m, [n in N_dc],
        pb_dc_node[n] >=  ec[n]*(qmin_tf[n]^2 + pmin_tf[n]^2)
    )

    m.ext[:constraints][:pb_dc_node_lower] = @constraint(m, [n in N_dc],
        pb_dc_node[n] >=  ec[n]*(qmax_tf[n]^2 + pmax_tf[n]^2)
    )


    # Power flow constraints in from and to direction
    m.ext[:constraints][:pbij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] ==  (gb[b] + gfr[b])*vm[i]^2/b_tap[b]^2 - (gb[b] * vm[i] * vm[j] * cos(va[i] - va[j] - b_shift[b]))/b_tap[b] - (bb[b] * vm[i] * vm[j] * sin(va[i] - va[j] - b_shift[b]))/b_tap[b]) # active power i to j
    m.ext[:constraints][:qbij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == -(bb[b] + bfr[b])*vm[i]^2/b_tap[b]^2 + (bb[b] * vm[i] * vm[j] * cos(va[i] - va[j] - b_shift[b]))/b_tap[b] - (gb[b] * vm[i] * vm[j] * sin(va[i] - va[j] - b_shift[b]))/b_tap[b]) # reactive power i to j
    m.ext[:constraints][:pbji] = @NLconstraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] ==  (gb[b] + gto[b])*(vm[j])^2 - (gb[b] * (vm[j]) * vm[i] * cos(va[j] - va[i] + b_shift[b]))/b_tap[b] - (bb[b] * vm[j] * vm[i] * sin(va[j] - va[i] + b_shift[b]))/b_tap[b]) # active power j to i
    m.ext[:constraints][:qbji] = @NLconstraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == -(bb[b] + bto[b])*(vm[j])^2 + (bb[b] * (vm[j]) * vm[i] * cos(va[j] - va[i] + b_shift[b]))/b_tap[b] - (gb[b] * vm[j] * vm[i] * sin(va[j] - va[i] + b_shift[b]))/b_tap[b]) # reactive power j to i

    # Thermal limits for the branches
    m.ext[:constraints][:sij] = @NLconstraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)]^2 + qb[(b, i, j)]^2 <= smax[b]^2)
    m.ext[:constraints][:sji] = @NLconstraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)]^2 + qb[(b, j, i)]^2 <= smax[b]^2)

    # Branch angle limits
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] >= angmin[b])
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] >= angmin[b])

    # DC cables power flow contrains
    m.ext[:constraints][:pbdcij] = @NLconstraint(m, [(n,i,j) = B_dc_to], pb_dc[(n, i, j)] ==  vm_dc[i] * (vm_dc[i] - vm_dc[j])/ rd[n] * ed[n]) # active power i to j
    m.ext[:constraints][:pbdcji] = @NLconstraint(m, [(n,j,i) = B_dc_fr], pb_dc[(n, j, i)] ==  vm_dc[j] * (vm_dc[j] - vm_dc[i])/ rd[n] * ed[n]) # active power i to j

    # transformer power flow contrains
    m.ext[:constraints][:pbtfij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], pb_tf[(n, i, j)] == (gtf[n]*vm[j]^2/tc_tf[n]^2 - (gtf[n] * vm[j] * vm_tf[i]/ tc_tf[n]^2 * cos(va[j] - va_tf[i])) - (btf[n] * vm[j] * vm_tf[i]/ tc_tf[n]^2 * sin(va[j] - va_tf[i])))* ec[n]) # active power i to j
    m.ext[:constraints][:qbtfij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], qb_tf[(n, i, j)] == (-btf[n]*vm[j]^2/tc_tf[n]^2 + (btf[n] * vm[j] * vm_tf[i]/ tc_tf[n]^2 * cos(va[j] - va_tf[i])) - (gtf[n] * vm[j] * vm_tf[i]/ tc_tf[n]^2 * sin(va[j] - va_tf[i]))) * ec[n]) # reactive power i to j
    m.ext[:constraints][:pbtfji] = @NLconstraint(m, [(n,j,i) = B_tf_dc_ac_to], pb_tf[(n, i, j)] == (gtf[n]*vm_tf[i]^2 - (gtf[n] * vm[j] * vm_tf[i]/ tc_tf[n] * cos(va_tf[i] - va[j])) - (btf[n] * vm[j] * vm_tf[i]/ tc_tf[n] * sin(va_tf[i] - va[j])))*ec[n])  # active power j to i
    m.ext[:constraints][:qbtfji] = @NLconstraint(m, [(n,j,i) = B_tf_dc_ac_to], qb_tf[(n, i, j)] == (-btf[n]*vm_tf[i]^2 + (btf[n] * vm[j] * vm_tf[i]/ tc_tf[n] * cos(va_tf[i] - va[j])) - (gtf[n] * vm[j] * vm_tf[i]/ tc_tf[n] * sin(va_tf[i] - va[j])))*ec[n]) # reactive power j to i


    # Filter

    m.ext[:constraints][:pbfij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], pb_tf[(n, j, i)]+ pb_ph[(n, i, j)] == 0) # active power i to j
    m.ext[:constraints][:qbfij] = @NLconstraint(m, [(b,i,j) = B_tf_dc_ac_fr, n = N_tf], qb_tf[(b, j, i)]+ qb_ph[(b, i, j)] - ec[n]*bf[n]* vm_tf[b]^2 == 0) # active power i to j
    
    
    # phase reactor power flow conatrains
    m.ext[:constraints][:pbphij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], pb_ph[(n, i, j)] == (gc[n]*vm_tf[i]^2 - (gc[n] * vm_tf[i] * vm_ph[i] * cos(va_tf[i] - va_ph[i])) - (bc[n] * vm_tf[i] * vm_ph[i]* sin(va_tf[i] - va_ph[i])))* ec[n]) # active power i to j
    m.ext[:constraints][:qbphij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], qb_ph[(n, i, j)] == (-bc[n]*vm_tf[i]^2 + (bc[n] * vm_tf[i] * vm_ph[i]* cos(va_tf[i] - va_ph[i])) - (gc[n] * vm_tf[i] * vm_ph[i] * sin(va_tf[i] - va_ph[i])))*ec[n]) # reactive power i to j
    m.ext[:constraints][:pbphji] = @NLconstraint(m, [(n,j,i) = B_tf_dc_ac_to], pb_ph[(n, j, i)] == (gc[n]*vm_ph[i]^2 - (gc[n] * vm_tf[i] * vm_ph[i] * cos(va_ph[i] - va_tf[i])) - (bc[n] * vm_tf[i] * vm_ph[i] * sin(va_ph[i] - va_tf[i])))*ec[n])  # active power j to i
    m.ext[:constraints][:qbphji] = @NLconstraint(m, [(n,j,i) = B_tf_dc_ac_to], qb_ph[(n, j, i)] == (-bc[n]*vm_ph[i]^2 + (bc[n] * vm_tf[i] * vm_ph[i] * cos(va_ph[i] - va_tf[i])) - (gc[n] * vm_tf[i] * vm_ph[i] * sin(va_ph[i] - va_tf[i])))*ec[n]) # reactive power j to i

    # convertere power flow constraints
    m.ext[:constraints][:pbcvij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], pb_dc_node[i] -  pb_ph[(n,j,i)] == a_loss_cv[n]* ec[n] + b_loss_cv[n] * i_cv[n] + c_loss_cv[n] * i_cv[n]^2) # 
    m.ext[:constraints][:pbcvij] = @NLconstraint(m, [(n,i,j) = B_tf_dc_ac_fr], pb_ph[(n,j,i)]^2 +  qb_ph[(n,i,j)]^2 == vm_cv[n]^2 * i_cv[n]^2) # 

    # DC node power flow contrains 
    m.ext[:constraints][:p_balance] = @constraint(m, [i in N], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )


    # Kirchhoff's current law, i.e., nodal power balance
    if isempty(S)
        m.ext[:constraints][:p_balance] = @constraint(m, [i in N], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
        m.ext[:constraints][:q_balance] = @constraint(m, [i in N], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    else
        m.ext[:constraints][:p_balance] = @NLconstraint(m, [i in N], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) - sum(gs[s]*vm[i]^2 for s in S_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
        m.ext[:constraints][:q_balance] = @NLconstraint(m, [i in N], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) + sum(bs[s]*vm[i]^2 for s in S_ac[i]) == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    end
    # Voltage angle on reference bus = 0, reference bus is bus 4 in this case
    m.ext[:constraints][:varef] = @constraint(m, [n_sl in N_sl], va[n_sl] == 0)
    
    return m 
end
