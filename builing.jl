function build_ac_opf!(m::Model)

    # This function builds the polar form of a linear AC power flow formulation

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
    B_tf_dc_ac = m.ext[:sets][:B_tf_dc_ac]
    B_tf_ac_dc = m.ext[:sets][:B_tf_ac_dc]

    B_dc_arcs = m.ext[:sets][:B_dc_arcs] 
    B_dc_aux_to_fr = m.ext[:sets][:B_dc_aux_to_fr]

    
    N_dc_index = m.ext[:sets][:N_dc_index]
    N_dc_ac_connected = m.ext[:sets][:N_dc_ac_connected]
    N_no_dc_connected = m.ext[:sets][:N_no_dc_connected]
    N_all_dc_connected = m.ext[:sets][:N_all_dc_connected]



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
    vmmin_dc = m.ext[:parameters][:vmmin_dc]
    vmmax_dc = m.ext[:parameters][:vmmax_dc]

    #branch DC parameters
    smax_dc = m.ext[:parameters][:smax_dc] # branch rated power in pu
    rd = m.ext[:parameters][:rd] # branch dc resistance
    branch_cost = m.ext[:parameters][:branch_cost] =  Dict(b => data["branchdc_ne"][b]["cost"] for b in B_dc) # cost branch
    vmax_dc =  m.ext[:parameters][:vmax_dc]
    vmin_dc =  m.ext[:parameters][:vmin_dc]

    # Transformer parameter
    rtf = m.ext[:parameters][:rtf] # branch resistance
    xtf = m.ext[:parameters][:xtf] # branch reactance
    gtf =  m.ext[:parameters][:gtf] # branch series conductance
    btf =  m.ext[:parameters][:btf] # branch series admittance


    tc_tf = m.ext[:parameters][:tc_tf]  # Transformer winding ratio

    vmax_tf =  m.ext[:parameters][:vmax_tf] # Maximum transformer voltage magnitude 
    vmin_tf =  m.ext[:parameters][:vmin_tf] # Minimum transformer voltage magnitude


    pmax_tf = m.ext[:parameters][:pmax_tf] # Maximum transformer active Power
    pmin_tf = m.ext[:parameters][:pmin_tf] # Minumun transformer active Power 
    qmax_tf = m.ext[:parameters][:qmax_tf] # Maximum transformer reactive Power
    qmin_tf = m.ext[:parameters][:qmin_tf]  # Minumun transformer reactive Power 


    # Phase reactor
    rc = m.ext[:parameters][:rc]  # phase resistance
    xc = m.ext[:parameters][:xc] # phase reactance
    gc =  m.ext[:parameters][:gc] # phase conductance
    bc =  m.ext[:parameters][:bc] # phase admittance


    # Filter  
    bf = m.ext[:parameters][:bf]  # branch resistance

    # converter 
    i_cv_lim = m.ext[:parameters][:i_cv_lim] # current limit converter
    a_loss_cv = m.ext[:parameters][:a_loss_cv] # zero order losses
    b_loss_cv = m.ext[:parameters][:b_loss_cv] # first order losses

    cv_cost = m.ext[:parameters][:cv_cost] # cost converter

     

 

    ##### Create variables 
    # integer variables
    ed = m.ext[:variables][:ed] = @variable(m, [n = B_dc], binary=true, base_name="DC cable") # DC on off parameter
    ec = m.ext[:variables][:ec] = @variable(m, [n = N_tf], binary=true, base_name="converter")   # Converter on off parameter

    ei = m.ext[:variables][:ei] = @variable(m, [n = N_tf], binary=true, base_name="I_cv abs auxilary") # DC on off parameter


    # # Bus variables
    phi_i = m.ext[:variables][:phi_i] = @variable(m, [n = N], base_name = "phi_i") #  phi i variable
    phi_i_aux = m.ext[:variables][:phi_i_aux] = @variable(m, [n = N], base_name = "phi_i_aux") #  phi i auxilary variable
    va = m.ext[:variables][:va] = @variable(m, [i=N], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va") # voltage angle
    va_aux = m.ext[:variables][:va_aux] = @variable(m, [i=N], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_aux") # voltage auxilary angle


    # Bus variables DC
    phi_dc = m.ext[:variables][:phi_dc] = @variable(m, [i=N_dc], lower_bound =  vmin_dc[i] - 1,upper_bound =  vmax_dc[i] - 1, base_name = "phi_dc") # voltage magnitude
    phi_dc_aux = m.ext[:variables][:phi_dc_aux] = @variable(m, [(n,i) = B_dc_aux_to_fr], vmin_dc[i] - 1,upper_bound =  vmax_dc[i] - 1, base_name = "phi_dc_aux") #auxilary voltage magnitude
   
    # voltage Transformer
    va_tf = m.ext[:variables][:va_tf] = @variable(m, [i=N_tf], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_tf") # voltage angle
    va_tf_aux = m.ext[:variables][:va_tf_aux] = @variable(m, [i=N_tf], lower_bound = vamin[i], upper_bound = vamax[i], base_name = "va_tf_aux") # voltage auxilary angle
    phi_tf = m.ext[:variables][:phi_tf] = @variable(m, [i=N_tf], lower_bound =  vmin_tf[i] - 1, upper_bound = vmax_tf[i] - 1, base_name = "phi_tf") # voltage  magnitude difference
    phi_tf_aux = m.ext[:variables][:phi_tf_aux] = @variable(m, [i=N_tf], lower_bound = vmin_tf[i] - 1, upper_bound = vmax_tf[i] - 1, base_name = "phi_tf_aux") # auxilary voltage magnitude difference


    #voltage phase reactor
    va_ph = m.ext[:variables][:va_ph] = @variable(m, [i=N_tf], lower_bound = vmin_tf[i] - 1, upper_bound = vmax_tf[i] - 1, base_name = "va_ph") # voltage angle
    phi_ph = m.ext[:variables][:phi_ph] = @variable(m, [i=N_tf], lower_bound = vmin_tf[i] - 1, upper_bound = vmax_tf[i] - 1, base_name = "phi_ph") # voltage magnitude difference



    # Generator variables
    pg = m.ext[:variables][:pg] = @variable(m, [g=G], lower_bound = pmin[g], upper_bound = pmax[g], base_name = "pg") # active and reactive
    qg = m.ext[:variables][:qg] = @variable(m, [g=G], lower_bound = qmin[g], upper_bound = qmax[g], base_name = "qg") # voltage angle

    # # Branch variables
    pb = m.ext[:variables][:pb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "pb") # from side active power flow (i->j)
    qb = m.ext[:variables][:qb] = @variable(m, [(b,i,j) in B_ac], lower_bound = -smax[b], upper_bound = smax[b], base_name = "qb") # from side reactive power flow (i->j)
   
    # # Branch variables DC
    pb_dc = m.ext[:variables][:pb_dc] = @variable(m, [(b,i,j) in B_dc_to_fr], lower_bound = -smax_dc[b], upper_bound = smax_dc[b], base_name = "pb_dc") # from side active power flow (i->j)

    # # Transformer
 
    pb_tf_dc_ac = m.ext[:variables][:pb_tf_dc_ac] = @variable(m, [(b,i,j) = B_tf_dc_ac], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_tf_dc_ac") # from side active power flow (i->j)
    pb_tf_ac_dc = m.ext[:variables][:pb_tf_ac_dc] = @variable(m, [(b,i,j) = B_tf_ac_dc], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_tf_ac_dc") # from side active power flow (i->j)
    qb_tf_dc_ac = m.ext[:variables][:qb_tf_dc_ac] = @variable(m, [(b,i,j) = B_tf_dc_ac], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_tf_dc_ac") # from side reactive power flow (i->j)
    qb_tf_ac_dc = m.ext[:variables][:qb_tf_ac_dc] = @variable(m, [(b,i,j) = B_tf_ac_dc], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_tf_ac_dc") # from side reactive power flow (i->j)


    # # phase 

    pb_ph_dc_ac = m.ext[:variables][:pb_ph_dc_ac] = @variable(m, [(b,i,j) = B_tf_dc_ac], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_ph_dc_ac") # from side active power flow (i->j)
    pb_ph_ac_dc = m.ext[:variables][:pb_ph_ac_dc] = @variable(m, [(b,i,j) = B_tf_ac_dc], lower_bound = pmin_tf[b], upper_bound = pmax_tf[b], base_name = "pb_ph_ac_dc") # from side active power flow (i->j)
    qb_ph_dc_ac = m.ext[:variables][:qb_ph_dc_ac] = @variable(m, [(b,i,j) = B_tf_dc_ac], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_ph_dc_ac") # from side reactive power flow (i->j)
    qb_ph_ac_dc = m.ext[:variables][:qb_ph_ac_dc] = @variable(m, [(b,i,j) = B_tf_ac_dc], lower_bound = qmin_tf[b], upper_bound = qmax_tf[b], base_name = "qb_ph_ac_dc") # from side reactive power flow (i->j)

    # # converter current 

    i_cv = m.ext[:variables][:i_cv] = @variable(m, [n = N_tf], lower_bound = -i_cv_lim[n], upper_bound = i_cv_lim[n], base_name = "i_cv") # current of converter


    # # DC node power
    pb_dc_node = m.ext[:variables][:pb_dc_node] = @variable(m, [n = N_dc], base_name = "pb_dc_node") # current of converter



    ###### cosine approx

    cos_aux_b_ij = m.ext[:variables][:cos_aux_b_ij] = @variable(m, [(b,i,j) = B_ac_fr], base_name = "cos_aux_b_ij") # current of converter
    cos_aux_b_ji = m.ext[:variables][:cos_aux_b_ji] = @variable(m, [(b,j,i) = B_ac_to], base_name = "cos_aux_b_ji") # current of converter
    
    cos_aux_tf_ij = m.ext[:variables][:cos_aux_tf_ij] = @variable(m, [(n,i,j) = B_tf_ac_dc], base_name = "cos_aux_tf_ij") # current of converter
    cos_aux_tf_ji = m.ext[:variables][:cos_aux_tf_ji] = @variable(m, [(n,i,j) = B_tf_dc_ac], base_name = "cos_aux_tf_ji") # current of converter
    
    cos_aux_ph_ij = m.ext[:variables][:cos_aux_ph_ij] = @variable(m, [(n,i,j) = B_tf_ac_dc], base_name = "cos_aux_ph_ij") # current of converter
    cos_aux_ph_ji = m.ext[:variables][:cos_aux_ph_ji] = @variable(m, [(n,i,j) = B_tf_dc_ac], base_name = "cos_aux_ph_ji") # current of converter
    


    ###### objective
    
    m.ext[:objective] = @objective(m, Min, sum(cv_cost[n]*ec[n] for n in N_tf) + sum(branch_cost[b]*ed[b] for (b,i,j) in B_dc_to))



    ######### on off contraints, auxilary contrains, cosinus constraints
 
    m.ext[:constraints][:pd_dc_upper] = @constraint(m, [(b,i,j) in B_dc_to_fr],
        pb_dc[(b,i,j)] <=  ed[b]*smax_dc[b]
    )

    m.ext[:constraints][:pd_dc_lower] = @constraint(m, [(b,i,j) in B_dc_to_fr],
        pb_dc[(b,i,j)] >=  ed[b]*(-smax_dc[b])
    )




    m.ext[:constraints][:pd_tf_dc_ac_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        pb_tf_dc_ac[(b,i,j)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pd_tf_dc_ac_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        pb_tf_dc_ac[(b,i,j)] >=  ec[b]*pmin_tf[b]
    )

    m.ext[:constraints][:pd_tf_ac_dc_upper] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        pb_tf_ac_dc[(b,i,j)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pd_tf_ac_dc_lower] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        pb_tf_ac_dc[(b,i,j)] >=  ec[b]*pmin_tf[b]
    )


    m.ext[:constraints][:qd_tf_dc_ac_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        qb_tf_dc_ac[(b,i,j)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qd_tf_dc_ac_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        qb_tf_dc_ac[(b,i,j)] >=  ec[b]*qmin_tf[b]
    )

    m.ext[:constraints][:qd_tf_ac_dc_upper] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        qb_tf_ac_dc[(b,i,j)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qd_tf_ac_dc_lower] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        qb_tf_ac_dc[(b,i,j)] >=  ec[b]*qmin_tf[b]
    )





    m.ext[:constraints][:i_cv_uper] = @constraint(m, [n = N_tf],
        i_cv[n] <=  ec[n]*i_cv_lim[n]
    )
    
    m.ext[:constraints][:i_cv_lower] = @constraint(m, [n = N_tf],
        i_cv[n] >=  ec[n]*-i_cv_lim[n]
    )






    m.ext[:constraints][:pd_ph_dc_ac_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        pb_ph_dc_ac[(b,i,j)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pd_ph_dc_ac_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        pb_ph_dc_ac[(b,i,j)] >=  ec[b]*pmin_tf[b]
    )

    m.ext[:constraints][:pd_ph_ac_dc_upper] = @constraint(m, [(b,j,i) in B_tf_ac_dc],
        pb_ph_ac_dc[(b,j,i)] <=  ec[b]*pmax_tf[b]
    )

    m.ext[:constraints][:pd_ph_ac_dc_lower] = @constraint(m, [(b,j,i) in B_tf_ac_dc],
        pb_ph_ac_dc[(b,j,i)] >=  ec[b]*pmin_tf[b]
    )







    m.ext[:constraints][:qd_ph_dc_ac_upper] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        qb_ph_dc_ac[(b,i,j)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qd_ph_dc_ac_lower] = @constraint(m, [(b,i,j) in B_tf_dc_ac],
        qb_ph_dc_ac[(b,i,j)] >=  ec[b]*qmin_tf[b]
    )

    m.ext[:constraints][:qd_ph_ac_dc_upper] = @constraint(m, [(b,j,i) in B_tf_ac_dc],
        qb_ph_ac_dc[(b,j,i)] <=  ec[b]*qmax_tf[b]
    )

    m.ext[:constraints][:qd_ph_ac_dc_lower] = @constraint(m, [(b,j,i) in B_tf_ac_dc],
        qb_ph_ac_dc[(b,j,i)] >=  ec[b]*qmin_tf[b]
    )

    
    

    m.ext[:constraints][:pb_dc_node_upper] = @constraint(m, [n in N_dc],
        pb_dc_node[n] <=  ec[n]*sqrt(qmax_tf[n]^2 + pmax_tf[n]^2)
    )

    m.ext[:constraints][:pb_dc_node_lower] = @constraint(m, [n in N_dc],
        pb_dc_node[n] >=  -ec[n]*sqrt(qmin_tf[n]^2 + pmin_tf[n]^2)
    )



    m.ext[:constraints][:phi_i_upper] = @constraint(m, [n in N],
        phi_i[n] <=  (vmmax[n] - 1)
    )

    m.ext[:constraints][:phi_i_lower] = @constraint(m, [n in N],
        phi_i[n] >=  (vmmin[n] - 1)
    )




    m.ext[:constraints][:phi_i_aux_upper1] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        phi_i_aux[j_ac] <=  ec[i_dc]*(vmmax[j_ac] - 1)
    )

    m.ext[:constraints][:phi_i_aux_lower1] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        phi_i_aux[j_ac] >=  ec[i_dc]*(vmmin[j_ac] - 1)
    )


    m.ext[:constraints][:phi_i_aux_upper2] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        phi_i_aux[j_ac] <=  phi_i[j_ac] - (1 - ec[i_dc])*(vmmin[j_ac] - 1)
    )

    m.ext[:constraints][:phi_i_aux_lower2] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        phi_i_aux[j_ac] >=  phi_i[j_ac] - (1 - ec[i_dc])*(vmmax[j_ac] - 1)
    )







    m.ext[:constraints][:va_aux_upper1] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        va_aux[j_ac] <=  ec[i_dc]*(vamax[j_ac])
    )

    m.ext[:constraints][:va_aux_lower1] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        va_aux[j_ac] >=  ec[i_dc]*(vamin[j_ac])
    )


    m.ext[:constraints][:va_aux_upper2] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        va_aux[j_ac] <=  va[j_ac] - (1-ec[i_dc])*(vamin[j_ac])
    )

    m.ext[:constraints][:va_aux_lower2] = @constraint(m, [(b,i_dc,j_ac) in B_tf_dc_ac_fr],
        va_aux[j_ac] >=  va[j_ac] - (1-ec[i_dc])*(vamax[j_ac])
    )





    m.ext[:constraints][:phi_tf_upper] = @constraint(m, [n in N_dc],
        phi_tf[n] <=  ec[n]*(vmmax[n] - 1)
    )   

    m.ext[:constraints][:phi_tf_lower] = @constraint(m, [n in N_dc],
        phi_tf[n] >=  ec[n]*(vmmin[n] - 1)
    )




    m.ext[:constraints][:phi_tf_aux_upper1] = @constraint(m, [n in N_dc],
        phi_tf_aux[n] <=  ec[n]*(vmmax[n] - 1)
    )

    m.ext[:constraints][:phi_tf_aux_lower1] = @constraint(m, [n in N_dc],
        phi_tf_aux[n] >=  ec[n]*(vmmin[n] - 1)
    )


    m.ext[:constraints][:phi_tf_aux_upper2] = @constraint(m, [n in N_dc],
        phi_tf_aux[n] <=  phi_tf[n] - (1 - ec[n])*(vmmin[n] - 1)
    )

    m.ext[:constraints][:phi_tf_aux_lower2] = @constraint(m, [n in N_dc],
        phi_tf_aux[n] >=  phi_tf[n] - (1 - ec[n])*(vmmax[n] - 1)
    )

    m.ext[:constraints][:phi_va_upper] = @constraint(m, [n in N_dc],
        va_tf[n] <=  ec[n]*(vamax[n])
    )   

    m.ext[:constraints][:phi_va_lower] = @constraint(m, [n in N_dc],
        va_tf[n] >=  ec[n]*(vamin[n])
    )
    
    
    m.ext[:constraints][:va_tf_aux_upper1] = @constraint(m, [n in N_dc],
        va_tf_aux[n] <=  ec[n]*(vamax[n])
    )

    m.ext[:constraints][:va_tf_aux_lower1] = @constraint(m, [n in N_dc],
        va_tf_aux[n] >=  ec[n]*(vamin[n])
    )


    m.ext[:constraints][:va_tf_aux_upper2] = @constraint(m, [n in N_dc],
        va_tf_aux[n] <=  phi_i[n] - (1 - ec[n])*(vamin[n])
    )   

    m.ext[:constraints][:va_tf_aux_lower2] = @constraint(m, [n in N_dc],
        va_tf_aux[n] >=  phi_i[n] - (1 - ec[n])*(vamax[n])
    )






    m.ext[:constraints][:phi_dc_upper] = @constraint(m, [n in N_dc],
        phi_dc[n] <= (vmmax_dc[n]-1)
    )

    m.ext[:constraints][:phi_dc_lower] = @constraint(m, [n in N_dc],
        phi_dc[n] >=  (vmmin_dc[n]-1)
    )





    m.ext[:constraints][:phi_dc_aux_upper1] = @constraint(m, [(n,i) in B_dc_to],
        phi_dc_aux[(n,i)] <=  ed[n]*(vmmax_dc[i]-1)
    )

    m.ext[:constraints][:phi_dc_aux_lower1] = @constraint(m, [(n,i) in B_dc_to],
        phi_dc_aux[(n,i)] >=  ed[n]*(vmmin_dc[i] - 1)
    )


    m.ext[:constraints][:phi_dc_aux_upper2] = @constraint(m, [(n,i) in B_dc_to],
        phi_dc_aux[(n,i)] <=  phi_dc[i] - (1 - ed[n])*(vmmin_dc[i]-1)
    )

    m.ext[:constraints][:phi_dc_aux_lower2] = @constraint(m, [(n,i) in B_dc_to],
        phi_dc_aux[(n,i)] >=  phi_dc[i] - (1 - ed[n])*(vmmax_dc[i]-1)
    )



    m.ext[:constraints][:cos_aux_b_ij_con] = @constraint(m, [(b,i,j) = B_ac_fr],
        cos_aux_b_ij[(b,i,j)] ==  (va[i] - va[j] - b_shift[b])
    )

    m.ext[:constraints][:cos_aux_b_ji_con] = @constraint(m, [(b,j,i) = B_ac_to],
        cos_aux_b_ji[(b,j,i)] ==  (va[i] - va[j] - b_shift[b])
    )


    m.ext[:constraints][:cos_aux_tf_ij_con] = @constraint(m, [(n,i,j) = B_tf_ac_dc],
        cos_aux_tf_ij[(n,i,j)] ==  (va_aux[i] - va_tf[j])
    )

    m.ext[:constraints][:cos_aux_tf_ji_con] = @constraint(m, [(n,j,i) = B_tf_dc_ac],
        cos_aux_tf_ji[(n,j,i)] ==  (va_aux[j] - va_tf[i])
    )


    m.ext[:constraints][:cos_aux_ph_ij_con] = @constraint(m, [(n,j,i) = B_tf_ac_dc],
        cos_aux_ph_ij[(n,j,i)] ==  (va_tf[j] - va_ph[j])
    )

    m.ext[:constraints][:cos_aux_ph_ji_con] = @constraint(m, [(n,i,j) = B_tf_dc_ac],
        cos_aux_ph_ji[(n,i,j)] ==  (va_tf[i] - va_ph[i])
    )



    m.ext[:constraints][:pd_ph_ac_dc_gen_stop] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        pb_ph_ac_dc[(b,i,j)] + pb_ph_dc_ac[(b,j,i)] >=  0
    )


    m.ext[:constraints][:pd_tf_ac_dc_gen_stop] = @constraint(m, [(b,i,j) in B_tf_ac_dc],
        pb_tf_ac_dc[(b,i,j)] + pb_tf_dc_ac[(b,j,i)] >=  0
    )


    m.ext[:constraints][:pd_tf_ac_dc_gen_stop] = @constraint(m, [(b,i,j) in B_ac],
        pb[(b,i,j)] + pb[(b,j,i)] >=  0
    )




   ####### cosinus approximation corner

    cos_approx_b_ij = m.ext[:variables][:cos_approx_b_ij] = JuMP.@variable(m, [(b,i,j) in B_ac_fr])
    for (b,i,j) in B_ac_fr
        
        f, partition = (x) -> cos(x ), collect(-pi/3:pi/60:pi/3)
        formulation_info_milp = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_b_ij[(b,i,j)] , cos_approx_b_ij[(b,i,j)], partition, false)
    end

    cos_approx_b_ji =  m.ext[:variables][:cos_approx_b_ji] = JuMP.@variable(m, [(b,j,i) in B_ac_to])
    for (b,j,i) in B_ac_to
        
        f, partition = (x) -> cos(x ),collect(-pi/3:pi/60:pi/3)
        formulation_info_milp_er = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_b_ji[(b,j,i)] , cos_approx_b_ji[(b,j,i)], partition, false)
    end



    cos_approx_tf_ij = m.ext[:variables][:cos_approx_tf_ij] = JuMP.@variable(m, [(n,i,j) in B_tf_ac_dc])
    for (n,i,j) in B_tf_ac_dc
        
        f, partition = (x) -> cos(x ),collect(-pi/3:pi/60:pi/3)
        formulation_info_milp = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_tf_ij[(n,i,j)] , cos_approx_tf_ij[(n,i,j)], partition, false)
    end

    cos_approx_tf_ji =  m.ext[:variables][:cos_approx_tf_ji] = JuMP.@variable(m, [(n,i,j) in B_tf_dc_ac])
    for (n,i,j) in B_tf_dc_ac
        
        f, partition = (x) -> cos(x ),collect(-pi/3:pi/60:pi/3)
        formulation_info_milp_er = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_tf_ji[(n,i,j)] , cos_approx_tf_ji[(n,i,j)], partition, false)
    end



    cos_approx_ph_ij = m.ext[:variables][:cos_approx_ph_ij] = JuMP.@variable(m, [(n,i,j) in B_tf_ac_dc])
    for (n,i,j) in B_tf_ac_dc
        
        f, partition = (x) -> cos(x ),collect(-pi/3:pi/60:pi/3)
        formulation_info_milp = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_ph_ij[(n,i,j)] , cos_approx_ph_ij[(n,i,j)], partition, false)
    end

    cos_approx_ph_ji =  m.ext[:variables][:cos_approx_ph_ji] = JuMP.@variable(m, [(n,i,j) in B_tf_dc_ac])
    for (n,i,j) in B_tf_dc_ac
        
        f, partition = (x) -> cos(x ), collect(-pi/3:pi/60:pi/3)
        formulation_info_milp_er = PolyhedralRelaxations.construct_univariate_relaxation!(m, f, cos_aux_ph_ji[(n,i,j)] , cos_approx_ph_ji[(n,i,j)], partition, false)
    end

    ################ Branch AC power flow


    m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] ==  (gb[b] + gfr[b])*(1 + 2 *phi_i[i]) /b_tap[b]^2 - (gb[b] * (cos_approx_b_ij[(b,i,j)] + phi_i[i] + phi_i[j]))/b_tap[b] - (bb[b] * (va[i] - va[j]))/b_tap[b]) # active power i to j
    m.ext[:constraints][:qbij] = @constraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == -(bb[b] + bfr[b])*(1 + 2 *phi_i[i]) /b_tap[b]^2 + (bb[b] * (cos_approx_b_ij[(b,i,j)] + phi_i[i] + phi_i[j]))/b_tap[b] - (gb[b] * (va[i] - va[j]))/b_tap[b]) # reactive power i to j
    m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] ==  (gb[b] + gto[b])*(1 + 2 *phi_i[j]) - (gb[b] * (cos_approx_b_ji[(b,j,i)] + phi_i[j] + phi_i[i]))/b_tap[b] - (bb[b] * (va[j] - va[i]))/b_tap[b]) # active power j to i
    m.ext[:constraints][:qbji] = @constraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == -(bb[b] + bto[b])*(1 + 2 *phi_i[j]) + (bb[b] * (cos_approx_b_ji[(b,j,i)] + phi_i[j] + phi_i[i]))/b_tap[b] - (gb[b] * (va[j]- va[i]))/b_tap[b]) # reactive power j to i



    # m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] ==  (gb[b] + gfr[b])*(1 + 2 *phi_i[i]) /b_tap[b]^2 - (gb[b] * (1 - 0*(va[i] - va[j] - b_shift[b]) + phi_i[i] + phi_i[j]))/b_tap[b] - (bb[b] * (va[i] - va[j]))/b_tap[b]) # active power i to j
    # m.ext[:constraints][:qbij] = @constraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == -(bb[b] + bfr[b])*(1 + 2 *phi_i[i]) /b_tap[b]^2 + (bb[b] * (1 - 0*(va[i] - va[j] - b_shift[b]) + phi_i[i] + phi_i[j]))/b_tap[b] - (gb[b] * (va[i] - va[j]))/b_tap[b]) # reactive power i to j
    # m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] ==  (gb[b] + gto[b])*(1 + 2 *phi_i[j]) - (gb[b] * (1 - 0*(va[i] - va[j] - b_shift[b]) + phi_i[j] + phi_i[i]))/b_tap[b] - (bb[b] * (va[j] - va[i]))/b_tap[b]) # active power j to i
    # m.ext[:constraints][:qbji] = @constraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == -(bb[b] + bto[b])*(1 + 2 *phi_i[j]) + (bb[b] * (1 - 0*(va[i] - va[j] - b_shift[b]) + phi_i[j] + phi_i[i]))/b_tap[b] - (gb[b] * (va[j]- va[i]))/b_tap[b]) # reactive power j to i

    # m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] == - (bb[b] * (va[i] - va[j]))/b_tap[b]) # active power i to j
    # m.ext[:constraints][:qbij] = @constraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == - (gb[b] * (va[i] - va[j]))/b_tap[b]) # reactive power i to j
    # m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] == - (bb[b] * (va[j] - va[i]))/b_tap[b]) # active power j to i
    # m.ext[:constraints][:qbji] = @constraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == - (gb[b] * (va[j] - va[i]))/b_tap[b]) # reactive power j to i

    # m.ext[:constraints][:pbij] = @constraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)] ==  0) # active power i to j
    # m.ext[:constraints][:qbij] = @constraint(m, [(b,i,j) = B_ac_fr], qb[(b, i, j)] == 0) # reactive power i to j
    # m.ext[:constraints][:pbji] = @constraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)] ==  0) # active power j to i
    # m.ext[:constraints][:qbji] = @constraint(m, [(b,j,i) = B_ac_to], qb[(b, j, i)] == 0) # reactive power j to i



    # # Thermal limits for the branches
    m.ext[:constraints][:sij] = @constraint(m, [(b,i,j) = B_ac_fr], pb[(b, i, j)]^2 + qb[(b, i, j)]^2 <= smax[b]^2)
    m.ext[:constraints][:sji] = @constraint(m, [(b,j,i) = B_ac_to], pb[(b, j, i)]^2 + qb[(b, j, i)]^2 <= smax[b]^2)

 
 
    # # DC cables power flow contrains
    m.ext[:constraints][:pbdcij] = @constraint(m, [(n,i,j) = B_dc_to], pb_dc[(n, i, j)] ==  (phi_dc_aux[(n,i)] - phi_dc_aux[(n,j)])/ (rd[n] + 0.001)) # active power i to j
    m.ext[:constraints][:pbdcji] = @constraint(m, [(n,i,j) = B_dc_fr], pb_dc[(n, i, j)] ==  (phi_dc_aux[(n,i)] - phi_dc_aux[(n,j)])/ (rd[n] + 0.001)) # active power i to j

    # m.ext[:constraints][:pbdcij] = @constraint(m, [(n,i,j) = B_dc_to], pb_dc[(n, i, j)] ==  0) # active power i to j
    # m.ext[:constraints][:pbdcji] = @constraint(m, [(n,i,j) = B_dc_fr], pb_dc[(n, i, j)] ==  0) # active power i to j

    ################transformer power flow contraint


    m.ext[:constraints][:pbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_tf_ac_dc[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[i])/tc_tf[n]^2 - gtf[n]*(ec[n]*cos_approx_tf_ij[(n,i,j)] + phi_i_aux[i] + phi_tf[j]) / tc_tf[n] - btf[n]*(va_aux[i] - va_tf[j])  /tc_tf[n])) # active power i to j
    m.ext[:constraints][:qbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_tf_ac_dc[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[i])/tc_tf[n]^2 + btf[n]*(ec[n]*cos_approx_tf_ij[(n,i,j)] + phi_i_aux[i] + phi_tf[j]) / tc_tf[n] - gtf[n]*(va_aux[i] - va_tf[j])  /tc_tf[n])) # reactive power i to j
    m.ext[:constraints][:pbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[j]) - gtf[n]*(ec[n]*cos_approx_tf_ji[(n,i,j)] + phi_i_aux[j] + phi_tf[i]) / tc_tf[n] - btf[n]*(va_tf[i] - va_aux[j])  /tc_tf[n]))  # active power j to i
    m.ext[:constraints][:qbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[j]) + btf[n]*(ec[n]*cos_approx_tf_ji[(n,i,j)] + phi_i_aux[j] + phi_tf[i]) / tc_tf[n] - gtf[n]*(va_tf[i] - va_aux[j])  /tc_tf[n])) # reactive power j to i

    # m.ext[:constraints][:pbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_tf_ac_dc[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[i]) - gtf[n]*(ec[n]*cos_approx_tf_ij[(n,i,j)] + phi_i_aux[i] + phi_tf[j]) - btf[n]*(va_aux[i] - va_tf[j]))) # active power i to j
    # m.ext[:constraints][:qbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_tf_ac_dc[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[i]) + btf[n]*(ec[n]*cos_approx_tf_ij[(n,i,j)] + phi_i_aux[i] + phi_tf[j]) - gtf[n]*(va_aux[i] - va_tf[j]))) # reactive power i to j
    # m.ext[:constraints][:pbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[j]) - gtf[n]*(ec[n]*cos_approx_tf_ji[(n,i,j)] + phi_i_aux[j] + phi_tf[i]) - btf[n]*(va_tf[i] - va_aux[j])))  # active power j to i
    # m.ext[:constraints][:qbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[j]) + btf[n]*(ec[n]*cos_approx_tf_ji[(n,i,j)] + phi_i_aux[j] + phi_tf[i]) - gtf[n]*(va_tf[i] - va_aux[j]))) # reactive power j to i


    

    # m.ext[:constraints][:pbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_tf_ac_dc[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[i])/tc_tf[n]^2 - gtf[n]*(  (ec[n] - 0*(va_aux[i] - va_tf[j])) + phi_i_aux[i] + phi_tf[j]) / tc_tf[n] - btf[n]*(va_aux[i] - va_tf[j])  /tc_tf[n])) # active power i to j
    # m.ext[:constraints][:qbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_tf_ac_dc[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[i])/tc_tf[n]^2 + btf[n]*(  (ec[n] - 0*(va_aux[i] - va_tf[j])) + phi_i_aux[i] + phi_tf[j]) / tc_tf[n] - gtf[n]*(va_aux[i] - va_tf[j])  /tc_tf[n])) # reactive power i to j
    # m.ext[:constraints][:pbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[j]) - gtf[n]*(  (ec[n] - 0*(va_aux[j] - va_tf[i])) + phi_i_aux[j] + phi_tf[i]) / tc_tf[n] - btf[n]*(va_tf[i] - va_aux[j])  /tc_tf[n]))  # active power j to i
    # m.ext[:constraints][:qbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[j]) + btf[n]*(  (ec[n] - 0*(va_aux[j] - va_tf[i])) + phi_i_aux[j] + phi_tf[i]) / tc_tf[n] - gtf[n]*(va_tf[i] - va_aux[j])  /tc_tf[n])) # reactive power j to i


    # m.ext[:constraints][:pbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_tf_ac_dc[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[i]) - gtf[n]*(ec[n] + phi_i_aux[i] + phi_tf[j]) - btf[n]*(va_aux[i] - va_tf[j])  )) # active power i to j
    # m.ext[:constraints][:qbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_tf_ac_dc[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[i]) + btf[n]*(ec[n] + phi_i_aux[i] + phi_tf[j]) - gtf[n]*(va_aux[i] - va_tf[j])  )) # reactive power i to j
    # m.ext[:constraints][:pbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] == ( gtf[n] * (ec[n] + 2*phi_i_aux[j]) - gtf[n]*(ec[n] + phi_i_aux[j] + phi_tf[i]) - btf[n]*(va_tf[i] - va_aux[j])  ))  # active power j to i
    # m.ext[:constraints][:qbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(n, i, j)] == (-btf[n] * (ec[n] + 2*phi_i_aux[j]) + btf[n]*(ec[n] + phi_i_aux[j] + phi_tf[i]) - gtf[n]*(va_tf[i] - va_aux[j])  )) # reactive power j to i



    # m.ext[:constraints][:pbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_tf_ac_dc[(n, i, j)] == 0) # active power i to j
    # m.ext[:constraints][:qbtfacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_tf_ac_dc[(n, i, j)] == 0) # reactive power i to j
    # m.ext[:constraints][:pbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] == 0)  # active power j to i
    # m.ext[:constraints][:qbtfdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(n, i, j)] == 0) # reactive power j to i


    # # # Filter

    m.ext[:constraints][:pbfij] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_tf_dc_ac[(n, i, j)] + pb_ph_ac_dc[(n, j, i)] == 0) # active power tf to ph
    m.ext[:constraints][:qbfij] = @constraint(m, [(b,i,j) = B_tf_dc_ac], qb_tf_dc_ac[(b, i, j)]+ qb_ph_ac_dc[(b, j, i)] - bf[b]*(ec[b] + 2 * phi_tf[b])== 0) # reactive power tf to ph
    
    
    # # # phase reactor power flow conatrains

    m.ext[:constraints][:pbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_ph_ac_dc[(n, i, j)] == ( gc[n] * (ec[n] + 2*phi_tf_aux[j]) - gc[n]*(ec[n]*cos_approx_ph_ij[(n, i, j)] + phi_tf_aux[j] + phi_ph[j]) - bc[n]*(va_tf_aux[j] - va_ph[j])  )) # active power i to j
    m.ext[:constraints][:qbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_ph_ac_dc[(n, i, j)] == (-bc[n] * (ec[n] + 2*phi_tf_aux[j]) + bc[n]*(ec[n]*cos_approx_ph_ij[(n, i, j)] + phi_tf_aux[j] + phi_ph[j]) - gc[n]*(va_tf_aux[j] - va_ph[j]) )) # reactive power i to j
    m.ext[:constraints][:pbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_ph_dc_ac[(n, i, j)] == ( gc[n] * (ec[n] + 2*phi_tf_aux[i]) - gc[n]*(ec[n]*cos_approx_ph_ji[(n, i, j)] + phi_tf_aux[i] + phi_ph[i]) - bc[n]*(va_ph[i] - va_tf_aux[i])  ))  # active power j to i
    m.ext[:constraints][:pbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_ph_dc_ac[(n, i, j)] == (-bc[n] * (ec[n] + 2*phi_tf_aux[i]) + bc[n]*(ec[n]*cos_approx_ph_ji[(n, i, j)] + phi_tf_aux[i] + phi_ph[i]) - gc[n]*(va_ph[i] - va_tf_aux[i])  )) # reactive power j to i


    # m.ext[:constraints][:pbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_ph_ac_dc[(n, i, j)] == ( gc[n] * (ec[n] + 2*phi_tf_aux[j]) - gc[n]*(  (ec[n] - 0*(va_tf[j] - va_ph[j])) + phi_tf_aux[j] + phi_ph[j]) - bc[n]*(va_tf_aux[j] - va_ph[j])  )) # active power i to j
    # m.ext[:constraints][:qbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_ph_ac_dc[(n, i, j)] == (-bc[n] * (ec[n] + 2*phi_tf_aux[j]) + bc[n]*(  (ec[n] - 0*(va_tf[j] - va_ph[j])) + phi_tf_aux[j] + phi_ph[j]) - gc[n]*(va_tf_aux[j] - va_ph[j]) )) # reactive power i to j
    # m.ext[:constraints][:pbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_ph_dc_ac[(n, i, j)] == ( gc[n] * (ec[n] + 2*phi_tf_aux[i]) - gc[n]*(  (ec[n] - 0*(va_tf[i] - va_tf[i])) + phi_tf_aux[i] + phi_ph[i]) - bc[n]*(va_ph[i] - va_tf_aux[i])  ))  # active power j to i
    # m.ext[:constraints][:pbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_ph_dc_ac[(n, i, j)] == (-bc[n] * (ec[n] + 2*phi_tf_aux[i]) + bc[n]*(  (ec[n] - 0*(va_tf[i] - va_tf[i])) + phi_tf_aux[i] + phi_ph[i]) - gc[n]*(va_ph[i] - va_tf_aux[i])  )) # reactive power j to i



    # # m.ext[:constraints][:pbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], pb_ph_ac_dc[(n, i, j)] == 0) # active power i to j
    # # m.ext[:constraints][:qbphacdcij] = @constraint(m, [(n,i,j) = B_tf_ac_dc], qb_ph_ac_dc[(n, i, j)] == 0) # reactive power i to j
    # # m.ext[:constraints][:pbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_ph_dc_ac[(n, i, j)] == 0)  # active power j to i
    # # m.ext[:constraints][:qbphdcacji] = @constraint(m, [(n,i,j) = B_tf_dc_ac], qb_ph_dc_ac[(n, i, j)] == 0) # reactive power j to i



    ####### converter power flow constraint

    m.ext[:constraints][:ei_lower] = @constraint(m, [n in N_tf], (1 - 2*(ei[n]))*(i_cv[n]) >= 0.9*(i_cv[n])) # ei is 1 when I_cv is negative 
    m.ext[:constraints][:ei_upper] = @constraint(m, [n in  N_tf], ei[n] <= ec[n]) # ei is 0 when I_cv is positive 



    # m.ext[:constraints][:ec_ed_con] =  @constraint(m, [i in N_dc_ac_connected], sum(ed[n] for (n,k,j) in B_dc_arcs[i]) >= ec[i]) 
    # m.ext[:constraints][:ec_ed_con] =  @constraint(m, [(n,i,j) in B_dc_to], ec[j] + ec[i] >= 1.1*ed[n]) 






    # m.ext[:constraints][:pbcvij] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_dc_node[i] -  pb_ph_dc_ac[(n,i,j)] == (a_loss_cv[n]* ec[n] + b_loss_cv[n] * (-pb_ph_dc_ac[(n,i,j)])))
    m.ext[:constraints][:pbcvij] = @constraint(m, [(n,i,j) = B_tf_dc_ac],- pb_dc_node[i] -  pb_ph_dc_ac[(n,i,j)] == (a_loss_cv[n]* ec[n] + b_loss_cv[n] * (1 - 2*(ei[n]))*(i_cv[n]))) # converter losses
    # m.ext[:constraints][:pbcvij] = @constraint(m, [(n,i,j) = B_tf_dc_ac], -pb_dc_node[i] -  pb_ph_dc_ac[(n,i,j)] == 0)
    m.ext[:constraints][:pqbcvij] = @constraint(m, [(n,i,j) = B_tf_dc_ac], pb_ph_dc_ac[(n,i,j)]^2 + qb_ph_dc_ac[(n,i,j)]^2 <= pmax_tf[n]^2 + qmax_tf[n]^2)



    # # DC node power flow contrains 
    m.ext[:constraints][:p_dc_balance] = @constraint(m, [n in N_dc_ac_connected],  sum(pb_dc_node[n_dc] for n_dc in N_dc_index[n]) - sum( pb_dc[(b,i,j)] for (b,i,j) in B_dc_arcs[n]) == 0)


 
    # Branch angle limits
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,i,j) = B_ac_fr], va[i] - va[j] >= angmin[b])
    m.ext[:constraints][:thetaij] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] <= angmax[b])
    m.ext[:constraints][:thetaji] = @constraint(m, [(b,j,i) = B_ac_to], va[j] - va[i] >= angmin[b])

    # m.ext[:constraints][:thetaij_tf] = @constraint(m, [(b,i,j) = B_tf_ac_dc], va[i] - va_tf[j] <= angmax[b])
    # m.ext[:constraints][:thetaji_tf] = @constraint(m, [(b,i,j) = B_tf_ac_dc], va[i] - va_tf[j] >= angmin[b])
    # m.ext[:constraints][:thetaij_tf] = @constraint(m, [(b,j,i) = B_tf_dc_ac], va_tf[j] - va[i] <= angmax[b])
    # m.ext[:constraints][:thetaji_tf] = @constraint(m, [(b,j,i) = B_tf_dc_ac], va_tf[j] - va[i] >= angmin[b])

    # m.ext[:constraints][:thetaij_ph] = @constraint(m, [(b,i,j) = B_tf_ac_dc], va_tf[i] - va_ph[j] <= angmax[b])
    # m.ext[:constraints][:thetaji_ph] = @constraint(m, [(b,i,j) = B_tf_ac_dc], va_tf[i] - va_ph[j] >= angmin[b])
    # m.ext[:constraints][:thetaij_ph] = @constraint(m, [(b,j,i) = B_tf_dc_ac], va_ph[j] - va_tf[i] <= angmax[b])
    # m.ext[:constraints][:thetaji_ph] = @constraint(m, [(b,j,i) = B_tf_dc_ac], va_ph[j] - va_tf[i] >= angmin[b])



    # # Kirchhoff's current law, i.e., nodal power balance

    m.ext[:constraints][:p_balance_dc_connected] = @constraint(m, [i in N_dc_ac_connected], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) -  sum(pb_tf_ac_dc[(j,k,j)] for (j,k) in N_all_dc_connected[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) ) #pb verwijderen en pg
    m.ext[:constraints][:q_balance_dc_connected] = @constraint(m, [i in N_dc_ac_connected], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) - sum(qb_tf_ac_dc[(j,k,j)] for (j,k) in N_all_dc_connected[i])  == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    m.ext[:constraints][:p_balance] = @constraint(m, [i in N_no_dc_connected], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    m.ext[:constraints][:q_balance] = @constraint(m, [i in N_no_dc_connected], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )

    # m.ext[:constraints][:p_balance] = @constraint(m, [i in N], sum(pg[g] for g in G_ac[i]) - sum(pd[l] for l in L_ac[i]) == sum(pb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
    # m.ext[:constraints][:q_balance] = @constraint(m, [i in N], sum(qg[g] for g in G_ac[i]) - sum(qd[l] for l in L_ac[i]) == sum(qb[(b,i,j)] for (b,i,j) in B_arcs[i]) )
 

    m.ext[:constraints][:varef] = @constraint(m, [n_sl in N_sl], va[n_sl] == 0)
    return m 
 
end