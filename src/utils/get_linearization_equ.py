import numpy as np
def equ_val_assignment(dc_network,dc_system,dc_syms,ac_ext_net_sol,system,bus,model,alg):
    '''
    This function collects ac-dc network details, and ac-dc power flow solutions, and assigns the equilibrium values 
    for each varaibales.
    '''

    # collect dc solutions
    Vdc_eq = dc_network.V
    psi_eq = 0
    Iinj = dc_network.Pdc/dc_network.V
    Idc_eq = Iinj

    # assign variable values
    values = {}

    # assign DC variable values
    values[dc_syms.psi] = psi_eq 
    for k in range(len(dc_system.dc_buses)):
        values[dc_syms.v_dc[k]] = Vdc_eq[k]
        values[dc_syms.i_dc[k]] = Idc_eq[k]

    values.update({dc_syms.c_dc[k]: dc_system.capacitances[k] for k in range(len(dc_system.dc_buses)) if dc_system.dc_pars['c'][0] == 'sym'})
    values.update({dc_syms.kp_dc[k]: dc_system.kp_dc[k] for k in range(len(dc_system.dc_buses)) if dc_system.dc_pars['kp'][0] == 'sym'})
    values.update({dc_syms.ki_dc[k]: dc_system.ki_dc[k] for k in range(len(dc_system.dc_buses)) if dc_system.dc_pars['ki'][0] == 'sym'})
    values.update({dc_syms.Vref[k]: dc_system.Vref[k] for k in range(len(dc_system.dc_buses)) if dc_system.dc_pars['Vref'][0] == 'sym'})


    # assign AC variable values
    for k in range(len(alg.bus_extended)):
        values[bus.a[k]] = ac_ext_net_sol.Aext[k]
        values[bus.v[k]] = ac_ext_net_sol.Vext[k]

    for k in range(len(system.gfl)):
        r = np.where(alg.extended_resource_name.astype(str)==list(system.gfl['name'])[k])[0]
        new_bus_idx = alg.extended_resource_name[r,3][0] - 1
        old_bus_idx = alg.extended_resource_name[r,2][0] - 1
        values[model.a_gfl[k]] = ac_ext_net_sol.Aext[old_bus_idx]
        values[model.w_gfl[k]] = 1
        values[model.gfl_id[k]] = ac_ext_net_sol.Idext[new_bus_idx]
        values[model.gfl_iq[k]] = ac_ext_net_sol.Iqext[new_bus_idx]
        values[model.Q_ref[k]] = ac_ext_net_sol.Qg[new_bus_idx]

        if (old_bus_idx+1) not in dc_system.dcac_interface_acbus:
            values[model.P_ref[k]] = ac_ext_net_sol.Pg[new_bus_idx]


    for k in range(len(system.gen)):
        if model.gen_order[k] == 4:
            r = np.where(alg.extended_resource_name.astype(str)==list(system.gen['name'])[k])[0]
            new_bus_idx = alg.extended_resource_name[r,3][0] - 1
            values[model.id[k]] = ac_ext_net_sol.Idext[new_bus_idx]
            values[model.iq[k]] = ac_ext_net_sol.Iqext[new_bus_idx]
        
    for k in range(len(system.bus)):
        values[bus.pl[k]] = ac_ext_net_sol.Pl[k]
        values[bus.ql[k]] = ac_ext_net_sol.Ql[k]

    for k in range(len(system.gen)):
        r = np.where(alg.extended_resource_name.astype(str)==list(system.gen['name'])[k])[0]
        old_bus_idx = alg.extended_resource_name[r,2][0] - 1
        new_bus_idx = alg.extended_resource_name[r,3][0] - 1
        values[model.Q_g[k]] = ac_ext_net_sol.Qg[new_bus_idx]  
        values[model.P_g[k]] = ac_ext_net_sol.Pg[new_bus_idx]
        values[model.w_gen[k]] = 1
    for k in model.order_4_gens:
        values[model.P_mech[k]] = ac_ext_net_sol.Pg[new_bus_idx]
        # values[model.P_mech[k]] = ac_ext_net_sol.Pg[new_bus_idx]

    values.update({model.M_gen[k]: list(system.gen['M'])[k] for k in range(len(system.gen)) if model.gen_pars['M'][0]== 'sym'})
    values.update({model.D_gen[k]: list(system.gen['D'])[k] for k in range(len(system.gen)) if model.gen_pars['D'][0]== 'sym'})
    values.update({model.xd[k]: list(system.gen['xd'])[k] for k in range(len(system.gen)) if model.gen_pars['xd'][0]== 'sym'})
    values.update({model.xq[k]: list(system.gen['xq'])[k] for k in range(len(system.gen)) if model.gen_pars['xq'][0]== 'sym'})
    values.update({model.xd1[k]: list(system.gen['xd1'])[k] for k in range(len(system.gen)) if model.gen_pars['xd1'][0]== 'sym'})
    values.update({model.xq1[k]: list(system.gen['xq1'])[k] for k in range(len(system.gen)) if model.gen_pars['xq1'][0]== 'sym'})
    values.update({model.Td01[k]: list(system.gen['Td10'])[k] for k in range(len(system.gen)) if model.gen_pars['Td10'][0]== 'sym'})
    values.update({model.Tq01[k]: list(system.gen['Tq10'])[k] for k in range(len(system.gen)) if model.gen_pars['Tq10'][0]== 'sym'})
    values.update({model.Ra[k]: list(system.gen['ra'])[k] for k in range(len(system.gen)) if model.gen_pars['ra'][0]== 'sym'})

    values.update({model.T_mech[k]: list(system.tgov['Tm'])[k] for k in range(len(system.gen)) if model.gen_pars['Tm'][0]== 'sym'})
    values.update({model.R_gov[k]: list(system.tgov['Rgov'])[k] for k in range(len(system.gen)) if model.gen_pars['Rgov'][0]== 'sym'})
    
    values.update({model.Te[k]: list(system.exc['Te'])[k] for k in range(len(system.gen)) if model.gen_pars['Te'][0]== 'sym'})
    values.update({model.Ke[k]: list(system.exc['Ke'])[k] for k in range(len(system.gen)) if model.gen_pars['Ke'][0]== 'sym'})
    values.update({model.vref[k]: list(system.exc['Vref'])[k] for k in range(len(system.gen)) if model.gen_pars['Vref'][0]== 'sym'})

    values.update({model.k_droop[k]: list(system.gfl['kdroop'])[k] for k in range(len(system.gfl)) if model.gfl_pars['kdroop'][0]== 'sym'})
    values.update({model.kpPLL[k]: list(system.gfl['KpPLL'])[k] for k in range(len(system.gfl)) if model.gfl_pars['KpPLL'][0]== 'sym'})
    values.update({model.kiPLL[k]: list(system.gfl['KiPLL'])[k] for k in range(len(system.gfl)) if model.gfl_pars['KiPLL'][0]== 'sym'})
    values.update({model.kpc[k]: list(system.gfl['Kp'])[k] for k in range(len(system.gfl)) if model.gfl_pars['Kp'][0]== 'sym'})
    values.update({model.kic[k]: list(system.gfl['Ki'])[k] for k in range(len(system.gfl)) if model.gfl_pars['Ki'][0]== 'sym'})
    values.update({model.gfl_cap[k]: list(system.gfl['Sl'])[k] for k in range(len(system.gfl)) if model.gfl_pars['Sl'][0]== 'sym'})
    values.update({model.Xl[k]: list(system.gfl['xq'])[k] for k in range(len(system.gfl)) if model.gfl_pars['xq'][0]== 'sym'})
    values.update({model.K[k]: list(system.gfl['kxl'])[k] for k in range(len(system.gfl)) if model.gfl_pars['kxl'][0]== 'sym'})
    values.update({model.gfl_Ra[k]: list(system.gfl['ra'])[k] for k in range(len(system.gfl)) if model.gfl_pars['ra'][0]== 'sym'})

    values.update({model.mp[k]: list(system.gfm['mp'])[k] for k in range(len(system.gfm)) if model.gfm_pars['mp'][0]== 'sym'})
    values.update({model.mq[k]: list(system.gfm['mq'])[k] for k in range(len(system.gfm)) if model.gfm_pars['mq'][0]== 'sym'})
    values.update({model.cap[k]: list(system.gfm['Sl'])[k] for k in range(len(system.gfm)) if model.gfm_pars['Sl'][0]== 'sym'})
    values.update({model.kpv[k]: list(system.gfm['kpv'])[k] for k in range(len(system.gfm)) if model.gfm_pars['kpv'][0]== 'sym'})
    values.update({model.kiv[k]: list(system.gfm['kiv'])[k] for k in range(len(system.gfm)) if model.gfm_pars['kiv'][0]== 'sym'})
    values.update({model.tau[k]: list(system.gfm['tau'])[k] for k in range(len(system.gfm)) if model.gfm_pars['tau'][0]== 'sym'})


    for k in range(len(system.gfm)):
        r = np.where(alg.extended_resource_name.astype(str)==list(system.gfm['name'])[k])[0]
        old_bus_idx = alg.extended_resource_name[r,2][0] - 1
        new_bus_idx = alg.extended_resource_name[r,3][0] - 1
        values[model.Q_set[k]] = ac_ext_net_sol.Qg[new_bus_idx]  
        values[model.P_set[k]] = ac_ext_net_sol.Pg[new_bus_idx]
        values[model.V_set[k]] = ac_ext_net_sol.Vext[old_bus_idx]
        values[model.a_gfm[k]] = ac_ext_net_sol.Aext[new_bus_idx]
        values[model.w_gfm[k]] = 1
        values[model.v_conv_err[k]] = 0

    return values