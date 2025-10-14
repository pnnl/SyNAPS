from src.components.ac_bus import *
from src.components.ac_dynamic_models import *
from src.components.ac_network import *
from src.dae.equation_form import *
from src.utils.cal_send_v import *
from src.IO.xcel import *
from src.acdc.inter_connect import *
from src.dae.equation_form import *
from src.core.sys import *
from src.components.dc_network import *
from src.utils.py_2_mat import *

import numpy as np
import shutil
import os
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)


def dc_wrapper(dc_case_name,update_dc_params={}):
    '''
    Based on the dc_case_name, this function combines 3 dc classes of dc_network.py
    '''
    dc_system = DC_system_initialization(dc_case_name,update_dc_params)
    dc_system.preprocess()

    dc_buses = dc_system.dc_buses
    branches = dc_system.branches
    Pinj = dc_system.pinj
    capacitances = dc_system.capacitances

    # Initialize the DC network
    dc_network = DC_network_solution(dc_buses, branches, Pinj, capacitances)
    Ybus_DC = dc_network.get_ybus()
    # 1. Setup initial conditions
    Udc = np.ones(len(dc_buses))  # Initial voltages (p.u.)
    Pdc = np.array([p for _, p in Pinj])  # Power injections/absorptions (MW)

    # 2. Calculate DC Power Flow using dcpf
    Udc_final, Pdc_final = dc_network.dcpf(Udc, Pdc)

    # 3. Display Results
    print("DC bus Voltages (p.u.): \n", Udc_final)
    print("DC bus Powers (1e3 MW): \n", Pdc_final)

    dc_syms = DC_network_equations(dc_system,dc_network)
    dc_syms.dc_vars()
    dc_syms.dc_params()
    dc_syms.dc_states()
    dc_syms.dc_dynamics()

    return dc_system,dc_network,dc_syms

def ac_pf_update(ac_case_name_base,ac_case_name_pf,sheet_name,dc_system,dc_network,system):
    '''
    1. Based on the updated DC injections update the **PQ** sheet of the **xxxxAC.xlsx**, and create a new file **xxxxAC_with_dc_mod.xlsx**
    2. (if separate GFL exists apart from DC-AC converter) Based on the updated GFL Pref/Qref update the **PQ** sheet of the **xxxxAC.xlsx**, 
    and create a new file **xxxxAC_with_dc_mod.xlsx** 
    3. Based on the updated GFM Pset/Qset update the **PQ** sheet of the **xxxxAC.xlsx**, 
    and create a new file **xxxxAC_with_dc_mod.xlsx** 
    
    '''
    if not os.path.exists(ac_case_name_pf):
        # print(f"{target_file} does not exist.")

        # Check if the source file exists
        if os.path.exists(ac_case_name_base):
            shutil.copy(ac_case_name_base, ac_case_name_pf)
        else:
            print(f"Source file {ac_case_name_base} does not exist. Cannot create {ac_case_name_pf}.")
    else:
        print(f"{ac_case_name_pf} already exists.")
    # import AC system data and update the DC/AC PCC injection and save as same file name
    Pdc_final = dc_network.Pdc
    df = excel_sheet_extraction(ac_case_name_base,sheet_name) 
    for j in range(len(dc_system.dcac_interface_acbus)):
        old_val = df.at[df[df['bus'] == dc_system.dcac_interface_acbus[j]].index[0], 'p0']
        df.at[df[df['bus'] == dc_system.dcac_interface_acbus[j]].index[0], 'p0']\
                = old_val - abs(Pdc_final[dc_system.dcac_interface_dcbus[j] - 1])

    # import AC system data and update the GFL (other than DC/AC) PCC injection and save as same file name
    for j in range(len(system.gfl)):
        if system.gfl['bus'][j] not in dc_system.dcac_interface_acbus:
            old_val_p = df.at[df[df['bus'] == system.gfl['bus'][j]].index[0], 'p0']
            old_val_q = df.at[df[df['bus'] == system.gfl['bus'][j]].index[0], 'q0']
            df.at[df[df['bus'] == system.gfl['bus'][j]].index[0], 'p0']\
                    = old_val_p - abs(system.gfl['Pref'][j]) 
            df.at[df[df['bus'] == system.gfl['bus'][j]].index[0], 'q0']\
                    = old_val_q - abs(system.gfl['Qref'][j]) 
            
    # import AC system data and update the GFM (other than DC/AC) PCC injection and save as same file name
    for j in range(len(system.gfm)):
        old_val_p = df.at[df[df['bus'] == system.gfm['bus'][j]].index[0], 'p0']
        old_val_q = df.at[df[df['bus'] == system.gfm['bus'][j]].index[0], 'q0']
        df.at[df[df['bus'] == system.gfm['bus'][j]].index[0], 'p0']\
                = old_val_p - abs(system.gfm['Pset'][j]) 
        df.at[df[df['bus'] == system.gfm['bus'][j]].index[0], 'q0']\
                = old_val_q - abs(system.gfm['Qset'][j]) 
    excel_sheet_replace(ac_case_name_pf,sheet_name,df)


def ac_wrapper(ac_case_name,updates_gen={},updates_gfl={},updates_gfm={}):
    '''
    Based on the ac_case_name, this function combines all ac classes.
    '''
    system = get_ac_model(ac_case_name,updates_gen,updates_gfl,updates_gfm)
    bus = Bus(system)
    model = AC_Dynamic_models(system,bus)

    model.SG_vars()
    model.SG_params()
    model.SG_dynamics()
    model.excitor_dynamics()
    model.GFL_vars()
    model.GFL_params()
    model.GFL_dynamics()
    model.GFM_vars()
    model.GFM_params()
    model.GFM_dynamics()
    model.gen_var_classification()
    model.conv_var_classification()

    alg = network_equations(system,bus,model)

    alg.branch_extension()
    alg.bus_extension()
    alg.alg_var_extension()
    alg.Y_bus()
    alg.injections()
    alg.include_load()
    alg.include_gen()

    return system,bus,model,alg

def ac_dc_connection(system,bus,model,alg,dc_system,dc_network,dc_syms):
    '''
    This function includes the dc injections in the symbolic form of algebraic equations.
    '''
    alg.P_inj,alg.Q_inj,dc_alg_eq,dc_alg_var = DC_injection_inclusion_v1(system,bus,model,alg,dc_system,dc_syms)
    alg.update_gen_conv_dynamics()
    # write DAE for ac system
    dae_ac = DAE(system,bus,model,alg)

    dae_ac.x_var()
    dae_ac.y_var()
    dae_ac.f_func()
    dae_ac.g_func()
    dae_ac.update_y_var()

    # combine ac and dc DAEs
    f,g,x,y = acdc_DAE_v1(dae_ac,dc_syms,dc_alg_eq,dc_alg_var)

    return f,g,x,y

def ac_connection(system,bus,model,alg):
    alg.update_gen_conv_dynamics()
    dae_ac = DAE(system,bus,model,alg)

    dae_ac.x_var()
    dae_ac.y_var()
    dae_ac.f_func()
    dae_ac.g_func()
    dae_ac.update_y_var()

    # combine ac and dc DAEs
    f,g,x,y = ac_DAE(dae_ac)
    return f,g,x,y