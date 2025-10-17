from sympy import *
import numpy as np

def DC_injection_inclusion_v1(system0,bus0,model0,alg0,dc_system0,dc_syms0):
    '''
    This function is to locate the dc-ac interface buses and update ac-side symbolic algebraic equations to 
    include the dc injections with droop parameters.
    '''
    for k in range(len(system0.gfl)):
        r = np.where(alg0.extended_resource_name.astype(str)==list(system0.gfl['name'])[k])[0]
        # old_bus_idx = alg0.extended_resource_name[r,2][0] - 1
        new_bus_idx = alg0.extended_resource_name[r,3][0] - 1

        old_bus_no = alg0.extended_resource_name[r,2][0]
        if old_bus_no in dc_system0.dcac_interface_acbus:
            model0.P_ref[k] = -dc_syms0.v_dc[dc_system0.dcac_interface_dcbus[k]-1]\
                        *dc_syms0.i_dc[dc_system0.dcac_interface_dcbus[k]-1]

        alg0.Q_inj[new_bus_idx] = alg0.Q_inj[new_bus_idx] - \
            model0.Q_ref[k]*cos(bus0.a[model0.gfl_nodes[k]-1] - model0.a_gfl[k])

        alg0.P_inj[new_bus_idx] = alg0.P_inj[new_bus_idx] - \
            (model0.P_ref[k]*cos(bus0.a[model0.gfl_nodes[k]-1] - model0.a_gfl[k])\
                    + (model0.k_droop[k]*model0.gfl_cap[k])*(system0.wb/system0.wb-model0.w_gfl[k])\
                        *cos(bus0.a[model0.gfl_nodes[k]-1] - model0.a_gfl[k]))             

    dcalgeq = []
    dcalgvar = []

    return alg0.P_inj, alg0.Q_inj, dcalgeq, dcalgvar


def acdc_DAE_v1(dae0,dc_syms0,dc_alg_eq,dc_alg_var):
    '''
    Forming ac-dc DAE with ac-dc f, g, x, y
    '''
    f = Matrix(list(dae0.f) + dc_syms0.dot_xdc)

    g = Matrix(list(dae0.g) + dc_alg_eq)

    x = Matrix(list(dae0.x) + dc_syms0.xdc)

    y = Matrix(list(dae0.y) + dc_alg_var)

    return f,g,x,y


def ac_DAE(dae0):
    '''
    Forming ac DAE with ac-only f, g, x, y
    '''
    
    f = dae0.f

    g = dae0.g

    x = dae0.x

    y = dae0.y

    return f,g,x,y




