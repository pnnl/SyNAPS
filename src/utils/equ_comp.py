import numpy as np
import os
from types import SimpleNamespace

from src.utils.cal_send_v import *
from src.utils.andes_wrapper import *
from src.IO.xcel import *
from src.IO.file_move import *


def base_equ_comp(base_path,pf_file_name,target_file_name,system,alg):  
    '''
    This function solves the power flow with the case file including dc injections, and saves the result as 
    .xlsx file.

    Inputs:
    -------
    - base_path: case file directory.
    - pf_file_name: dc_injection updated ac test case file.
    - target_file_name: .xlsx file where the power flow results will be saved.
    - system: created from core/model.py
    - alg: created from componenets/ac_network.py
    '''
    new_loc = os.path.join(base_path, target_file_name)

    Gval, Bval = alg.Y_bus_base_network()

    ss = PF_sol(base_path, pf_file_name)

    Veq = ss.Bus.v.v
    Aeq = ss.Bus.a.v
    P_inj = []
    Q_inj = []
    for k in range(len(system.bus)):
        p = 0
        q = 0
        for j in range(len(system.bus)):
            p += Veq[k]*Veq[j]*(Gval[k,j]*np.cos(Aeq[k]-Aeq[j])\
                                                + Bval[k,j]*np.sin(Aeq[k]-Aeq[j]))
            q += Veq[k]*Veq[j]*(Gval[k,j]*np.sin(Aeq[k]-Aeq[j])\
                                                - Bval[k,j]*np.cos(Aeq[k]-Aeq[j]))            
        P_inj.append(p)
        Q_inj.append(q)
    data = {}
    data['Bus No'] = list(system.bus['idx'])
    data['Pinj'] = P_inj
    data['Qinj'] = Q_inj
    data['Pload'] = list(np.zeros(len(list(system.bus['idx'])),))
    data['Qload'] = list(np.zeros(len(list(system.bus['idx'])),))
    for k in range(len(list(system.pq['bus']))):
        data['Pload'][list(system.pq['bus'])[k]-1] = list(system.pq['p0'])[k]
        data['Qload'][list(system.pq['bus'])[k]-1] = list(system.pq['q0'])[k]

    data['Pg'] = list(np.zeros(len(list(system.bus['idx'])),))
    data['Qg'] = list(np.zeros(len(list(system.bus['idx'])),))
    for k in range(len(list(system.bus['idx']))):
        if data['Pinj'][k] > 0:
            data['Pg'][k] = data['Pinj'][k] + data['Pload'][k]
        else:
            data['Pg'][k] = data['Pinj'][k] + data['Pload'][k]

        if data['Qinj'][k] > 0:
            data['Qg'][k] = data['Qinj'][k] + data['Qload'][k]
        else:
            data['Qg'][k] = data['Qinj'][k] + data['Qload'][k]

    data['V'] = Veq
    data['A'] = Aeq

    new_df = pd.DataFrame(data)
    write_excel(new_loc,new_df)

def extended_equ_comp(pf_saved_data,pf_data_extended,system,alg): 
    '''
    This function calculates the node variables for extended network (including generating resources internal nodes).

    Inputs:
    -------
    - pf_saved_data: .xlsx file where the power flow results are saved using base_equ_comp() function.
    - pf_data_extended: .xlsx file where extended system results will be saved.
    - system: created from core/model.py.
    - alg: created from componenets/ac_network.py.

    Outputs:
    -------
    - returns the bus vol, angle, gen, load for the extended network.  
    '''

    df = read_excel(pf_saved_data,sheet_name=None)
    Veq = df['V']
    Aeq = df['A']

    Vext = list(np.zeros(len(alg.bus_extended))) 
    Aext = list(np.zeros(len(alg.bus_extended))) 

    Idext = list(np.zeros(len(alg.bus_extended)))
    Iqext = list(np.zeros(len(alg.bus_extended)))

    Pg = list(np.zeros(len(alg.bus_extended)))
    Qg = list(np.zeros(len(alg.bus_extended)))

    Vext[:len(list(system.bus['idx']))] = Veq
    Aext[:len(list(system.bus['idx']))] = Aeq

    Pl = list(df['Pload'])
    Ql = list(df['Qload'])

    for k in range(len(alg.branch_info_extended[len(system.line):,:])):
        orig_bus_idx = (alg.branch_info_extended[len(system.line):,:][k,0]).astype(int) - 1
        new_bus_idx = (alg.branch_info_extended[len(system.line):,:][k,1]).astype(int) - 1
        ra = alg.branch_info_extended[len(system.line):,:][k,2] 
        xq = alg.branch_info_extended[len(system.line):,:][k,3]
        V_r = cmath.rect(Veq[orig_bus_idx], Aeq[orig_bus_idx])
        P_r = list(df['Pg'])[list(df['Bus No'])[orig_bus_idx]-1]*alg.extended_resource_name[k,1]
        Q_r = list(df['Qg'])[list(df['Bus No'])[orig_bus_idx]-1]*alg.extended_resource_name[k,1]
        Vext[new_bus_idx],Aext[new_bus_idx],Idext[new_bus_idx],Iqext[new_bus_idx],\
            Pg[new_bus_idx],Qg[new_bus_idx] = \
                calculate_sending_end_voltage(V_r, P_r, Q_r, ra, xq)

        display_data = {}
        display_data['Bus_no'] = alg.bus_extended
        display_data['V_ext'] = np.array(Vext)
        display_data['A_ext'] = np.array(Aext)
        display_data['Idext'] = np.array(Idext)
        display_data['Iqext'] = np.array(Iqext)
        display_data['Pg'] = np.array(Pg)
        display_data['Qg'] = np.array(Qg)
        display_df = pd.DataFrame(display_data)

        write_excel(pf_data_extended,display_df)

        # Create the namespace
        ac_net_sol = SimpleNamespace(
            Vext=Vext,
            Aext=Aext,
            Idext=Idext,
            Iqext=Iqext,
            Pg=Pg,
            Qg=Qg,
            Pl=Pl,
            Ql=Ql
        )
    return ac_net_sol
    
    