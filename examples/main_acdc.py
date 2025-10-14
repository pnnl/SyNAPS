import os 
import sys
from sympy import *
import pandas as pd
from scipy.integrate import odeint
from sympy import symbols, init_printing
import matplotlib.pyplot as plt
import cmath
import numpy as np

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(parent_dir)

from src.IO.xcel import *
from src.core.wrapper import *
from src.IO.get_file_names import *
from src.utils.equ_comp import *
from numpy.linalg import inv
from src.utils.get_linearization_equ import *

# Define base path relative to the project directory
case_path = os.path.join(parent_dir, "cases\\Case_5bus")

# (1) 'sym': Symbolic - parameters are treated as symbolic variables for analytical computations or symbolic modeling.
# (2) 'num': Numeric - parameters are assigned specific numeric values for direct numerical simulations.
# (3) 'user': User-defined - parameters are provided by the user as custom inputs for tailored simulations or system-specific settings.

# SG params: M	D	xd	xq	xd1	xq1	Td10	Tq10	ra	Tm	Rgov	Te	Ke	Vref
user_updates_gen_params = {
    "M": 'num',   
    "D": 'num'
}

# GFL params: KpPLL	KiPLL	kdroop	xq	ra	kxl	Kp	Ki	Sl
user_updates_gfl_params = {
    "xq": 'sym',   
    "ra": 'sym'
}

# DC params: c	kp	ki	Vref
user_updates_dc_params = {
    "c": 'sym'
}

dc_case_file,ac_base_file,ac_pf,\
    ac_base_sym, eq_data_file, eq_data_file_extended = case_files(case_path)

dc_system,dc_network,dc_syms = dc_wrapper(dc_case_file,user_updates_dc_params)

system,bus,model,alg = ac_wrapper(ac_base_sym,user_updates_gen_params,user_updates_gfl_params)

f,g,x,y = ac_dc_connection(system,bus,model,alg,dc_system,dc_network,dc_syms)

ac_pf_update(ac_base_file,ac_pf,'PQ',dc_system,dc_network,system)

# compute jacobians
fx_acdc = f.jacobian(x)
fy_acdc = f.jacobian(y)
gx_acdc = g.jacobian(x)
gy_acdc = g.jacobian(y)

fx_acdc.simplify()
fy_acdc.simplify()
gx_acdc.simplify()
gy_acdc.simplify()


base_equ_comp(case_path,ac_pf,eq_data_file,system,alg)

ac_ext_net_sol = extended_equ_comp(eq_data_file,eq_data_file_extended,system,alg)

equ_value = equ_val_assignment(dc_network,dc_system,dc_syms,ac_ext_net_sol,system,bus,model,alg)


fx_sub = fx_acdc.subs(equ_value)
fy_sub = fy_acdc.subs(equ_value)
gx_sub = gx_acdc.subs(equ_value)
gy_sub = gy_acdc.subs(equ_value)

inv_gy = inv(np.array(gy_sub,dtype=float))
if np.linalg.det(inv_gy) != 0:
    print("Correct Gy")

A = fx_sub - fy_sub@inv_gy@gx_sub

print(A)

# update some parameter values with user-defined ones
mp_hat = 0.010
values_ud = {}
for k in range(len(system.gfl)):
    values_ud[model.gfl_cap[k]] = 1.0
    values_ud[model.k_droop[k]] = 1/mp_hat 
    values_ud[model.kiPLL[k]] = 1.7 #5  
    values_ud[model.kpPLL[k]] = 0.1 #0.2 


A_matrix1 = A.subs(values_ud)
A1 = np.array(A_matrix1,dtype=float)


eigenvalues,_ = np.linalg.eig(A1) 

print('Eigenvalues:')
for k in range(A_matrix1.shape[0]):
    print(np.round(eigenvalues[k],4))