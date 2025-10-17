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
from src.utils.py_2_mat import *

# Define base path relative to the project directory
case_path = os.path.join(parent_dir, "cases\\Case_5bus_v1")

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
    # "xq": 'num',   
    # "ra": 'num',
    "Sl": 'num',
}
# DC params: c	kp	ki	Vref
user_updates_dc_params = {
    "c": 'num'
}

dc_case_file,ac_base_file,ac_pf,\
    ac_base_sym, eq_data_file, eq_data_file_extended = case_files(case_path)

dc_system,dc_network,dc_syms = dc_wrapper(dc_case_file,user_updates_dc_params)

system,bus,model,alg = ac_wrapper(ac_base_sym,user_updates_gen_params,user_updates_gfl_params)

f,g,x,y = ac_dc_connection(system,bus,model,alg,dc_system,dc_network,dc_syms)

ac_pf_update(ac_base_file,ac_pf,'PQ',dc_system,dc_network,system)

base_equ_comp(case_path,ac_pf,eq_data_file,system,alg)

ac_ext_net_sol = extended_equ_comp(eq_data_file,eq_data_file_extended,system,alg)

equ_value = equ_val_assignment(dc_network,dc_system,dc_syms,ac_ext_net_sol,system,bus,model,alg)

if max(abs(g.subs(equ_value))) < 1e-5 and max(abs(f.subs(equ_value))) < 1e-5:
    print("Correct Equilibrium")
else:
    print("Incorrect Equilibrium")

# compute jacobians
fx_acdc = f.jacobian(x)
fy_acdc = f.jacobian(y)
gx_acdc = g.jacobian(x)
gy_acdc = g.jacobian(y)

fx_acdc.simplify()
fy_acdc.simplify()
gx_acdc.simplify()
gy_acdc.simplify()

fx_sub = fx_acdc.subs(equ_value)
fy_sub = fy_acdc.subs(equ_value)
gx_sub = gx_acdc.subs(equ_value)
gy_sub = gy_acdc.subs(equ_value)

inv_gy = inv(np.array(gy_sub,dtype=float))
if np.linalg.det(np.array(gy_sub,dtype=float)) != 0:
    print("Correct Gy")

A = fx_sub - fy_sub@inv_gy@gx_sub

print(A)

# Getting the details of user-defined parameters:
# - Check which variables are defined as user-defined parameters.
# - Assign values for those variables (please map the correct symbolic notation)
    
dataframes = {'SG': system.gen_pars, 'GFL': system.gfl_pars, 'GFM': system.gfm_pars, 'DC': dc_system.dc_pars}
results = {}
for name, df in dataframes.items():  
    user_columns = []  
    
    for column in df.columns:
        if df[column].isin(['user']).any():  
            user_columns.append(column)
    if user_columns:
        results[name] = user_columns

print("User-defined parameters")
for dataframe, columns in results.items():
    print(f"{dataframe}: {columns}")

# update parameter values with user-defined ones
mp_hat = 0.01
values_ud = {}
for k in range(len(system.gfl)):
    values_ud[model.k_droop[k]] = 1/mp_hat 
    values_ud[model.kiPLL[k]] = 1.7 #5  
    values_ud[model.kpPLL[k]] = 0.1 #0.2 

# update some parameter values with user-defined ones
A = A.subs(values_ud)

A1 = np.array(A,dtype=float)

eigenvalues,_ = np.linalg.eig(A1) 

print('Eigenvalues:')
for k in range(A1.shape[0]):
    print(np.round(eigenvalues[k],4))

# ### Saving Parametric A matrix for Stability characterization:
# - Save the Parametric A matrix as a MATLAB script.
# - Exporting parameters and their values for verification.

variables = model.kpPLL + model.kiPLL + model.k_droop  # Combined symbolic variables (modify this line according to thr user defined var setting)
var_vals = values_ud

# Define the directory to check/create
matlab_dir = os.path.join(parent_dir, "examples\\matlab_files")

# Check if the folder exists; if not, create the folder
if not os.path.exists(matlab_dir):
    os.makedirs(matlab_dir)
    print(f"Folder '{matlab_dir}' created.")
else:
    print(f"Folder '{matlab_dir}' already exists.")

# Combine symbolic MATLAB expressions with the numeric matrix

var_name = 'A'  # Variable name for the matrix
output_matrix = A  # A matrix with high precision
line_limit = 80  # Set line limit for MATLAB code formatting

# Generate MATLAB code
matlab_code = (f"clear all;\n" +
    f"clc;\n\n" +
    matlab_syms(variables) +  # Generate symbolic definitions
    matrix_to_matlab(output_matrix, var_name, line_limit=80) + "\n\n" +  # Matrix definition (high precision)
    matlab_value_assignments(var_vals) + "\n" +  # Variable assignments
    f"disp('{var_name} with substituted values:');\ndisp({var_name}_sub);\n" +
    f"disp('Numeric result:');\ndisp({var_name}_numeric);\n" +
    f"disp('Eigenvalues:');\ndisp(eig({var_name}_numeric));"
)

# Output file location
out_file = os.path.join(matlab_dir, 'symbolic_matrix_ACDC_'+os.path.split(case_path)[1]+'.m')

# Write the MATLAB code to the .m file
with open(out_file, 'w') as f1:
    f1.write(matlab_code)

print(f"MATLAB code saved to '{out_file}' successfully.")