from sympy import *
import numpy as np
# from sympy import init_printing #(for notebook)
# # Initialize pretty printing
# init_printing(use_unicode=True)
import numpy as np
from src.IO.xcel import *
from src.utils.config_upd import *


class DC_system_initialization(object):
    def __init__(self,case_name,update_dc_params={}):
        self.bus = read_excel(case_name, sheet_name='Bus')
        self.line = read_excel(case_name, sheet_name='Line')
        self.dcac_interface_dcbus = list(read_excel(case_name, sheet_name='Interface')['DC Bus'])
        self.dcac_interface_acbus = list(read_excel(case_name, sheet_name='Interface')['AC Bus'])
        self.dc_slack_bus = list(read_excel(case_name, sheet_name='Slack')['dc_bus_no'])[0]
        self.ac_side_of_dc_slack = list(read_excel(case_name, sheet_name='Slack')['ac_bus_no'])[0]
        self.dc_pars = read_excel(case_name, sheet_name='DC_params_status')

        print('Default status DC-net params', self.dc_pars)

        self.dc_pars = update_config(self.dc_pars,update_dc_params)

        print('Current status DC-net params', self.dc_pars)

    def preprocess(self):
        self.dc_buses = [list(self.bus['name'])[k]-1 for k in range(len(list(self.bus['name'])))]
        self.branches = [(branch[0]-1,branch[1]-1,branch[2])\
                          for branch in list(self.line[['bus1','bus2','r']].itertuples(index=False, name=None))]
        self.pinj = list(self.bus[['name','Pinj']].itertuples(index=False, name=None))
        self.capacitances = list(self.bus['c'])
        self.kp_dc = list(self.bus['kp'])
        self.ki_dc = list(self.bus['ki'])
        self.Vref = list(self.bus['Vref'])

class DC_network_solution:
    '''
    The class models a DC network for steady-state and dynamic power flow analysis.  
    It supports steady-state DC power flow calculations and dynamic simulations with  
    capacitor-based voltage dynamics and PI control for regulating slack bus voltage.

    Inputs:
    -------
    - buses: list  
        List of integers representing bus indices in the DC network.  
    - branches: list of tuples  
        Each tuple (from_bus, to_bus, resistance) specifies a branch between two buses and its resistance (in ohms).  
    - Pdc: numpy array  
        Initial active power injections (in per-unit) for each bus.  
    - capacitances: numpy array  
        Capacitance values (in per-unit) connected at each DC bus for dynamic modeling.  
    '''

    def __init__(self, buses, branches, Pdc, capacitances):
        self.buses = buses
        self.branches = branches
        self.P = Pdc
        self.capacitances = capacitances
        self.Y_bus = self.get_ybus()
        self.V = None

    def get_ybus(self):
        '''
        Compute the Y-bus matrix based on branch resistances.
        '''
        n = len(self.buses)
        Y_bus = np.zeros((n, n))
        for branch in self.branches:
            i, j, r = branch
            y = 1 / r
            Y_bus[i, i] += y
            Y_bus[j, j] += y
            Y_bus[i, j] -= y
            Y_bus[j, i] -= y
        return Y_bus
    
    def dcpf(self, Udc, Pdc):
        '''
        Perform steady-state DC power flow analysis using Newton-Raphson method.
        Parameters:
        - Udc: Initial guesses for voltages.
        - Pdc: Active power injections.

        Returns:
        - Updated voltages (Udc) and power injections (Pdc).
        '''
        n = len(self.Y_bus)
        f = np.ones(n)
        it = 0
        while max(abs(f)) > 1e-12 and it < 10:
            it += 1
            f = np.matmul(self.Y_bus, Udc)
            for jj in range(n):
                f[jj] -= Pdc[jj] / Udc[jj]

            Y11 = self.Y_bus[0:n-1, 0:n-1]
            I = np.identity(n-1)
            for ii in range(n-1):
                I[ii, ii] = Pdc[ii] / (Udc[ii] ** 2)
            J11 = Y11 + I
            J12 = np.zeros([n-1, 1])
            J1 = np.concatenate((J11, J12), axis=1)

            J21 = np.array([self.Y_bus[n-1, 0:n-1]])
            J22 = np.zeros([1, 1])
            J22[0, 0] = -1 / Udc[n-1]
            J2 = np.concatenate((J21, J22), axis=1)

            J = np.concatenate((J1, J2), axis=0)
            x = np.concatenate((Udc[0:n-1], Pdc[n-1:n]), axis=0)
            x_new = x - np.matmul(np.linalg.inv(J), f)

            for ii in range(n-1):
                Udc[ii] = x_new[ii]
            Pdc[n-1] = x_new[n-1]

        self.V = Udc
        self.Pdc = Pdc
        return Udc, Pdc
    

class DC_network_equations(object):
    def __init__(self,dc_system,dc_network):
        self.dc_buses = dc_system.dc_buses
        self.Vref = dc_system.Vref
        self.dcac_interface_dcbus = dc_system.dcac_interface_dcbus
        self.dcac_interface_acbus = dc_system.dcac_interface_acbus
        self.dc_slack_bus = dc_system.dc_slack_bus
        # self.r_dc = [var(f"R_dc{self.dc_buses[k]+1}") for k in range(len(self.dc_buses))]
        self.Ybus = dc_network.get_ybus()
        self.capacitances = dc_system.capacitances
        self.kp_dc = dc_system.kp_dc
        self.ki_dc = dc_system.ki_dc
        self.Vref = dc_system.Vref
        self.dc_pars = dc_system.dc_pars

    def dc_vars(self):
        self.v_dc = [var(f"V_dc{self.dc_buses[k]+1}") for k in range(len(self.dc_buses))]
        self.i_dc = [var(f"I_dc{self.dc_buses[k]+1}") for k in range(len(self.dc_buses))]
        self.psi = var('psi')

    def dc_params(self):
        self.c_dc = [
            var(f"C_dc{self.dc_buses[k]+1}") if self.dc_pars['c'][0] in ['sym', 'user']
            else self.capacitances[k] for k in range(len(self.dc_buses))
            ]
        self.kp_dc = [
            var(f"K_pdc{self.dc_buses[k]+1}") if self.dc_pars['kp'][0] in ['sym', 'user']
            else self.kp_dc[k] for k in range(len(self.dc_buses))
            ]

        self.ki_dc = [
            var(f"K_idc{self.dc_buses[k]+1}") if self.dc_pars['ki'][0] in ['sym', 'user']
            else self.ki_dc[k] for k in range(len(self.dc_buses))
            ]
        
        self.Vref = [
            var(f"V_ref{self.dc_buses[k]+1}") if self.dc_pars['Vref'][0] in ['sym', 'user']
            else self.Vref[k] for k in range(len(self.dc_buses))
            ]

    def dc_states(self):
        self.xdc = [self.v_dc[k] for k in range(len(self.dc_buses))] 
        self.xdc.append(self.psi)

    def dc_dynamics(self):
        numDCBus = len(self.dc_buses)
        self.dot_xdc = [None] * (numDCBus+1)
        for k in range(numDCBus):
            self.dot_Vdc = (1/self.c_dc[k])*self.i_dc[k]
            for j in range(numDCBus):
                self.dot_Vdc += (1/self.c_dc[k])*(-(self.v_dc[j]-self.v_dc[k])*self.Ybus[k,j]) 
            self.dot_xdc[k] = self.dot_Vdc

        self.dot_xdc[self.dc_slack_bus-1] = self.dot_xdc[self.dc_slack_bus-1] + \
                                        (1/self.c_dc[self.dc_slack_bus-1])*((-self.ki_dc[k])*self.psi \
                                                                               + (-self.kp_dc[k])*(self.v_dc[self.dc_slack_bus-1]- self.Vref[k]))
        self.dot_xdc[numDCBus] = self.v_dc[self.dc_slack_bus-1] - self.Vref[k]