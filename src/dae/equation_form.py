from sympy import *
import pandas as pd
from sympy import init_printing

# Initialize pretty printing
init_printing(use_unicode=True)

import andes
from andes.utils.paths import get_case
andes.config_logger()

pd.options.display.max_columns = None
pd.options.display.max_rows = None

class DAE(object):
    
    """Class for Differential Algebraic Equations (dae)"""

    def __init__(self, system, bus, model, network_equations):
        self.system = system
        self.bus = bus
        self.model = model
        self.network_equations = network_equations

    def x_var(self):
        self.x = Matrix(self.model.gen_x + self.model.conv_x)
    def y_var(self):
        self.y = self.bus.a + self.bus.v
        if len(self.model.order_4_gens) != 0:
            self.y = self.y + [self.model.iq[k] for k in range(len(self.model.gen_idx))\
                                                              if self.model.gen_order[k] == 4]
            self.y = self.y + [self.model.id[k] for k in range(len(self.model.gen_idx))\
                                                              if self.model.gen_order[k] == 4]
        self.y = Matrix(self.y)


    def f_func(self):
        self.conv_f = []
        self.gen_f = self.model.dot_ag + self.model.dot_wg 

        if len(self.model.order_4_gens) != 4:
            self.gen_f += self.model.dot_ed1 + self.model.dot_eq1 + self.model.dot_vf

        if self.model.add_governor_dynamics is True:
            self.gen_f += self.model.dot_pm

        if self.model.no_gfl_flag is False:
            self.conv_f += self.model.dot_ap_gfl + self.model.dot_wp_gfl 

        if self.model.no_gfm_flag is False:
            self.conv_f += self.model.dot_ap_gfm + self.model.dot_wp_gfm + self.model.dot_Ve + self.model.dot_E

        self.f = Matrix(self.gen_f + self.conv_f)
    
    def g_func(self):
        self.g1_eq = [self.model.vq[k] + self.model.Ra[k]*self.model.iq[k] \
                      + self.model.xd1[k]*self.model.id[k] - self.model.eq1[k]
                            for k in range(len(self.model.gen_idx)) if self.model.gen_order[k] == 4]
                                                                   
        self.g2_ed = [self.model.vd[k] + self.model.Ra[k]*self.model.id[k]\
                       - self.model.xq1[k]*self.model.iq[k] - self.model.ed1[k]
                                for k in range(len(self.model.gen_idx)) if self.model.gen_order[k] == 4]  

        # assuming Pg and Pl are all numbers g(.) = Pinj + Pl - Pg = 0 
        self.g = self.network_equations.P_inj + self.network_equations.Q_inj 
 
        if len(self.model.order_4_gens) != 0:
            self.g = self.g + self.g1_eq + self.g2_ed 
        

        self.g = Matrix(self.g) 

    def update_y_var(self):
        self.y = Matrix([i for i in self.y if i not in self.model.a_gen])
        if self.model.no_gfm_flag is False:
            self.y = Matrix([i for i in self.y if i not in self.model.a_gfm])
            self.y = Matrix([i for i in self.y if i not in self.model.v_conv])
