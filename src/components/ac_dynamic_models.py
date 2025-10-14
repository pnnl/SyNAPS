from sympy import *
import pandas as pd
from sympy import init_printing
# Initialize pretty printing
init_printing(use_unicode=True)

pd.options.display.max_columns = None
pd.options.display.max_rows = None

class AC_Dynamic_models(object):
    '''
    The class models the AC dynamic components -- SG, GFL, and GFM.
    
    Inputs:
    -------
    - system: AC-network dataframes.
    - bus: bus model. 

    retruns differential equations for SG, GFL and GFMs.
    '''
    
    def __init__(self, system, bus):
        self.system = system
        self.bus = bus
        self.nodes = list(self.system.bus['idx'])    
        self.gen_name = list(self.system.gen['name'])
        self.gen_nodes = list(self.system.gen['bus'])
        self.gen_idx = list(self.system.gen['idx'])
        self.gen_uid = list(self.system.gen['uid'])
        self.gfm_idx = list(self.system.gfm['idx']) 
        self.gfl_idx = list(self.system.gfl['idx'])
        self.conv_idx = list(range(len(self.gfl_idx)+len(self.gfm_idx)))
        self.gfl_nodes = list(self.system.gfl['bus'])
        self.gfm_nodes = list(self.system.gfm['bus'])
        self.gen_pars = self.system.gen_pars
        self.gfl_pars = self.system.gfl_pars
        self.gfm_pars = self.system.gfm_pars

        '''
        The current SG implementation supports order-2 (GENCLS) and order-4 (GENROU) models.
        Exciter and governor dynamics are added for order-4 model.
        order-2 model includes only the swing equation. 
        '''
        self.order_4_gens = [self.gen_uid[k] for k in range(len(self.gen_idx)) \
                             if 'GENROU' in self.gen_name[k].split('_')[0]]
        self.order_2_gens = [self.gen_uid[k] for k in range(len(self.gen_idx)) \
                        if 'GENCLS' in self.gen_name[k].split('_')[0]]
        self.gen_order = [2 if x in self.order_2_gens else 4\
                           for x in self.order_2_gens + self.order_4_gens]
        
        self.add_governor_dynamics = True # only for 4th order model, not reqd for 2nd order, always True
        

        
        if len(self.gfm_idx) == 0:
            self.no_gfm_flag = True
        else:
            self.no_gfm_flag = False

        if len(self.gfl_idx) == 0:
            self.no_gfl_flag = True
        else:
            self.no_gfl_flag = False


    def SG_vars(self):
        '''
        Defining SG vars with symbolic variables.
        '''
        self.a_gen = [var(f"delta_g{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.w_gen = [var(f"omega_g{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.v_gen = [var(f"e_g{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.P_mech = [var(f"P_m{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.P_e = [var(f"P_e{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.P_g = [var(f"P_g{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.Q_g = [var(f"Q_g{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.eq1 = [var(f"e'_q{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.ed1 = [var(f"e'_d{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.iq = [var(f"i_q{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.id = [var(f"i_d{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.vf = [var(f"v_f{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.vq = [var(f"v_q{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
        self.vd = [var(f"v_d{self.gen_idx[k]}") for k in range(len(self.gen_idx))]
    
    def SG_params(self):
        '''
        Defining SG params with symbolic variables for 'sym' and 'user' or with numeric values
        for 'num' setting.
        '''
        self.M_gen = [
            var(f"M_g{self.gen_idx[k]}") if self.gen_pars['M'][0] in ['sym', 'user']
            else list(self.system.gen['M'])[k]
            for k in range(len(self.gen_idx))
        ]
        self.D_gen = [
            var(f"D_g{self.gen_idx[k]}") if self.gen_pars['D'][0] in ['sym', 'user']
            else list(self.system.gen['D'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.xd = [
            var(f"x_d{self.gen_idx[k]}") if self.gen_pars['xd'][0] in ['sym', 'user']
            else list(self.system.gen['xd'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.xq = [
            var(f"x_q{self.gen_idx[k]}") if self.gen_pars['xq'][0] in ['sym', 'user']
            else list(self.system.gen['xq'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.xd1 = [
            var(f"x'_d{self.gen_idx[k]}") if self.gen_pars['xd1'][0] in ['sym', 'user']
            else list(self.system.gen['xd1'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.xq1 = [
            var(f"x'_q{self.gen_idx[k]}") if self.gen_pars['xq1'][0] in ['sym', 'user']
            else list(self.system.gen['xq1'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.Td01 = [
            var(f"T'_d0{self.gen_idx[k]}") if self.gen_pars['Td10'][0] in ['sym', 'user']
            else list(self.system.gen['Td10'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.Tq01 = [
            var(f"T'_q0{self.gen_idx[k]}") if self.gen_pars['Tq10'][0] in ['sym', 'user']
            else list(self.system.gen['Tq10'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.Ra = [
            var(f"r_a{self.gen_idx[k]}") if self.gen_pars['ra'][0] in ['sym', 'user']
            else list(self.system.gen['ra'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.T_mech = [
            var(f"T_m{self.gen_idx[k]}") if self.gen_pars['Tm'][0] in ['sym', 'user']
            else list(self.system.tgov['Tm'])[k]
            for k in range(len(self.gen_idx))
        ]

        self.R_gov = [
            var(f"R_g{self.gen_idx[k]}") if self.gen_pars['Rgov'][0] in ['sym', 'user']
            else list(self.system.tgov['Rgov'])[k]
            for k in range(len(self.gen_idx))
        ]
        
    def SG_dynamics(self):
        '''
        This function creates the SG dynamic equations.
        '''
        self.dot_ag = [(self.w_gen[k] - 1)*self.system.wb for k in range(len(self.gen_idx))]

        self.dot_wg = [(1/self.M_gen[k])*(self.D_gen[k]*(1-self.w_gen[k])\
                                        + self.P_g[k] - self.P_e[k]) \
                                            for k in range(len(self.gen_idx))]  
                    
        self.dot_ed1 = [(1/self.Td01[k]) * (-self.eq1[k] - \
                                           (self.xd[k] - self.xd1[k])*self.id[k] + self.vf[k])\
                                                        for k in range(len(self.gen_idx))\
                                                              if self.gen_order[k] == 4]
        self.dot_eq1 = [(1/self.Tq01[k]) * (-self.ed1[k] - \
                                           (self.xq[k] - self.xq1[k])*self.iq[k])\
                                           for k in range(len(self.gen_idx))\
                                                              if self.gen_order[k] == 4]  
 

        for k in self.order_4_gens:
            self.dot_wg[k] = (1/self.M_gen[k])*(self.D_gen[k]*(1-self.w_gen[k])\
                                            + self.P_mech[k] \
                                                - (self.id[k]*(self.vd[k]+self.Ra[k]*self.id[k])\
                                                    + self.iq[k]*(self.vq[k]+self.Ra[k]*self.iq[k])))
        self.dot_pm = [(1/self.T_mech[k]) * ( self.P_g[k] \
                                    - self.P_mech[k] - (1/self.R_gov[k])\
                                            * (self.w_gen[k] - 1) )\
                                                for k in range(len(self.gen_idx))\
                                                        if self.gen_order[k] == 4]     
 
    def excitor_dynamics(self):
        '''
        This function creates the excitor dynamic equations.
        '''
        self.Te = [
        var(f"T_e{self.gen_idx[k]}") if self.gen_pars['Te'][0] in ['sym', 'user']
        else list(self.system.exc['Te'])[k]
        for k in range(len(self.gen_idx))
        ]
        self.Ke = [
        var(f"K_e{self.gen_idx[k]}") if self.gen_pars['Ke'][0] in ['sym', 'user']
        else list(self.system.exc['Ke'])[k]
        for k in range(len(self.gen_idx))
        ]
        self.vref = [
        var(f"V_ref{self.gen_idx[k]}") if self.gen_pars['Vref'][0] in ['sym', 'user']
        else list(self.system.exc['Vref'])[k]
        for k in range(len(self.gen_idx))
        ]            

        self.dot_vf = [(1/self.Te[k])*(self.Ke[k]*(self.vref[k] \
                                    - self.bus.v[self.gen_nodes[k]-1]) - self.vf[k])\
                                            for k in range(len(self.gen_idx))\
                                                        if self.gen_order[k] == 4] 

    def GFL_vars(self):
        '''
        Defining GFL vars with symbolic variables.
        '''
        self.a_gfl = [var(f"delta_pll{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.w_gfl = [var(f"omega_pll{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.a_gfl_E = [var(f"delta_c{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.e_gfl_E = [var(f"e_c{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.P_ref = [var(f"P_ref_{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.Q_ref = [var(f"Q_ref_{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.gfl_iq = [var(f"i_cq{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.gfl_id = [var(f"i_cd{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.zeta_gfl = [var(f"zeta{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]
        self.kappa_gfl = [var(f"kappa{self.gfl_idx[k]}") for k in range(len(self.gfl_idx))]

    def GFL_params(self):
        '''
        Defining GFL params with symbolic variables for 'sym' and 'user' or with numeric values
        for 'num' setting.
        '''
        self.k_droop = [
            var(f"k_{self.gfl_idx[k]}") if self.gfl_pars['kdroop'][0] in ['sym', 'user']
            else list(self.system.gfl['kdroop'])[k]
            for k in range(len(self.gfl_idx))
            ]
        self.kpPLL = [
            var(f"K_ppll{self.gfl_idx[k]}") if self.gfl_pars['KpPLL'][0] in ['sym', 'user']
            else list(self.system.gfl['KpPLL'])[k]
            for k in range(len(self.gfl_idx))
            ]
        self.kiPLL = [
            var(f"K_ipll{self.gfl_idx[k]}") if self.gfl_pars['KiPLL'][0] in ['sym', 'user']
            else list(self.system.gfl['KiPLL'])[k]
            for k in range(len(self.gfl_idx))
            ]
        self.kpc = [
            var(f"K_pc{self.gfl_idx[k]}") if self.gfl_pars['Kp'][0] in ['sym', 'user']
            else list(self.system.gfl['Kp'])[k]
            for k in range(len(self.gfl_idx))
            ]   
 
        self.kic = [
            var(f"K_ic{self.gfl_idx[k]}") if self.gfl_pars['Ki'][0] in ['sym', 'user']
            else list(self.system.gfl['Ki'])[k]
            for k in range(len(self.gfl_idx))
            ] 
        self.gfl_cap = [
            var(f"S_l{self.gfl_idx[k]}") if self.gfl_pars['Sl'][0] in ['sym', 'user']
            else list(self.system.gfl['Sl'])[k]
            for k in range(len(self.gfl_idx))
            ]          

        self.Xl = [
            var(f"X_l{self.gfl_idx[k]}") if self.gfl_pars['xq'][0] in ['sym', 'user']
            else list(self.system.gfl['xq'])[k]
            for k in range(len(self.gfl_idx))
            ]  
        
        self.K = [
            var(f"K_c{self.gfl_idx[k]}") if self.gfl_pars['kxl'][0] in ['sym', 'user']
            else list(self.system.gfl['kxl'])[k]
            for k in range(len(self.gfl_idx))
            ]  

        self.gfl_Ra = [
            var(f"r_a{self.gfl_idx[k]}") if self.gfl_pars['ra'][0] in ['sym', 'user']
            else list(self.system.gfl['ra'])[k]
            for k in range(len(self.gfl_idx))
            ]  

    def GFL_dynamics(self):
        '''
        This function creates the GFL dynamic equations.
        '''
        self.dot_ap_gfl = [(self.w_gfl[k] -1)*self.system.wb for k in range(len(self.gfl_idx))]
        self.dot_wp_gfl = [self.kpPLL[k]*self.bus.v[self.gfl_nodes[k]-1]\
                    *cos(self.bus.a[self.gfl_nodes[k]-1]-\
                        self.a_gfl[k])*(1-self.w_gfl[k])*self.system.wb \
                        + self.kiPLL[k]*self.bus.v[self.gfl_nodes[k]-1]*sin(self.bus.a[self.gfl_nodes[k]-1]-\
                            self.a_gfl[k]) for k in range(len(self.gfl_idx))]


    def GFM_vars(self):
        '''
        Defining GFM vars with symbolic variables.
        '''
        self.a_gfm = [var(f"delta_i{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.w_gfm = [var(f"omega_i{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.v_conv_err = [var(f"V_ei{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.v_conv = [var(f"v_i{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.Q_c = [var(f"Q_c{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.P_c = [var(f"P_c{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.P_set = [var(f"P_set{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.Q_set = [var(f"Q_set{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]
        self.V_set = [var(f"V_set{self.gfm_idx[k]}") for k in range(len(self.gfm_idx))]


    def GFM_params(self):
        '''
        Defining GFM params with symbolic variables for 'sym' and 'user' or with numeric values
        for 'num' setting.
        '''
        self.mp = [
            var(f"m_p{self.gfm_idx[k]}") if self.gfm_pars['mp'][0] in ['sym', 'user']
            else list(self.system.gfm['mp'])[k]
            for k in range(len(self.gfm_idx))
            ]
        self.mq = [
            var(f"m_q{self.gfm_idx[k]}") if self.gfm_pars['mq'][0] in ['sym', 'user']
            else list(self.system.gfm['mq'])[k]
            for k in range(len(self.gfm_idx))
            ]        
        self.cap = [
            var(f"S_c{self.gfm_idx[k]}") if self.gfm_pars['Sl'][0] in ['sym', 'user']
            else list(self.system.gfm['Sl'])[k]
            for k in range(len(self.gfm_idx))
            ]
        self.kpv = [
            var(f"K_pv{self.gfm_idx[k]}") if self.gfm_pars['kpv'][0] in ['sym', 'user']
            else list(self.system.gfm['kpv'])[k]
            for k in range(len(self.gfm_idx))
            ]    
        self.kiv = [
            var(f"K_iv{self.gfm_idx[k]}") if self.gfm_pars['kiv'][0] in ['sym', 'user']
            else list(self.system.gfm['kiv'])[k]
            for k in range(len(self.gfm_idx))
            ]  
        self.tau = [
            var(f"tau_{self.gfm_idx[k]}") if self.gfm_pars['tau'][0] in ['sym', 'user']
            else list(self.system.gfm['tau'])[k]
            for k in range(len(self.gfm_idx))
            ]   

    def GFM_dynamics(self):
        '''
        This function creates the REGFM-A1 dynamic equations.
        '''
        self.dot_ap_gfm = [(self.w_gfm[k] -1)*self.system.wb for k in range(len(self.gfm_idx))]
        self.dot_wp_gfm =  [(1/self.tau[k])*(1 - self.w_gfm[k] + \
                                            (self.mp[k]/self.cap[k])*(self.P_set[k] - self.P_c[k])) for k in range(len(self.gfm_idx))]
        self.dot_Ve = [(1/self.tau[k]) * (self.V_set[k] - self.v_conv_err[k] - self.bus.v[self.gfm_nodes[k]-1] + (self.mq[k]/self.cap[k]) * \
                                        (self.Q_set[k] - self.Q_c[k]) ) for k in range(len(self.gfm_idx))]
        
        self.dot_E  = [self.kpv[k] * (1/self.tau[k]) * (self.V_set[k] - self.v_conv_err[k] - self.bus.v[self.gfm_nodes[k]-1] + (self.mq[k]/self.cap[k]) * \
                                        (self.Q_set[k] - self.Q_c[k]) ) + self.kiv[k] * self.v_conv_err[k] for k in range(len(self.gfm_idx))]
                           
    def gen_var_classification(self):
        '''
        This function creates the SG's x vars for DAE formation
        '''
        self.gen_x = self.a_gen + self.w_gen \
                + [self.ed1[i] for i in self.order_4_gens] \
                        +[self.eq1[i] for i in self.order_4_gens] + [self.vf[i] for i in self.order_4_gens] \
                            + [self.P_mech[i] for i in self.order_4_gens if self.add_governor_dynamics is True]  

    def conv_var_classification(self):
        '''
        This function creates the Inverter's x vars for DAE formation
        '''
        self.conv_x = []
        self.conv_y = []
        if self.no_gfl_flag is False:
            self.conv_x = self.a_gfl + self.w_gfl 
            
        if self.no_gfm_flag is False:
            self.conv_x = self.conv_x + self.a_gfm + self.w_gfm + self.v_conv_err + self.v_conv