from sympy import *
import numpy as np
from src.utils.custom_YBus import *

class network_equations(object):
    '''
    The class models the AC network equations for power-flow studies.
    
    Inputs:
    -------
    - system: AC-network dataframes.
    - bus: bus model. 
    - model: SG, GFM, GFL models.

    retruns algebraic equations for the AC network.
    '''
    def __init__(self, system, bus, model):
        self.system = system
        self.bus = bus
        self.model = model
        self.baseMVA = system.baseMVA
        self.branch_info = np.array(self.system.line[['bus1','bus2','r','x','b']])
        self.gen_idx = self.model.gen_idx
        self.gfm_idx = self.model.gfm_idx
        self.gfl_idx = self.model.gfl_idx
        self.conv_idx = self.model.conv_idx
        self.branch_base_list = self.branch_info.copy() 
        self.conv_old_bus_no = []

    
    def Y_bus_base_network(self):
        '''
        Compute the Y_bus matrix for the base network
        '''
        bus_base_network = np.array(self.system.bus['idx']).reshape(len(self.system.bus['idx']),1)
        branch_base_list = self.branch_base_list.copy()
        branch_base_list[:,0] = branch_base_list[:,0] - 1
        branch_base_list[:,1] = branch_base_list[:,1] - 1
        Ybus_base,_,_ = makeYbus_custom(self.baseMVA, bus_base_network-1, branch_base_list)
        self.G_base = np.real(Ybus_base.toarray())
        self.B_base = np.imag(Ybus_base.toarray())
        return self.G_base,self.B_base

    def bus_extension(self):
        '''
        For buses with gen resources (SG/GFM/GFL), the internal buses are added to the 
        given network.
        '''
        self.bus_extended = np.concatenate((np.array(self.system.bus['idx']),\
                                            self.gen_new_bus_no,self.conv_new_bus_no),axis=0)
        self.bus_extended_oldno = np.concatenate((np.array(self.system.bus['idx']),\
                                            self.gen_old_bus_no,self.conv_old_bus_no),axis=0)       
        self.bus_mapping = np.zeros([len(self.bus_extended),2])
        self.bus_mapping[:,0] = self.bus_extended
        self.bus_mapping[:,1] = self.bus_extended_oldno

    def branch_extension(self):
        '''
        Adding branches for the newly added bus and old bus using SG/GFM/GFL r and x.
        '''
        self.gen_new_bus_no = [len(self.system.bus)+k+1 for k in range(len(self.system.gen['idx']))]
        self.conv_new_bus_no = [len(self.system.bus)+len(self.system.gen['idx'])+k+1 \
                           for k in range(len(self.model.conv_idx))]  
        self.gen_old_bus_no = list(self.system.gen['bus'])

        if self.model.no_gfl_flag is False:
            self.conv_old_bus_no = self.conv_old_bus_no + list(self.system.gfl['bus']) 

        if self.model.no_gfm_flag is False:
            self.conv_old_bus_no = self.conv_old_bus_no + list(self.system.gfm['bus']) 

        extended_branch = np.zeros([len(self.system.gen['idx']) + len(self.system.gfl['idx']) + len(self.system.gfm['idx']),5])

        self.extended_resource_name = np.empty([len(self.system.gen['idx']) + len(self.system.gfl['idx']) + len(self.system.gfm['idx']),4], dtype=object)

        for k in range(len(self.system.gen['idx'])):
            extended_branch[k,0] = self.gen_old_bus_no[k]
            extended_branch[k,1] = self.gen_new_bus_no[k]
            extended_branch[k,2] = self.system.gen['ra'][k]
            extended_branch[k,3] = self.system.gen['xd1'][k]
            self.extended_resource_name[k,0] = self.system.gen['name'][k]
            self.extended_resource_name[k,1] = self.system.gen['gamma'][k]
            self.extended_resource_name[k,2] = self.gen_old_bus_no[k]
            self.extended_resource_name[k,3] = self.gen_new_bus_no[k]

        if self.model.no_gfl_flag is False:
            for k in range(len(self.system.gen['idx']),len(self.system.gen['idx'])\
                        +len(self.system.gfl['idx'])):
                extended_branch[k,0] = self.conv_old_bus_no[k-len(self.system.gen['idx'])]
                extended_branch[k,1] = self.conv_new_bus_no[k-len(self.system.gen['idx'])]
                extended_branch[k,2] = self.system.gfl['ra'][k-len(self.system.gen['idx'])]
                extended_branch[k,3] =  self.system.gfl['xq'][k-len(self.system.gen['idx'])]*\
                    self.system.gfl['kxl'][k-len(self.system.gen['idx'])]
                self.extended_resource_name[k,0] = self.system.gfl['name'][k-len(self.system.gen['idx'])]
                self.extended_resource_name[k,1] = self.system.gfl['gamma'][k-len(self.system.gen['idx'])]
                self.extended_resource_name[k,2] = self.conv_old_bus_no[k-len(self.system.gen['idx'])]
                self.extended_resource_name[k,3] = self.conv_new_bus_no[k-len(self.system.gen['idx'])]


        if self.model.no_gfm_flag is False:
            for k in range(len(self.system.gen['idx'])+len(self.system.gfl['idx']),len(self.system.gen['idx'])\
                        + len(self.system.gfl['idx']) +len(self.system.gfm['idx'])):
                extended_branch[k,0] = self.conv_old_bus_no[k-len(self.system.gen['idx'])]
                extended_branch[k,1] = self.conv_new_bus_no[k-len(self.system.gen['idx'])]
                extended_branch[k,2] = self.system.gfm['ra'][k-len(self.system.gen['idx'])-len(self.system.gfl['idx'])]
                extended_branch[k,3] = self.system.gfm['xq'][k-len(self.system.gen['idx'])-len(self.system.gfl['idx'])]
                self.extended_resource_name[k,0] = self.system.gfm['name'][k-len(self.system.gen['idx'])-len(self.system.gfl['idx'])]
                self.extended_resource_name[k,1] = self.system.gfm['gamma'][k-len(self.system.gen['idx'])-len(self.system.gfl['idx'])]
                self.extended_resource_name[k,2] = self.conv_old_bus_no[k-len(self.system.gen['idx'])]
                self.extended_resource_name[k,3] = self.conv_new_bus_no[k-len(self.system.gen['idx'])]

        self.branch_info_extended = np.concatenate((self.branch_info,extended_branch),axis=0)


    def alg_var_extension(self):
        '''
        Adding additioanl vars for new buses.
        '''
        self.bus.v = self.bus.v + [self.model.v_gen[k] for k in range(len(self.gen_idx))]
        self.bus.a = self.bus.a + [self.model.a_gen[k] for k in range(len(self.gen_idx))]

        if self.model.no_gfl_flag is False:
            self.bus.v = self.bus.v + [self.model.e_gfl_E[k] for k in range(len(self.gfl_idx))]
            self.bus.a = self.bus.a + [self.model.a_gfl_E[k] for k in range(len(self.gfl_idx))]
        
        if self.model.no_gfm_flag is False:
            self.bus.v = self.bus.v + [self.model.v_conv[k] for k in range(len(self.gfm_idx))]
            self.bus.a = self.bus.a + [self.model.a_gfm[k] for k in range(len(self.gfm_idx))]

    def Y_bus(self):
        '''
        Creating Ybus for the extended network.
        '''        
        bus_reshaped = self.bus_extended.reshape(len(self.bus_extended),1)
        branch_info_extended = self.branch_info_extended.copy()
        branch_info_extended[:,0] = branch_info_extended[:,0] - 1
        branch_info_extended[:,1] = branch_info_extended[:,1] - 1
        Ybus_SG,_,_ = makeYbus_custom(self.baseMVA, bus_reshaped-1, branch_info_extended)
        self.G = np.real(Ybus_SG.toarray())
        self.B = np.imag(Ybus_SG.toarray())
        return self.G,self.B   

    def injections(self):
        '''
        Creating active and reactive power injection for forming algebraic equations.
        ''' 
        self.P_inj = []
        self.Q_inj = []
        for k in range(len(self.bus_extended)):
            p = 0
            q = 0
            for j in range(len(self.bus_extended)):
                p += self.bus.v[k]*self.bus.v[j]*(self.G[k,j]*cos(self.bus.a[k]-self.bus.a[j])\
                                                + self.B[k,j]*sin(self.bus.a[k]-self.bus.a[j]))
                q += self.bus.v[k]*self.bus.v[j]*(self.G[k,j]*sin(self.bus.a[k]-self.bus.a[j])\
                                                - self.B[k,j]*cos(self.bus.a[k]-self.bus.a[j]))            
            self.P_inj.append(p)
            self.Q_inj.append(q)
 
    def include_load(self):
        '''
        Including load as variable.
        ''' 
        for k in range(len(self.bus.nodes)):
            self.P_inj[k] = self.P_inj[k] + self.bus.pl[k] 
            self.Q_inj[k] = self.Q_inj[k] + self.bus.ql[k] 


    def include_gen(self):
        '''
        Including generation as variable.
        ''' 
        for k in range(len(self.system.gen)):
            r = np.where(self.extended_resource_name.astype(str)==list(self.system.gen['name'])[k])[0]
            new_bus_idx = self.extended_resource_name[r,3][0] - 1
            self.Q_inj[new_bus_idx] = self.Q_inj[new_bus_idx] - self.model.Q_g[k]


    def update_gen_conv_dynamics(self):
        '''
        Updating SG and GFM dynamics for including P,Q and remove respective qts. from algebraic eqns.
        ''' 
        self.model.vq = [self.bus.v[self.model.gen_nodes[k]-1]*cos(self.model.a_gen[k]\
                            -self.bus.a[self.model.gen_nodes[k]-1]) for k in range(len(self.gen_idx))]
        
        self.model.vd = [self.bus.v[self.model.gen_nodes[k]-1]*sin(self.model.a_gen[k]\
                            -self.bus.a[self.model.gen_nodes[k]-1]) for k in range(len(self.gen_idx))]
                                                                        
            
        self.model.P_e = [self.P_inj[k+len(list(self.system.bus['idx']))] for k in range(len(self.gen_idx))] #this P_inj is only generation with no load
    
        # eliminate gen Pinjs
        self.P_inj = [i for i in self.P_inj if i not in self.model.P_e]  

        # update gen P_e for 4th order
        for k in range(len(self.gen_idx)):
            if self.model.gen_order[k]==4:
                self.model.P_e[k] = self.model.id[k]*(self.model.vd[k] + self.model.Ra[k]*self.model.id[k])\
                + self.model.iq[k]*(self.model.vq[k] + self.model.Ra[k]*self.model.iq[k])

 
        self.model.dot_wg = [(1/self.model.M_gen[k])*(self.model.D_gen[k]*(1-self.model.w_gen[k])\
                                        + self.model.P_g[k] - self.model.P_e[k]) \
                                            for k in range(len(self.gen_idx))]       
        for k in self.model.order_4_gens:
            self.model.dot_wg[k] = (1/self.model.M_gen[k])*(self.model.D_gen[k]*(1-self.model.w_gen[k])\
                            + self.model.P_mech[k] - self.model.P_e[k]) 

        if self.model.no_gfm_flag is False:
            self.model.P_c = [self.P_inj[k+len(list(self.system.bus['idx']))+len(list(self.gfl_idx))]\
                            for k in range(len(self.gfm_idx))] #this P_inj is only generation with no load
                        
            self.model.Q_c = [self.Q_inj[k+len(list(self.system.bus['idx']))+len(list(self.gen_idx))+len(list(self.gfl_idx))]\
                            for k in range(len(self.gfm_idx))]
            
            # eliminate gfm Pinj and Qinj
            self.P_inj = [i for i in self.P_inj if i not in self.model.P_c]
            self.Q_inj = [i for i in self.Q_inj if i not in self.model.Q_c]
            
            self.model.dot_wp_gfm =  [(1/self.model.tau[k])*(1 - self.model.w_gfm[k] + \
                                    (self.model.mp[k]/self.model.cap[k])*(self.model.P_set[k] - self.model.P_c[k])) for k in range(len(self.gfm_idx))]

            self.model.dot_Ve = [(1/self.model.tau[k]) * (self.model.V_set[k] - self.model.v_conv_err[k] - self.bus.v[self.model.gfm_nodes[k]-1] + (self.model.mq[k]/self.model.cap[k])* \
                                         (self.model.Q_set[k] - self.model.Q_c[k]) ) for k in range(len(self.gfm_idx))]
            
            self.model.dot_E  = [self.model.kpv[k] * (1/self.model.tau[k]) * (self.model.V_set[k] - self.model.v_conv_err[k] - self.bus.v[self.model.gfm_nodes[k]-1] + (self.model.mq[k]/self.model.cap[k])* \
                                         (self.model.Q_set[k] - self.model.Q_c[k]) ) + self.model.kiv[k] * self.model.v_conv_err[k] for k in range(len(self.gfm_idx))]


