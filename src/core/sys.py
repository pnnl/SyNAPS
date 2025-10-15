import numpy as np
from src.IO.xcel import *
from src.utils.config_upd import *

class get_ac_model(object):
    '''
    The class reads the data file representing the AC network. The data file contains the different
    data sheets. 
    
    If user wants to change the default config setting, this class updates the gen/gfl/gfm parameters 
    default configs -- "sym", "num" and "user", as per user's requirement.
    
    baseMVA is set 100 as default.

    Inputs:
    -------
    - case_name: 
        excel file with AC network details.  
    - updates_gen: dict  
        user provided configs for gen parameters. 
    - updates_gfl: dict  
        user provided configs for gfl parameters. 
    - updates_gfm: dict  
        user provided configs for gfm parameters. 
    - baseMVA: 100 (default)

    returns AC network data sheets as pd (pandas) dataframes.
    '''
    
    def __init__(self,case_name,updates_gen={},updates_gfl={},updates_gfm={},baseMVA=100):
        self.bus = read_excel(case_name, sheet_name='Bus')
        self.pq = read_excel(case_name, sheet_name='PQ')
        self.pv = read_excel(case_name, sheet_name='PV')
        self.slack = read_excel(case_name, sheet_name='Slack')
        self.line = read_excel(case_name, sheet_name='Line')
        self.gen = read_excel(case_name, sheet_name='GEN')
        self.exc = read_excel(case_name, sheet_name='EXST1')
        self.tgov = read_excel(case_name, sheet_name='TGOV1')
        self.gen_pars = read_excel(case_name, sheet_name='GEN_params_status')
        self.gfl_pars = read_excel(case_name, sheet_name='Conv_GFL_params_status')
        self.gfm_pars = read_excel(case_name, sheet_name='Conv_GFM_params_status')
        self.gfl = read_excel(case_name, sheet_name='Conv_GFL')
        self.gfm = read_excel(case_name, sheet_name='Conv_GFM')
        self.wb = 2*np.pi*60
        self.baseMVA = baseMVA 

        print('Default status of SG params:', self.gen_pars)
        print('Default status of GFL params:', self.gfl_pars)
        print('Default status of GFM params:', self.gfm_pars)

        self.gen_pars = update_config(self.gen_pars,updates_gen)
        self.gfl_pars = update_config(self.gfl_pars,updates_gfl)
        self.gfm_pars = update_config(self.gfm_pars,updates_gfm)

        print('Current status of SG params:', self.gen_pars)
        print('Current status of GFL params:', self.gfl_pars)
        print('Current status of GFM params:', self.gfm_pars)

        self.gen_base_conv_factor = [list(self.gen['Sn'])[k]/self.baseMVA for k in range(len(self.gen))]
        print('SG correction factor',self.gen_base_conv_factor)
        self.gen['M'] = [list(self.gen['M'])[k]*self.gen_base_conv_factor[k] for k in range(len(self.gen))]
        self.gen['D'] = [list(self.gen['D'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]
        self.gen['xd'] = [list(self.gen['xd'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))] 
        self.gen['xq'] = [list(self.gen['xq'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]
        self.gen['xd1'] = [list(self.gen['xd1'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]   
        self.gen['xq1'] = [list(self.gen['xq1'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]  
        self.gen['xd2'] = [list(self.gen['xd2'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]   
        self.gen['xq2'] = [list(self.gen['xq2'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]  
        self.gen['ra'] = [list(self.gen['ra'])[k]/self.gen_base_conv_factor[k] for k in range(len(self.gen))]  

        # self.gfl_base_conv_factor = [list(self.gfl['Sn'])[k]/self.baseMVA for k in range(len(self.gfl))]
        # self.gfl['xq'] = [list(self.gfl['xq'])[k]/self.gfl_base_conv_factor[k] for k in range(len(self.gfl))]
        # self.gfl['ra'] = [list(self.gfl['ra'])[k]//self.gfl_base_conv_factor[k] for k in range(len(self.gfl))]
