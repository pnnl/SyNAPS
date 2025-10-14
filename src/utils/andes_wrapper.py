import andes
from andes.utils.paths import get_case
andes.config_logger()

from src.IO.file_move import *


def PF_sol(base_path,pf_file_name):
    '''
    The function acts as a wrapper for Andes PF solver.
    
    Inputs:
    -------
    - base_path: path of the case directory.
    - pf_file_name: ac case file name to be solved.

    retruns dict with pf solutions.
    '''
    ac_case_name_pf = os.path.join(base_path, pf_file_name)
    ss = andes.run(get_case(ac_case_name_pf),
               default_config=True)  
    andes_report_transfer(base_path)

    return ss

def TD_run(base_path,ss,tend):
    '''
    The function acts as a wrapper for Andes TD simulation.
    
    Inputs:
    -------
    - base_path: path of the case directory.
    - ss: output dict of PF solution.

    retruns delta,omega and time data as per our specific verification need. (this can be modified)
    '''

    ''' Constant P-Q load '''
    ss.PQ.config.p2p = 1.0
    ss.PQ.config.p2i = 0
    ss.PQ.config.p2z = 0

    ss.PQ.config.q2q = 1.0
    ss.PQ.config.q2i = 0
    ss.PQ.config.q2z = 0

    ''' Simulation time'''
    ss.TDS.config.tf = tend # simulate for 5 seconds
    # ss.TDS.config.method= "backeuler"

    ss.TDS.run()

    data_v = ss.TDS.plt.get_values((1,2,3,4,5,6,7,8,9,10))
    data_t = ss.TDS.plt.t.reshape(len(ss.TDS.plt.t),1)

    andes_report_transfer(base_path)

    return data_v,data_t