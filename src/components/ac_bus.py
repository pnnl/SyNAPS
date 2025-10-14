from sympy import *
import numpy as np
import pandas as pd

class Bus(object):
    '''
    The class creates the symbolic variables for all the buses.
    
    Inputs:
    -------
    - system: AC-network dataframes. 
    '''
    def __init__(self, system):
        self.system = system
        self.nodes = list(self.system.bus['idx'])    
        self.a = [var(f"delta{self.nodes[k]}") for k in range(len(self.system.bus))]
        self.v = [var(f"v{self.nodes[k]}") for k in range(len(self.system.bus))]
        self.pl = [var(f"P_l{self.nodes[k]}") for k in range(len(self.system.bus))]
        self.ql = [var(f"Q_l{self.nodes[k]}") for k in range(len(self.system.bus))]
        self.pVG = [var(f"P_VG{self.nodes[k]}") for k in range(len(self.system.bus))]
        self.qVG = [var(f"Q_VG{self.nodes[k]}") for k in range(len(self.system.bus))]