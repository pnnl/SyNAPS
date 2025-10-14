from sympy import *
import numpy as np
import cmath

def calculate_sending_end_voltage(V, Pg, Qg, ra, xq):
    '''
    This function computes the internal node V,theta,P,Q values for generating resources based 
    on the external bus values.
    '''
    # Step 1: Calculate the current at the receiving end (I_r)
    I_conj = (Pg + 1j * Qg) / V

    I = np.conjugate(I_conj)

    # Step 2: Calculate the voltage drop across the line (V_drop)
    V_drop = I * (ra + 1j*xq)

    # Step 3: Calculate the sending end voltage (V_s)
    Eq = V + V_drop

    Eq_mag = abs(Eq)
    Eq_angle = cmath.phase(Eq)

    V_mag = abs(V)
    V_angle = cmath.phase(V)


    Id0 = abs(I)*cos(V_angle - cmath.phase(I))
    Iq0 = - abs(I)*sin(V_angle - cmath.phase(I))

    Vd0 = Eq_mag*cos(Eq_angle - V_angle)
    Vq0 = Eq_mag*sin(Eq_angle - V_angle)

    Qg = Vq0*Id0 - Vd0*Iq0
    Pg = Vq0*Iq0 + Vd0*Id0

    return Eq_mag,Eq_angle,Id0,Iq0,Pg,Qg







