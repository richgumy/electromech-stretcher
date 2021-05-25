# -*- coding: utf-8 -*-
"""
Created on Mon May 24 17:19:22 2021

@author: rel80
"""

import numpy as np


# Resistance relaxation model
def res_relax(t, b0, b1, tau_1, b2, tau_2, b3, tau_3):
    return b0 + b1 * np.exp(-t/tau_1) + b2 * np.exp(-t/tau_2) + b3 * np.exp(-t/tau_3)

# Resistance strain model without hysteresis
def res_strain(strain, a0, a1, a2):
    return 
    

def res_strain_nonlinear()