"""
FILE: Resistance_model.py
AUTHOR: R Ellingham
DATE: June 2021
PROGRAM DESC: Piecing together a few different components of the resistance-strain model

TODO:
1) Quasi static case (i.e. looking at steady state resistances in material and certain strains)
2)
"""

import numpy as np

# Quasi-static resistance model
def res_qstatic(m,x,c):
    return m*x + c

# Resistance relaxation model
def res_relax(t, b0, b1, tau_1, b2, tau_2, b3, tau_3):
    return b0 + b1 * np.exp(-t/tau_1) + b2 * np.exp(-t/tau_2) + b3 * np.exp(-t/tau_3)

# Resistance strain model without hysteresis
def res_strain(strain, a0, a1, a2):
    return a0*strain**2 + a1*strain + a2
    
## Resistance strain model with hysteresis
# Resistance with strain loading
def res_strain_load(strain, t, a0, a1, a2):
    return 

# Resistance without strain loading
def res_strain_unload(strain, t, a0, a1, a2):
    return
    
# When the material changes states the resistance of the material seems to dip
#   with a 2nd order polynomial relationship which is only valid for a certain 
#   time.
def res_strain_between_states(strain, t, a0, a1, a2):
    return
    
# Combining both loading and unloading models
def res_strain_nonlinear():
    # Determine where on the strain-resistance curve we are
    # Apply the relaxation history to material
    
