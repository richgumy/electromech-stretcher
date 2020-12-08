"""
FILE: clamp_calcs.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: A model for the spring loaded clamp mechanism
TODO: All
"""
import numpy as np
import matplotlib.pyplot as plt

# Spring models:
#   (1) Plastic end-loaded cantilever spring
num_springs = 4 # on each side of clamp
theta_pla = 0.608 #[rad] spring angle relative to compression (|theta/ <--F)
E_pla = 2.5e9 #[Pa] PLA elestic modulus generally 2.5 - 3.5 GPa (but unknown variation with 3d printed structures in different planes)
I_pla = (6e-3*1e-3*3)/12 #[m^4] second (cross-sectional) moment of area
L_pla = 7e-3 #[m] cantilever length
K_pla = (num_springs*3*E_pla*I_pla)/(L_pla*np.cos(theta_pla)) #[N/m]
x_pla_i = 4e-3
# dx_pla #[m] clamp spring compression

#   (2) Elastomer specimen spring
E_spec = 150e3 #[Pa] specimen elastic modulus
A_spec = 20e-3*20e-3 #[m^2] clamped area of specimen
L_spec = 3e-3 #[m]
K_spec = A_spec*E_spec/L_spec #[N/m]
x_spec_i = 4e-3
# dx_spec #[m] specimen spring compression

K_series = 1/(1/K_spec + 1/K_pla + 1/K_pla) #[N/m] series equivalent spring const

x_tot_i = x_pla_i*2 + x_spec_i
dx_tot = 4e-3
F_cl = K_series*dx_tot

dx_spec = F_cl/K_spec
dx_pla = F_cl/K_pla

print(F_cl,dx_spec,dx_pla, dx_spec+2*dx_pla)

# Tensile stress
F_t_max = 100 #[N]
A_spec_sec = 10e-3*4e-3 #[m^2]
possion_rat = 0.498

for F_t in range(1,F_t_max):
    stress_t = F_t/A_spec_sec
    strain_t = stress_t/E_spec
    strain_perp = possion_rat*strain_t
    x_spec_f = x_spec_i*(1+strain_perp)
    print(x_spec_f)
