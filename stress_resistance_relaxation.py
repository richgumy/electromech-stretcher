"""
FILE: stress_resistance_relaxation.py
AUTHOR: R Ellingham
DATE MODIFIED: Jun 2021
DATE CREATED: Jun 2021
PROGRAM DESC: Matlab style script for fitting a relaxation model to stress and 
apparent resistance relaxation data.
NOTES: Made for SPyder IDE

TODO:
1)
2)
"""

import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
import scipy.stats as stats
import matplotlib
import data_ems

#%% Setup parameters and determine all stress and strain values

# Test specimen dimensions (in m):
spec_length = 40e-3
spec_length_inner_electrodes = 20e-3
spec_width = 10e-3
spec_thickness = 4e-3

input_filename = "2_7-5_E4pin_20mm_v6.csv"

# Extract experimental csv data
Ro,Ri,tR,P,tP,F,tF = data_ems.extract_csv(input_filename)

# Correct any erroneously gathered position values and scale to meters
P = data_ems.pos_outlier_corrector(P) * 1e-3 

# Interpolate data to get a constant smaple frequency across all measurements
t_lin, F_t, P_t, Ro_t, Ri_t = data_ems.interpolate_RFS_data(Ro,Ri,tR,P,tP,F,tF)

# Use a MAF on the resistance data if it is an 'AC' measurement
Ri_t = data_ems.MAF(Ri_t,12)

## Strain ##
# Calc engineering strain from displacement
Strain_eng_t = P_t/(spec_length) #dx/x
# Calc true strain
Strain_log_t = np.log((P_t+spec_length)/(spec_length))

## Stress ##
# Calc engineering stress
A_0 = spec_width * spec_thickness # initial cross-sectional area
Stress_eng_t = F_t/A_0
# Calc stress from force and changing cross-sectional area
poisson_ratio = 0.29 # Poisson's ratio. Found experimentally using #2_7.5%dragonskin10NV specimen
A_t = ((spec_width*spec_thickness)*(-Strain_eng_t*poisson_ratio+1)*(-Strain_eng_t*poisson_ratio+1))
Stress_pois_t = F_t/A_t # My approximation of true stress, with Poisson's ratio
# Common calc for the approximation of true strain
Stress_true_t = (F_t/A_0) * (1+Strain_eng_t)

# Write interpolated and partially processed data to csv ##
# data_ems.write_processed_data(input_filename, Ro_t, Ri_t, resistivity_t, Strain_eng_t, Strain_log_t, Stress_pois_t, Stress_eng_t, Stress_true_t, t_lin)
## ^UN/COMMENT OUT WHEN NOT REQUIRED^ ##

# Plot measurements over time
fig1, axs1 = plt.subplots(3, 1, constrained_layout=True,figsize = (10, 12))

ax = axs1[0]
ax.plot(t_lin, Ri_t,'r-')
ax.set_title('')
ax.set_ylabel('Resistance [Ohm]')
ax.grid(True)

ax = axs1[1]
ax.plot(t_lin, Stress_eng_t,'r-',t_lin, Stress_pois_t,'b-',t_lin, Stress_true_t,'g-')
ax.set_title('')
ax.legend(("Eng","Poisson","True"), loc='upper right')
ax.set_ylabel('Stress [Pa]')
ax.grid(True)

ax = axs1[2]
ax.plot(t_lin, Strain_eng_t,'r-',t_lin, Strain_log_t,'b-')
ax.set_title('')
ax.legend(("Eng","True"), loc='upper right')
ax.set_ylabel('Strain')
ax.grid(True)

ax.set_xlabel('Time [s]')

#%% Resistivity investigation

## Resistivity calc ##
resistivity_t = (Ri_t*A_t)/((1+Strain_log_t)*spec_length_inner_electrodes)

# Plot measurements over time
fig2, axs2 = plt.subplots(2, 1, constrained_layout=True,figsize = (10, 6))

ax = axs1[0]
ax.plot(t_lin, Ri_t,'r-')
ax.set_title('')
ax.set_ylabel('Resistance [Ohm]')
ax.grid(True)

ax = axs1[1]
ax.plot(t_lin, resistivity_t,'r-')
ax.set_title('')
ax.set_ylabel('Resistivity [Ohm.m]') # Res-strain time delay may affect data significantly?
ax.grid(True)


#%% Split data into separate pulses so each relaxation/loading-cycle can be analysed separately

strain_splits = data_ems.split_ramp_data(Strain_log_t)
for i in range(len(strain_splits)):
    print(t_lin[strain_splits[i]])

# take chunks of stress and resistance relaxing values from index i1 to i2
i1 = []
i2 = []
for i in range(int(len(strain_splits)/4)):
    i1.append(int(4*i + 2))
    i2.append(int(4*i + 3))

#%% Fitting a line to the stress-strain loading and unloading data

# Plot relaxation and fit curve

# Simple linear function to fit the stress-strain relationship
def lin_func(x,m,c):
    return m*x + c

consts_lin_fit = []
p_init = [1,1]

constsu_lin_fit = []
pu_init = [1,1]

fig3, ax3 = plt.subplots(figsize=(12,6))
# for i in range(1,int(len(strain_splits))):
for i in range(5):
    Strain_load = Strain_log_t[int(strain_splits[4*i+1]):int(strain_splits[4*i+2])] # Using log strain as a better representation of the strain of the material
    Strain_unload = Strain_log_t[int(strain_splits[4*i+3]):int(strain_splits[4*i+4])]

    Stress_load = Stress_pois_t[int(strain_splits[4*i+1]):int(strain_splits[4*i+2])]
    # Stress_load_min = np.min(Stress_load)
    # Stress_load = Stress_load - Stress_load_min
    Stress_unload = Stress_pois_t[int(strain_splits[4*i+3]):int(strain_splits[4*i+4])]
    # Stress_unload_min = np.min(Stress_load)
    # Stress_unload = Stress_load - Stress_load_min
    
    # Stress-strain loading fitting
    poptS, pcovS = optimize.curve_fit(lin_func, Strain_load, Stress_load, p0=p_init, maxfev=50000)
    p_init = poptS
    Strain_load_lin = np.linspace(min(Strain_load),max(Strain_load) , 100)
    Stress_load_lin = lin_func(Strain_load_lin,*poptS)
    consts_lin_fit.append(poptS) # Store fitted parameters of each relaxation
    
    # Stress-strain unloading fitting
    poptSu, pcovSu = optimize.curve_fit(lin_func, Strain_unload, Stress_unload, p0=pu_init, maxfev=50000)
    pu_init = poptSu
    Strain_unload_lin = np.linspace(min(Strain_unload),max(Strain_unload) , 100)
    Stress_unload_lin = lin_func(Strain_load_lin,*poptSu)
    constsu_lin_fit.append(poptSu) # Store fitted parameters of each relaxation
    
    ## Plot stress-strain loading/unloading of specimen
    # fig3, ax3 = plt.subplots(figsize=(16,8))
    # Loading
    ax3.plot(100*Strain_load,Stress_load,color='r',marker='x',ls='')
    # ax3.plot(Strain_load_lin,Stress_load_lin,color='y',ls='-')
    # Unloading
    ax3.plot(100*Strain_unload,Stress_unload,color='b',marker='x',ls='')
    # ax3.plot(Strain_unload_lin,Stress_unload_lin,color='g',ls='-')
    ax3.legend(["Loading", "Unloading"])
    ax3.set_ylabel('Stress [Pa]')
    ax3.set_xlabel('Strain [%]')


#%% Fitting models to the stress relaxation data

# Curve fitting code (curve_fit func using non-lin lstsqr)
e0 = .10 # We should b
# e able to input the elastic modulus parameter found in the previous section

# Viscoelastic models to fit data to 
def SLS_relax(t,E1,E2,C,mu):
    return (E1*E2)/(E1+E2) * e0 + C*np.exp(-((E1+E2)/mu)*t)

def SLS_relax_simple(t,E1,E2,C,mu):
    return E1 * e0 + C*np.exp(-(E2/mu)*t)

def generalised_SLS_2e(x, a0, a1, tau_1, a2, tau_2): # generalised Kelvin SLS relaxation model for n = 2
    return a0 + a1 * np.exp(-x/tau_1) + a2 * np.exp(-x/tau_2) 

def generalised_SLS_3e(x, a0, a1, tau_1, a2, tau_2, a3, tau_3): # generalised Kelvin SLS relaxation model for n = 3
    return a0 + a1 * np.exp(-x/tau_1) + a2 * np.exp(-x/tau_2) + a3 * np.exp(-x/tau_3)

## Fit the stress relaxation single element SLS model

# Stress modelling
consts_SLS = []
pSLS_init = [1,1,1,1] # initial guess of parameters

start_offset = 6 # index offset to compensate for time lag between strain change and resistance/stress
end_offset = 1

for i in range(1,int(len(strain_splits)/4)): 
    
        Strain_load = Strain_log_t[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]

        Stress_load = Stress_pois_t[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]
        Stress_load_min = np.min(Stress_load)
        Stress_load = Stress_load - Stress_load_min

        t_load = t_lin[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]
        t_load = t_load - t_load[0]
        
        ## Levenberg–Marquardt algorithm for non-linear leastsq, fitting the stress data to generalised SLS relaxation models
        # Stress fitting
        poptS_SLS, pcovS_SLS = optimize.curve_fit(SLS_relax_simple, t_load, Stress_load, p0=pSLS_init, maxfev=50000)
        
        pSLS_init = poptS_SLS # Have the next initial guess of the next relaxation equal the previous relaxation's fitted parameters
        
        # Make curve for fitted model
        t_load_lin = np.linspace(min(t_load),max(t_load) , 100)
        Stress_load_lin_SLS = SLS_relax_simple(t_load_lin,*poptS_SLS) + Stress_load_min
        poptS_SLS[0] = poptS_SLS[0] + Stress_load_min
        consts_SLS.append(poptS_SLS) # Store fitted parameters of each relaxation
        
        # Determine the goodness of fit
        Error = Stress_load - SLS_relax_simple(t_load,*poptS_SLS) + Stress_load_min
        
        fig4, ax4 = plt.subplots(figsize = (10, 5))
        ax.plot(t_load,Error,'r-')
        ax.set_xlabel('Time[s]')
        ax.set_ylabel('Stress [Pa]',color='r')
        
        # Plot relaxation data against
        fig5, ax5 = plt.subplots(figsize = (10, 5))
        ax.plot(t_load,Stress_load + Stress_load_min,'rx',t_load_lin,Stress_load_lin_SLS,'y-')
        ax.set_xlabel('Time[s]')
        ax.set_ylabel('Stress [Pa]',color='r')
