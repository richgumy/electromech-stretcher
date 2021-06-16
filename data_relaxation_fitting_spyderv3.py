"""
FILE: data_relaxation_fitting_spyderv2.py
AUTHOR: R Ellingham
DATE MODIFIED: Jun 2021
DATE CREATED: Dec 2020
PROGRAM DESC: Fit a relaxation model to stress relaxation data.

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

def extract_csv(input_filename)
    if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input

    with open(input_filename) as f:
        length_csv = sum(1 for line in f)
    
    Ro = np.zeros(length_csv)
    Ri = np.zeros(length_csv)
    tR = np.zeros(length_csv)
    P = np.zeros(length_csv)
    tP = np.zeros(length_csv)
    F = np.zeros(length_csv)
    tF = np.zeros(length_csv)
        
    with open(input_filename, 'r', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        data = list(data)
        i = 0
        for row in data:
                Ro[i] = float(row[0])
                Ri[i] = float(row[1])
                tR[i] = float(row[2])
                try:
                    P[i] = float(row[3])
                except ValueError:
                    print("Empty string found and approximated!")
                    P[i] = P[i-1]
                tP[i] = float(row[4])
                F[i] = float(row[5])
                tF[i] = float(row[6])
                print("%.2f" % (i/len(data)))
                i = i + 1
                # print(row[0],row[1],row[2],row[3],row[4],row[5])

# Test specimen dimensions (in m):
spec_length = 40e-3
spec_length_inner_electrodes = 20e-3
spec_width = 10e-3
spec_thickness = 4e-3

input_filename = "2_7-5_4Epin_20mm_quasistatic_v2.csv"

# Extract from csv
if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input

with open(input_filename) as f:
    length_csv = sum(1 for line in f)

Ro = np.zeros(length_csv)
Ri = np.zeros(length_csv)
tR = np.zeros(length_csv)
P = np.zeros(length_csv)
tP = np.zeros(length_csv)
F = np.zeros(length_csv)
tF = np.zeros(length_csv)


with open(input_filename, 'r', newline='') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    data = list(data)
    i = 0
    for row in data:
            Ro[i] = float(row[0])
            Ri[i] = float(row[1])
            tR[i] = float(row[2])
            try:
                P[i] = float(row[3])
            except ValueError:
                print("Empty string found and approximated!")
                P[i] = P[i-1]
            tP[i] = float(row[4])
            F[i] = float(row[5])
            tF[i] = float(row[6])
            print("%.2f" % (i/len(data)))
            i = i + 1
            # print(row[0],row[1],row[2],row[3],row[4],row[5])
            
t_lin, F_t, P_t, Ro_t, Ri_t = data_ems.interpolate_RFS_data(Ro,Ri,tR,P,tP,F,tF)

# Correct any erroneously gathered position values and scale to meters
P = data_ems.pos_outlier_corrector(P) * 1e-3 

# Use a MAF on the resistance data if it is an AC measurement
Ri_t = data_ems.MAF(Ri_t,12)

## Strain
# Calc engineering strain from displacement
Strain_t = P_t/(spec_length) #dx/x
# Calc true strain
Strain_true_t = np.log((P_t+spec_length)/(spec_length))

## Stress
# Calc engineering stress
A_0 = spec_width * spec_thickness
Stress_eng_t = F_t/A_0
# Calc stress from force and changing strain
poisson_ratio = 0.29 # Poisson's ratio. Found experimentally using #2_7.5%dragonskin10NV specimen
A_t = ((spec_width*spec_thickness)*(-Strain_t*poisson_ratio+1)*(-Strain_t*poisson_ratio+1))
Stress_pois_t = F_t/A_t # My approximation of true stress, with Poisson's ratio
# Common calc for the approximation of true strain
Stress_true_t = (F_t/A_0) * (1+Strain_t)

## Resistivity calc
resistivity_t = (Ri_t*A_t)/((1+Strain_true_t)*spec_length_inner_electrodes)

# Write interpolated and partially processed data to csv ##
# COMMENT OUT WHEN NOT REQUIRED ##
data_ems.write_processed_data(input_filename, Ro_t, Ri_t, resistivity_t, Strain_t, Strain_true_t, Stress_pois_t, Stress_eng_t, Stress_true_t, t_lin)

# Plot measurements over time
fig1, axs1 = plt.subplots(3, 1, constrained_layout=True,figsize = (10, 10))

ax = axs1[0]
ax.plot(t_lin, Ri_t,'r-')#),t_lin, Ri_t_fil,'b-')
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
ax.plot(t_lin, Strain_t,'r-',t_lin, Strain_true_t,'b-')
ax.set_title('')
ax.legend(("Eng","True"), loc='upper right')
ax.set_ylabel('Strain')
ax.set_xlabel('Time [s]')
ax.grid(True)

# ax = axs1[3]
# ax.plot(t_lin, resistivity_t,'r-')
# ax.set_title('')
# ax.set_ylabel('Resistivity [Ohm.m]') # Res-strain time delay may affect data significantly?
# ax.grid(True)

# Plot relaxation and fit curve
# Split data into piece-wise data
strain_splits = data_ems.split_ramp_data(Strain_t)
for i in range(len(strain_splits)):
    print(t_lin[strain_splits[i]])

# take chunks of stress and resistance relaxing values from index i1 to i2
i1 = []
i2 = []
for i in range(int(len(strain_splits)/4)):
    i1.append(int(4*i + 2))
    i2.append(int(4*i + 3))


# Curve fitting code (curve_fit func using non-lin lstsqr)
e0 = .10
# Viscoelastic models to fit data to 
def SLS_relax(t,E1,E2,C,mu):
    return (E1*E2)/(E1+E2) * e0 + C*np.exp(-((E1+E2)/mu)*t)

def SLS_relax_simple(t,E1,E2,C,mu):
    return E1 * e0 + C*np.exp(-(E2/mu)*t)

def generalisedi3(x, a0, a1, tau_1, a2, tau_2, a3, tau_3): # generalised Kelvin SLS relaxation model for n = 3
    return a0 + a1 * np.exp(-x/tau_1) + a2 * np.exp(-x/tau_2) + a3 * np.exp(-x/tau_3)

def generalisedi2(x, a0, a1, tau_1, a2, tau_2): # generalised Kelvin SLS relaxation model for n = 2
    return a0 + a1 * np.exp(-x/tau_1) + a2 * np.exp(-x/tau_2) 

# Stress modelling
consts_geni3 = []
pGen_init = [1,1,1,1,1]
# Resistance modelling
consts_geni3_R = []
pGen_init_R = [1,1,1,1,1]

start_offset = 6 # index offset to compensate for time lag between strain change and resistance/stress
end_offset = 1

quasi_stat_R = []
quasi_stat_S = []

try:
    # fig1, ax1 = plt.subplots()
    for i in range(1,int(len(strain_splits)/4)): 
        Res_load = Ri_t[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]
        Res_load_min = np.min(Res_load)
        Res_load = Res_load - Res_load_min

        Strain_load = Strain_t[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]

        Stress_load = Stress_pois_t[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]
        Stress_load_min = np.min(Stress_load)
        Stress_load = Stress_load - Stress_load_min

        t_load = t_lin[int(strain_splits[i1[i]])+start_offset:int(strain_splits[i2[i]])-end_offset]
        t_load = t_load - t_load[0]


        ## Levenberg–Marquardt algorithm for non-linear leastsq, fitting the stress data to generalised SLS relaxation models
        # Stress fitting
        poptS_gen, pcovS_gen = optimize.curve_fit(generalisedi2, t_load, Stress_load, p0=pGen_init, maxfev=50000)
        pGen_init = poptS_gen
        t_load_lin = np.linspace(min(t_load),max(t_load) , 100)
        Stress_load_lin_gen = generalisedi2(t_load_lin,*poptS_gen) + Stress_load_min
        poptS_gen[0] = poptS_gen[0] + Stress_load_min
        consts_geni3.append(poptS_gen) # Store fitted parameters of each relaxation
        
        # Resistance fitting
        poptS_gen_R, pcovS_gen_R = optimize.curve_fit(generalisedi2, t_load, Res_load, p0=pGen_init_R, maxfev=50000)
        pGen_init_R = poptS_gen_R
        Res_load_lin_gen = generalisedi2(t_load_lin,*poptS_gen_R) + Res_load_min
        poptS_gen_R[0] = poptS_gen_R[0] + Res_load_min
        consts_geni3_R.append(poptS_gen_R) # Store fitted parameters of each relaxation
                
        fig2, ax = plt.subplots(figsize = (10, 5))
        ax2 = ax.twinx()
        ax.plot(t_load,Stress_load + Stress_load_min,'rx',t_load_lin,Stress_load_lin_gen,'y-')
        ax2.plot(t_load,Res_load + Res_load_min,'bx',t_load_lin,Res_load_lin_gen,'y-')
        ax.set_xlabel('Time[s]')
        ax.set_ylabel('Stress [Pa]',color='r')
        ax2.set_ylabel('Resistance [Ohm]',color='b')
        plt.tight_layout()
        
        quasi_stat_R.append(Res_load[-1] + Res_load_min)
        quasi_stat_S.append(Strain_load[-1])


        
    
finally:
    consts_geni3 = np.array(consts_geni3)
    consts_geni3 = np.transpose(consts_geni3)
    print(consts_geni3)
    
    # a_indices = [0,1,3,5] # a_x indices
    tau_indices = [2,4] # tau_x indices
    fig1, ax1 = plt.subplots()
    for i in (range(len(consts_geni3))):
    # for i in range(len(tau_indices)):
        const = consts_geni3[i]
        # const = const[0:17]
        # fig1, ax1 = plt.subplots()
        ax1.plot(const,'rx')
        ax1.grid(True)
        const_mean = sum(const)/len(const)
        print("mean",const_mean)
        const_std = stats.tstd(const)
        print("std",const_std)
        const_var = stats.tvar(const)
        print("std %",100*const_std/const_mean)
        
        # x_axis = np.arange(const_mean - const_mean, const_mean + const_mean, 0.01)
        # fig, ax = plt.subplots(1, 1, constrained_layout=True)
        # ax.plot(x_axis, stats.norm.pdf(x_axis,const_mean,const_std))
        # ax.grid(True)
        
    plt.show()
