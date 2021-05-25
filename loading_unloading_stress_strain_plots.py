"""
FILE: data_relaxation_fitting_spyderv2.py
AUTHOR: R Ellingham
DATE: Dec 2020
PROGRAM DESC: Fit a relaxation model to stress relaxation data.

TODO:
1)
2)
"""

import csv
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

def split_ramp_data(data):
    """
    DESCR: Returns an array of indices showing where all of strain rates reach zero
    EG for below ramp returns [3,5,9,11]
                      ...
                    /|
                ...  |
              /|   | |
    data->...  |   | |
             | |   | |
    indx->---3-5---9-11--
    IN_PARAMS: Ramp profile data
    NOTES:
    TODO:
    """
    index_splits = [0]
    n = len(data)
    for i in range(1,n-1):
        if (data[i-1] != data[i] and data[i] == data[i+1]) or (data[i-1] == data[i] and data[i] != data[i+1]):
            index_splits.append(i)
    index_splits.append(n-1)
    return index_splits

def MAF(x,N):
    """
    DESCR: Moving average filter for a set of data averaging +/- N/2 each point
    IN_PARAMS: x data, sample period to average over
    NOTES: Requires numpy library.
    TODO:
    """
    x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
    x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')

    return x_smooth

def pos_outlier_corrector(strain_arr):
    new_strain_arr = strain_arr
    for i in range(1,len(strain_arr)):
        if abs(strain_arr[i] - strain_arr[i-1]) > 3:
            new_strain_arr[i] = strain_arr[i-1]
    return new_strain_arr

### MAIN CODE STARTS HERE ###

# Test specimen dimensions:
spec_length = 40e-3
spec_width = 10e-3
spec_thickness = 4e-3

input_filename = "1_CB0_v1_0.3Strain.csv"
# input_filename = "2_7-5_E4pin_20mm_v10_0.3Strain.csv"
# input_filename = "1_10_E4pin_20mm_v9_0.3Strain.csv"

# Extract from csv
if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input

with open(input_filename) as f:
    length_csv = sum(1 for line in f)

R = np.zeros(length_csv)
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
            R[i] = float(row[0])
            Ri[i] = float(row[1])
            tR[i] = float(row[2])
            P[i] = float(row[3])
            tP[i] = float(row[4])
            F[i] = float(row[5])
            tF[i] = float(row[6])
            print("%.2f" % (i/len(data)))
            i = i + 1
            # print(row[0],row[1],row[2],row[3],row[4],row[5])
            
P = pos_outlier_corrector(P)

# Interpolate all data
Fintp_P = np.interp(tP,tF,F)
Fintp_R = np.interp(tR,tF,F)

Pintp_F = np.interp(tF,tP,P)
Pintp_R = np.interp(tR,tP,P)

# timing data for resistance for inital testing is averaged (i.e. R and Ri measurements assumed to be taken at same time)
Rintp_P = np.interp(tP,tR,R)
Rintp_F = np.interp(tF,tR,R)

Riintp_P = np.interp(tP,tR,Ri)
Riintp_F = np.interp(tF,tR,Ri)

# Generate new interpolated data for the '_tot' arrays
# interpolation order is [force,pos,res]
F_tot = np.append(F,[Fintp_P,Fintp_R])
tF_tot = np.append(tF,[tP,tR])

P_tot = np.append(Pintp_F,[P,Pintp_R])
tP_tot = np.append(tF,[tP,tR])

R_tot = np.append(Rintp_F,[Rintp_P,R])
tR_tot = np.append(tF,[tP,tR])

Ri_tot = np.append(Riintp_F,[Riintp_P,Ri])

## Sort all _tot data into time order
# Make new arrays for sorted data
F_ = np.zeros(len(tF_tot))
P_ = np.zeros(len(tF_tot))
R_ = np.zeros(len(tF_tot))
Ri_ = np.zeros(len(tF_tot))
t_ = np.zeros(len(tF_tot))

ind = np.argsort(tF_tot)
for i in range(len(tF_tot)):
    F_[i] = F_tot[ind[i]]
    P_[i] = P_tot[ind[i]]
    R_[i] = R_tot[ind[i]]
    Ri_[i] = Ri_tot[ind[i]]
    t_[i] = tF_tot[ind[i]]

# Interpolate again to get time in a linear form
t_lin = np.linspace(min(t_),max(t_),int(len(t_)/3)) # divide by 3 to get original size of data
F_t = np.interp(t_lin,t_,F_)
P_t = np.interp(t_lin,t_,P_)
R_t = np.interp(t_lin,t_,R_)
Ri_t = np.interp(t_lin,t_,Ri_)

# Use a MAF on the resistance data if it is an AC measurement
Ri_t = MAF(Ri_t,8)

# Calc strain from displacement
Strain_t = -P_t/(spec_length*1000) # dx/x

# Calc stress from force and changing strain
poisson_ratio = 0.29 # Poisson's ratio. Found experimentally using #2_7.5%dragonskin10NV specimen
Stress_t = F_t/((spec_width*spec_thickness)*(-Strain_t*poisson_ratio+1)*(-Strain_t*poisson_ratio+1)) # force/cross-sectional-area

# Determine strain gauge factor dRes/dStrain
Gauge_factor = ((Ri_t - min(Ri_t))/Ri_t)
Gauge_factor_t = np.zeros(len(Ri_t))
for i in range(len(Gauge_factor)):
    if Strain_t[i] == 0:
        Gauge_factor_t[i] = 0
    else:
        Gauge_factor_t[i]  = Gauge_factor[i] / Strain_t[i]

# Plot measurements over time
fig1, axs1 = plt.subplots(4, 1, constrained_layout=True)

ax = axs1[0]
ax.plot(t_lin, R_t,'r-')
ax.set_title('')
ax.set_ylabel('Resistance outer [Ohm]')
ax.grid(True)

ax = axs1[1]
ax.plot(t_lin, Ri_t,'b-')#),t_lin, Ri_t_fil,'b-')
ax.set_title('')
ax.set_ylabel('Resistance inner [Ohm]')
ax.set_xlabel('Time [s]')
ax.grid(True)

ax = axs1[2]
ax.plot(t_lin, Stress_t,'r-')
ax.set_title('')
ax.set_ylabel('Stress [Pa]')
ax.grid(True)

ax = axs1[3]
ax.plot(t_lin, Strain_t,'r-')
ax.set_title('')
ax.set_ylabel('Strain')
ax.grid(True)

## Plot relaxation and fit curve
# Split data into piece-wise data
strain_splits = split_ramp_data(Strain_t)
for i in range(len(strain_splits)):
    print(t_lin[strain_splits[i]])

# take chunks of stress and resistance relaxing values from index i1 to i2
i1 = []
i2 = []
# i1 = [1,7,11,15,21] # MANUALLY ADD VALUES
# i2 = [2,8,12,16,22]
# i3 = [5,9,13,19,25]
# i4 = [6,10,14,20,26]
for i in range(int(len(strain_splits)/4)):
    i1.append(int(4*i + 1))
    i2.append(int(4*i + 2))

## Plotting and fitting the loading and unloading of the different specimens

Strain_load_arr = []
Stress_load_arr = []

def lin_func(x,m,c):
    return m*x + c

consts_lin_fit = []
p_init = [1,1]

constsu_lin_fit = []
pu_init = [1,1]

try:
    fig3, ax3 = plt.subplots(figsize=(12,6))
    # for i in range(1,int(len(strain_splits))):
    for i in range(5):
        Strain_load = Strain_t[int(strain_splits[4*i+1]):int(strain_splits[4*i+2])]
        Strain_unload = Strain_t[int(strain_splits[4*i+3]):int(strain_splits[4*i+4])]

        Stress_load = Stress_t[int(strain_splits[4*i+1]):int(strain_splits[4*i+2])]
        # Stress_load_min = np.min(Stress_load)
        # Stress_load = Stress_load - Stress_load_min
        Stress_unload = Stress_t[int(strain_splits[4*i+3]):int(strain_splits[4*i+4])]
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
    

finally:

    plt.show()
