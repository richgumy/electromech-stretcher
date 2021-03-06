"""
FILE: data_relaxation_fitting_spyder.py
AUTHOR: R Ellingham
DATE: Dec 2020
PROGRAM DESC: Fit a relaxation model to stress relaxation data.

TODO:
1)
2)
"""

import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize

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
    NOTES: Requires numpy library. Use an odd numbered N size window.
    TODO:
    """
    x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
    x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')

    return x_smooth

def E1_E2_solver(x,y):
    """
    DESCR: Solves quadratic and determines realistic elastic moduli values for
    SLM (solid linear model).
    From equation:
       stress(t) =  x*e0 + C*np.exp(-(y/mu)*t)
    converting to ...
       stress(t) = (E1*E2)/(E1+E2) * e0 + C*np.exp(-((E1+E2)/mu)*t)
    IN_PARAMS: x,y
    NOTES:
    TODO:
    """
    # Finding actual E1 and E2
    E1 = 0
    E2 = 0
    E2_p = (y+np.sqrt(y**2-4*x*y))/2
    E2_n = (y-np.sqrt(y**2-4*x*y))/2
    if E2_p > 0 and E2_n < 0:
        E2 = E2_p
        E1 = E2 - y
    elif E2_p > 0 and E2_n >= 0:
        E1_n = y - E2_n
        E1_p = y - E2_p
        if E1_p > 0 and E1_n < 0:
            E1 = E1_p
            E2 = E2_p
        else:
            E1 = E1_n
            E2 = E2_n
            print("Other positive solutions for E1 & E2 may exist.")
    if E1 == 0 or E2 == 0:
        E1 = 1
        E2 = 1
        print("E1 or E2 equal zero. Non SLM model!")
    return E1, E2

def find_Rsqr(f, x_data, y_data):
    """
    DESCR: Gives the Rsquared value for x/y data fitted to function 'f(x)=y'
    IN_PARAMS: f(x_data), x data, y data
    NOTES:
    TODO: Alter for use with non-linear function 'f'
    """
    SS_tot = 0
    SS_resid = 0
    y_avg = sum(y_data)/len(y_data)
    for i in range(len(x_data)):
        SS_tot = SS_tot + (y_data[i] - y_avg)**2
        SS_resid = SS_resid + (y_data[i] - f[i])**2
    R_sqr = 1 - SS_resid/SS_tot
    return R_sqr



# Test specimen dimensions:
spec_length = 40e-3
spec_width = 10e-3
spec_thickness = 4e-3

input_filename = "2_7-5_E4pin_20mm_v13_0.1_0.2_0.3_strain.csv"

R = np.array([])
Ri = np.array([])
tR = np.array([])
P = np.array([])
tP = np.array([])
F = np.array([])
tF = np.array([])

# Extract from csv
if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input
with open(input_filename, 'r', newline='') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    data = list(data)
    line_count = 0
    i = 0
    for row in data:
        if line_count == 0:
            # print(f'Column names are {", ".join(row)}')
            line_count = 1
        else:
            R = np.append(R,float(row[0]))
            Ri = np.append(Ri,float(row[1]))
            tR = np.append(tR,float(row[2]))
            P = np.append(P,float(row[3]))
            tP = np.append(tP,float(row[4]))
            F = np.append(F,float(row[5]))
            tF = np.append(tF,float(row[6]))
            print("%.2f" % (i/len(data)))
            i = i + 1
            # print(row[0],row[1],row[2],row[3],row[4],row[5])

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

# Calc engineering strain from displacement
Strain_t = -P_t/(spec_length*1000) #dx/x

# Calc actual strain
Strain_true_t = np.log((-P_t/1000+spec_length)/(spec_length))

# Calc actual stress from force and changing strain
poisson_ratio = 0.29 # Poisson's ratio. Found experimentally using #2_7.5%dragonskin10NV specimen
Stress_t = F_t/((spec_width*spec_thickness)*(-Strain_t*poisson_ratio+1)*(-Strain_t*poisson_ratio+1)) # force/cross-sectional-area
Stress_true_t = (F_t/spec_width*spec_thickness) * (1+Strain_t)

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
ax.plot(t_lin, Stress_true_t,'r-')
ax.set_title('')
ax.set_ylabel('True stress [Ohm]')
ax.grid(True)

ax = axs1[1]
ax.plot(t_lin, Convert,'r-')
ax.set_title('')
ax.set_ylabel('Eng stress [Ohm]')
ax.grid(True)

ax = axs1[2]
ax.plot(t_lin, Strain_true_t,'r-')
ax.set_title('')
ax.set_ylabel('True strain')
ax.grid(True)

ax = axs1[3]
ax.plot(t_lin, Strain_t,'r-')
ax.set_title('')
ax.set_ylabel('Eng strain')
ax.grid(True)


## Plot relaxation and fit curve
# Split data into piece-wise data
strain_splits = split_ramp_data(Strain_t)

# take chunks of stress and resistance relaxing values from index i1 to i2
i1 = []
i2 = []
for i in range(int(len(strain_splits)/4)):
    i1.append(int(4*i + 1))
    i2.append(int(4*i + 2))


# Curve fitting code (curve_fit func using non-lin lstsqr)
e0 = .10
def SLS_relax(t,E1,E2,C,mu):
    return (E1*E2)/(E1+E2) * e0 + C*np.exp(-((E1+E2)/mu)*t)

def SLS_relax_simple(t,E1,E2,C,mu):
    return E1 * e0 + C*np.exp(-(E2/mu)*t)

# def f(t, a, b, c, d):
#     return a * np.exp(-b * (t-c)) + d

def generalisedi3(x, a0, a1, tau_1, a2, tau_2): # generalised Kelvin SLS relaxation model for n = 3
    return a0 + a1 * np.exp(-x/tau_1) + a2 * np.exp(-x/tau_2)

def generalisedi2(x, a0, a1, tau_1): # generalised Kelvin SLS relaxation model for n = 3
    return a0 + a1 * np.exp(-x/tau_1)

consts_SLSrelax = []
consts_geni2 = []
pGen_init = [1,1,1]

try:
    for i in range(int(len(strain_splits)/4)):
        Res_load = Ri_t[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Strain_load = Strain_t[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Stress_load = Stress_t[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]
        Stress_load_min = np.min(Stress_load)
        Stress_load = Stress_load - Stress_load_min

        t_load = t_lin[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]
        t_load = t_load - t_load[0]


        # Levenberg–Marquardt algorithm for non-linear leastsq, fitting the stress data to SLS relaxation models
        poptS, pcovS = optimize.curve_fit(SLS_relax, t_load, Stress_load, maxfev=50000)
        print("SLS_relax(t,E1,E2,C,mu):",poptS)
        poptS_simple, pcovS_simple = optimize.curve_fit(SLS_relax_simple, t_load, Stress_load, maxfev=50000)
        poptS_gen, pcovS_gen = optimize.curve_fit(generalisedi2, t_load, Stress_load, p0=pGen_init, maxfev=50000)
        print("generalisedi2(x, a0, a1, tau_1):",poptS_gen)

        pGen_init = poptS_gen

        t_load_lin = np.linspace(min(t_load),max(t_load) , 100)
        Stress_load_lin = SLS_relax(t_load_lin,*poptS)
        Stress_load_lin_simple = SLS_relax_simple(t_load_lin,*poptS_simple)
        Stress_load_lin_f = generalisedi2(t_load_lin,*poptS_gen)

        fig, axs = plt.subplots(1, 1, constrained_layout=True)

        ax = axs
        ax.plot(t_load,Stress_load+Stress_load_min,'bx', t_load_lin,
                Stress_load_lin_f+Stress_load_min, 'r-', t_load_lin,
                Stress_load_lin+Stress_load_min, 'y-', t_load_lin,
                Stress_load_lin_simple+Stress_load_min,'m-')
        ax.set_title('')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Stress [Pa]')
        ax.legend(["Data", "Generalised SLM i=3", "SLS fit"])
        ax.grid(True)


        # ax = axs3[1]
        # ax.plot(t_load,Stress_load,'rx',t_load_lin,Stress_load_lin_simple,'y-')
        # ax.set_title('')
        # ax.set_xlabel('Time[s]')
        # ax.set_ylabel('Resistance [Ohm]')
        # ax.grid(True)

        poptS[2] = poptS[0] + Stress_load_min
        poptS_gen[2] = poptS_gen[0] + Stress_load_min
        consts_SLSrelax.append(poptS)
        consts_geni2.append(poptS_gen)

finally:
    print(consts_geni2)

    for j in (range(len(consts_geni2[0]))):
        const = []
        for i in (range(len(consts_geni2))):
            const.append(consts_geni2[i][j])
        fig, axs = plt.subplots(1, 1, constrained_layout=True)
        ax = axs
        ax.plot(const,'bx')
        ax.grid(True)

    plt.show()
