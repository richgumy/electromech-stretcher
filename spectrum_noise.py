# -*- coding: utf-8 -*-
"""
FILE: spectrum_noise.py
AUTHOR: R Ellingham
DATE MODIFIED: Jun 2021
PROGRAM DESC: Process stationary sensor data from a conductive elastomer 
DATE CREATED: Oct 2020

Parameters measured:
-> Resistance, Force, and their respective time stamps (for each individual measurement)

TODO:
1) Hypothesize the type of noise seen
2) Determine validity of data and analysis

"""

import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize
import scipy.signal as sig
import scipy.stats as stats
import matplotlib

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
    """
    DESCR: Remove any erroneously recorded zero values
    IN_PARAMS: x data array
    NOTES: Requires numpy library.
    TODO:
    """
    new_strain_arr = strain_arr
    for i in range(1,len(strain_arr)):
        if abs(strain_arr[i] - strain_arr[i-1]) > 1:
            new_strain_arr[i] = strain_arr[i-1]
    return new_strain_arr

def PSD(x, fs):
    """Estimate PSD of signal x with sampling frequency fs"""
    
    # Estimate PSD using Welch's method
    return sig.welch(x, fs)

def ASD(x, fs):
    """Estimate ASD of signal x with sampling frequency fs"""
    
    f, Px = PSD(x, fs)

    return f, np.sqrt(Px)

# Test specimen dimensions (in m):
spec_length = 40e-3
spec_length_inner_electrodes = 20e-3
spec_width = 10e-3
spec_thickness = 4e-3

# input_filename = "2_7-5_Epin_20mm_v4.csv" # AC
input_filename = "2_7-5_Epin_20mm_v5.csv" # DC

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

# Determine the mean and noise of resistance
Ri_zero = Ri - np.average(Ri)
Ro_zero = Ro - np.average(Ro)

# Plot R measurements over time
fig1, axs1 = plt.subplots(2, 1, constrained_layout=True,figsize = (10, 10))

ax = axs1[0]
ax.plot(tR, Ri_zero,'r-')
ax.set_title('')
ax.set_ylabel('Resistance i [Ohm]')

ax = axs1[1]
ax.plot(tR, Ro_zero,'r-')
ax.set_title('')
ax.set_ylabel('Resistance o [Ohm]')
ax.grid(True)


# Histogram plot of noise of Ri and Ro
bins = 200
plt.figure()
plt.hist(Ri_zero,bins)
plt.title("Ri noise distrib")
plt.figure()
plt.hist(Ro_zero,bins)
plt.title("Ro noise distrib")

# Split Ri to remove hacky AC artifact
# Ri_zero1 = []
# Ri_zero2 = []
# for i in range(int(len(Ri_zero)/2)):
#     Ri_zero1.append(Ri_zero[2*i])
#     Ri_zero2.append(Ri_zero[2*i+1])
    
# Ri_mean1 = np.average(Ri_zero1)
# Ri_mean2 = np.average(Ri_zero2)
    
# for i in range(len(Ri_zero1)):
#     Ri_zero1[i] = Ri_zero1[i] - Ri_mean1
#     Ri_zero2[i] = Ri_zero2[i] - Ri_mean2
    
# plt.figure()
# plt.hist(Ri_zero1,bins)
# plt.hist(Ri_zero2,bins)
# plt.figure()
# plt.hist(Ri_zero2,bins)
# plt.hist(Ri_zero1,bins)

#%% 

# Interpolate the Ri and Ro data so that there is a constant sample frequency
Rointp_P = np.interp(tP,tR,Ro)
Rointp_F = np.interp(tF,tR,Ro)

Riintp_P = np.interp(tP,tR,Ri)
Riintp_F = np.interp(tF,tR,Ri)

Ro_tot = np.append(Rointp_F,[Rointp_P,Ro])
Ri_tot = np.append(Riintp_F,[Riintp_P,Ri])

tR_tot = np.append(tF,[tP,tR])

fs_R = tR_tot[1] - tR_tot[0]

# Determine the mean and noise of interpolated resistance
Ri_zero_int = Ri_tot
Ro_zero_int = Ro_tot

# Plot R measurements over time
fig1, axs1 = plt.subplots(2, 1, constrained_layout=True,figsize = (10, 10))

ax = axs1[0]
ax.plot(tR_tot, Ri_zero_int,'r.')
ax.set_title('Interpolated Electrode Resistance Data')
ax.set_ylabel('Resistance i [Ohm]')

ax = axs1[1]
ax.plot(tR_tot, Ro_zero_int,'r.')
ax.set_title('')
ax.set_ylabel('Resistance o [Ohm]')
ax.grid(True)


# Histogram plot of noise of Ri and Ro
bins = 200
plt.figure()
plt.hist(Ri_zero,bins)
plt.title("$Ri_{interp}$ noise distrib")
plt.figure()
plt.hist(Ro_zero,bins)
plt.title("$Ro_{interp}$ noise distrib")

# Estimate the amplitude spectral density of the R signals
f, A_Ro = ASD(Ro_tot, fs_R)
f, A_Ri = ASD(Ri_tot, fs_R)

fig, axes = plt.subplots(1)
axes.plot(f, A_Ro, label='$A_{Ro}$')
axes.set_title('Resistance ASD')
axes.set_xlabel('Frequency [Hz]')
axes.set_ylabel('Resistance o ASD [Ohm]')
axes.legend(loc='upper right')
axes.grid(True)

axes2 = axes.twinx()
axes2.plot(f, A_Ri, label='$A_{Ri}$', color='r')
axes2.legend(loc='right')
axes2.set_ylabel('Resistance i ASD [Ohm]')
plt.show()


