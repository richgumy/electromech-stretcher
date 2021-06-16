"""
FILE: post_processing.py
AUTHOR: R Ellingham
DATE: June 2021
PROGRAM DESC: 
    Convert: position data into strain %, and force data into stress.
    Filter: Outliers in the strain data, AC input effects on resistance, noise seen in stress
    Plot: Resistance, strain, stress, in time.

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
    DESCR: Smooths any erroneous zero values in data (i.e. strain measurement drops 
                                                      to zero erroneously and aperiodically)
    IN_PARAMS: strain array
    NOTES:
    TODO:    
    """
    new_strain_arr = strain_arr
    for i in range(1,len(strain_arr)):
        if abs(strain_arr[i] - strain_arr[i-1]) > 3:
            new_strain_arr[i] = strain_arr[i-1]
    return new_strain_arr

font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

# Test specimen dimensions (in m):
spec_length = 40e-3
spec_length_inner_electrodes = 20e-3
spec_width = 10e-3
spec_thickness = 4e-3

input_filename = ".csv"

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
            

P = pos_outlier_corrector(P) * 1e-3 # Correct any erroneously gathered position values and scale to meters

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
Ri_t = -np.interp(t_lin,t_,Ri_)

# Use a MAF on the resistance data if it is an AC measurement
Ri_t = MAF(Ri_t,12)

