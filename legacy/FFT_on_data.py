"""
FILE: FFT_on_data.py
AUTHOR: R Ellingham
DATE: Jan 2021
PROGRAM DESC: Investigation into whether there are resonant frequencies within
the data gathered. In particular determining noise within the load cell readings

TODO:
1)
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

def main(input_filename, spec_length=40e-3, spec_width=10e-3, spec_thickness=4e-3):
    """
    DESCR: Completes data manipulation
    IN_PARAMS: filename for raw csv data, specimen dimensions in m
    NOTES: Not well tested.
    TODO:
    """
    R = np.array([])
    tR = np.array([])
    P = np.array([])
    tP = np.array([])
    F = np.array([])
    tF = np.array([])

    # Extract from csv
    if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input
    with open(input_filename, 'r', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        line_count = 0
        for row in data:
            if line_count == 0:
                # print(f'Column names are {", ".join(row)}')
                line_count = 1
            else:
                R = np.append(R,float(row[0]))
                tR = np.append(tR,float(row[1]))
                P = np.append(P,float(row[2]))
                tP = np.append(tP,float(row[3]))
                F = np.append(F,float(row[4]))
                tF = np.append(tF,float(row[5]))
                # print(row[0],row[1],row[2],row[3],row[4],row[5])

    # Interpolate all data
    Fintp_P = np.interp(tP,tF,F)
    Fintp_R = np.interp(tR,tF,F)

    Pintp_F = np.interp(tF,tP,P)
    Pintp_R = np.interp(tR,tP,P)

    Rintp_P = np.interp(tP,tR,R)
    Rintp_F = np.interp(tF,tR,R)

    # Generate new interpolated data for the '_tot' arrays
    # interpolation order is [force,pos,res]
    F_tot = np.append(F,[Fintp_P,Fintp_R])
    tF_tot = np.append(tF,[tP,tR])

    P_tot = np.append(Pintp_F,[P,Pintp_R])
    tP_tot = np.append(tF,[tP,tR])

    R_tot = np.append(Rintp_F,[Rintp_P,R])
    tR_tot = np.append(tF,[tP,tR])

    ## Sort all _tot data into time order
    # Make new arrays for sorted data
    F_ = np.zeros(len(tF_tot))
    P_ = np.zeros(len(tF_tot))
    R_ = np.zeros(len(tF_tot))
    t_ = np.zeros(len(tF_tot))

    ind = np.argsort(tF_tot)
    for i in range(len(tF_tot)):
        F_[i] = F_tot[ind[i]]
        P_[i] = P_tot[ind[i]]
        R_[i] = R_tot[ind[i]]
        t_[i] = tF_tot[ind[i]]

    # Interpolate again to get time in a linear form
    t_lin = np.linspace(min(t_),max(t_),int(len(t_)/3)) # divide by 3 to get original size of data
    F_t = np.interp(t_lin,t_,F_)
    P_t = np.interp(t_lin,t_,P_)
    R_t = np.interp(t_lin,t_,R_)

    print(tF)
    L = len(tF) # length of the data of significance
    T_s = np.mean(np.diff(tF)) # average sample period
    print(T_s)
    F_s = 1/T_s # averge sample freq
    F_n = F_s/2 # Nyquist freq
    F_v = np.linspace(0, 1, L)*F_n # Freq vector
    F_fft = np.fft.fft(F)/L # FFT
    F_fft_mag = np.sqrt(F_fft.real**2 + F_fft.imag**2)

    plt.plot(np.ones(len(tF)), tF, 'x')
    plt.show()

    return 0




# The real main driver
if __name__ == "__main__":
	# Default input parameters
	input_filename = "Preliminary tests/first_test_num12.csv"
    # input2
	import sys
	if len(sys.argv)>1: input_filename   = (sys.argv[1])
	# if len(sys.argv)>2: input2   =int(sys.argv[2])

	main(input_filename)
