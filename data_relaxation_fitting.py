"""
FILE: data_relaxation_fitting.py
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
    flat_flag = 0
    n = len(data)
    for i in range(1,n-1):
        if flat_flag == 0 and (data[i-1] != data[i] and data[i] == data[i+1]) or (data[i-1] == data[i] and data[i] != data[i+1]):
            index_splits.append(i)
    index_splits.append(n-1)
    return index_splits


def main(input_filename):
    R = np.array([])
    tR = np.array([])
    P = np.array([])
    tP = np.array([])
    F = np.array([])
    tF = np.array([])

    # Write to csv
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
    F_tot = np.append(F,[Fintp_P,Fintp_P])
    tF_tot = np.append(tF,[tP,tR])

    P_tot = np.append(P,[Pintp_F,Pintp_R])
    tP_tot = np.append(tP,[tF,tR])

    R_tot = np.append(R,[Rintp_P,Rintp_F])
    tR_tot = np.append(tR,[tP,tF])

    # Calc stress and strain from F_tot and P_tot
    Stress = F/(10e-3*4e-3)
    Stress_tot = F_tot/(10e-3*4e-3) # force/cross-sectional-area
    Strain = -P/30
    Strain_tot = -P_tot/30 # dx/x

    # Plot measurements over time
    fig1, axs1 = plt.subplots(3, 1, constrained_layout=True)

    ax = axs1[0]
    ax.plot(tR_tot, R_tot,'rx')
    ax.set_title('')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs1[1]
    ax.plot(tP_tot, Strain_tot,'rx')
    ax.set_title('')
    ax.set_ylabel('Strain')
    ax.grid(True)

    ax = axs1[2]
    ax.plot(tF_tot, Stress_tot,'rx')
    ax.set_title('F')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Stress [Pa]')
    ax.grid(True)

    ## Plot Res vs XX measurements
    fig3, axs3 = plt.subplots(2, 1, constrained_layout=True)

    ax = axs3[0]
    ax.plot(Strain_tot, R_tot,'r-')
    ax.set_title('')
    ax.set_xlabel('Strain')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs3[1]
    ax.plot(Stress_tot, R_tot,'ro')
    ax.set_title('')
    ax.set_xlabel('Stress[Pa]')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ## Plot Res vs strain (loading and unloading) measurements (specific for first_test_num12.csv)
    strain_splits = split_ramp_data(Strain_tot)
    print(strain_splits)
    print(len(strain_splits))
    # take a chunk of relaxing values from index i1 to i2
    i1 = [1,5,9,13]
    i2 = [2,6,10,14]
    for i in range(len(i1)):
        R_load = R_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Strain_load = Strain[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Stress_load = Stress_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        t_load = tR_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]
        t_load = t_load - t_load[0]
        # Curve fitting code (curve_fit func using non-lin lstsqr)
        def f(t, a, b, c, d):
            return a * np.exp(-b * (t-c)) + d
        popt, pcov = optimize.curve_fit(f, t_load, Stress_load)

        print("Formula:%.2f * exp(%.2f*(t-%.2f)) + %.2f" % (popt[0],popt[1],popt[2],popt[3]))

        t_load_lin = np.linspace(min(t_load),max(t_load) , 20)
        Stress_load_lin = f(t_load_lin,*popt)

        plt.figure()
        ax3 = plt.plot(t_load,Stress_load,'x',t_load_lin,Stress_load_lin,'-')
        plt.xlabel('t')
        plt.ylabel('stress')

    plt.show()

# The real main driver
if __name__ == "__main__":
	# Default input parameters
	input_filename = "Preliminary tests/first_test_num12.csv"
    # input2
	import sys
	if len(sys.argv)>1: input_filename   = (sys.argv[1])
	# if len(sys.argv)>2: input2   =int(sys.argv[2])

	main(input_filename)
