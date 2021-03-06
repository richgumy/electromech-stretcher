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
    n = len(data)
    for i in range(1,n-1):
        if (data[i-1] != data[i] and data[i] == data[i+1]) or (data[i-1] == data[i] and data[i] != data[i+1]):
            index_splits.append(i)
    index_splits.append(n-1)
    return index_splits

def MAF(x,dx):
    """
    DESCR: Moving average filter for a set of data averaging +/- each point (Bar
    the first and last 'dx' samples)
    IN_PARAMS: x data, sample period to average over
    NOTES:
    TODO:
    """
    xn = []
    x_sum = 0
    for i in range(dx,len(x)-dx):
        x_sum = x[i]
        for j in range(1,dx+1):
            x_sum = x_sum + x[i-j] + x[i+j]
        xn.append(x_sum/(2*dx+1))
    return xn


def main(input_filename, spec_length=40e-3, spec_width=10e-3, spec_thickness=4e-3):
    """
    DESCR: Completes data manipulation
    IN_PARAMS: filename for raw csv data, specimen dimensions in m
    NOTES: Not well tested.
    TODO:
    """
    Rout = np.array([])
    Rin = np.array([])
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
                Rout = np.append(R,float(row[0]))
                Rin = np.append(R,float(row[1]))
                tR = np.append(tR,float(row[2]))
                P = np.append(P,float(row[3]))
                tP = np.append(tP,float(row[4]))
                F = np.append(F,float(row[5]))
                tF = np.append(tF,float(row[6]))
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

    # Calc strain from displacement
    # Strain_tot = -P_tot/(spec_length*1000) # dx/x
    Strain_tot = -P_t/(spec_length*1000) # dx/x
    # Calc stress from force and changing strain
    possion_ratio = 0.25 # Poisson's ratio
    # Stress_tot = F_tot/((spec_width*spec_thickness)*(-Strain_tot*possion_ratio+1)*(-Strain_tot*possion_ratio+1)) # force/cross-sectional-area
    Stress_tot = F_t/((spec_width*spec_thickness)*(-Strain_tot*possion_ratio+1)*(-Strain_tot*possion_ratio+1)) # force/cross-sectional-area
    # # Calc stress from force (neglecting the changing cross-sectional area)
    # Stress = F/(spec_width*spec_thickness)
    # Stress_tot = F_tot/(spec_width*spec_thickness) # force/cross-sectional-area

    R_tot = R_t

    # Apply moving average filter to stress data
    filter_num = 5
    stress_fil = MAF(Stress_tot,filter_num)
    res_fil = MAF(R_tot,filter_num)
    strain_fil = MAF(Strain_tot,filter_num)
    t_fil = t_lin[filter_num:-filter_num]

    # Plot measurements over time
    fig1, axs1 = plt.subplots(3, 1, constrained_layout=True)

    ax = axs1[0]
    ax.plot(t_lin, R_tot,'rx',t_fil, res_fil,'b-')
    ax.set_title('')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs1[1]
    ax.plot(t_lin, Strain_tot,'rx',t_fil, strain_fil,'b-')
    ax.set_title('')
    ax.set_ylabel('Strain')
    ax.grid(True)

    ax = axs1[2]
    ax.plot(t_lin, Stress_tot,'rx',t_fil, stress_fil,'b-')
    ax.set_title('F')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Stress [Pa]')
    ax.grid(True)


    ## Plot relaxation and fit curve
    # Split data into piece-wise data
    strain_splits = split_ramp_data(Strain_tot)

    # take a chunk of relaxing values from index i1 to i2
    i1 = [1,5,9,13]
    i2 = [2,6,10,14]

    # Curve fitting code (curve_fit func using non-lin lstsqr)
    # stress(t) = Y * strain * exp^(-(Y/mu)*t)
    def f(t, a, b, c, d):
        return a * np.exp(-b * (t-c)) + d

    for i in range(4):
        Res_load = R_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Strain_load = Strain_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]

        Stress_load = Stress_tot[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]
        Stress_load_fil = stress_fil[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])-6] # MAF data. Leads other data by 6 samples ()

        t_load = t_lin[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])]
        t_load = t_load - t_load[0]
        t_load_fil = t_lin[int(strain_splits[i1[i]]):int(strain_splits[i2[i]])-6] # match size of MAF data
        t_load_fil = t_load_fil - t_load_fil[0]

        # Levenberg–Marquardt algorithm for non-linear leastsq
        poptS, pcovS = optimize.curve_fit(f, t_load, Stress_load)
        poptS_fil, pcovS_fil = optimize.curve_fit(f, t_load_fil, Stress_load_fil,maxfev=50000)
        poptR, pcovR = optimize.curve_fit(f, t_load, Res_load)

        print("Stress Formula:%.2f * exp(%.2f*(t-%.2f)) + %.2f" % (poptS[0],poptS[1],poptS[2],poptS[3]))
        print("Stress Formula(filtered):%.2f * exp(%.2f*(t-%.2f)) + %.2f" % (poptS_fil[0],poptS_fil[1],poptS_fil[2],poptS_fil[3]))
        print("Resistance Formula:%.2f * exp(%.2f*(t-%.2f)) + %.2f" % (poptR[0],poptR[1],poptR[2],poptR[3]))
        # print("Covar:")
        # print(pcov)

        t_load_lin = np.linspace(min(t_load),max(t_load) , 20)
        Stress_load_lin = f(t_load_lin,*poptS)
        Stress_load_fil_lin = f(t_load_lin,*poptS_fil)
        Res_load_lin = f(t_load_lin,*poptR)

        fig, axs3 = plt.subplots(3, 1, constrained_layout=True)

        ax = axs3[0]
        ax.plot(t_load,Stress_load,'bx',t_load_lin,Stress_load_lin,'y-')
        ax.set_title('')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Stress [Pa]')
        ax.grid(True)

        ax = axs3[1]
        ax.plot(t_load_fil,Stress_load_fil,'bx',t_load_lin,Stress_load_fil_lin,'y-')
        ax.set_title('')
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Stress [Pa]')
        ax.grid(True)

        ax = axs3[2]
        ax.plot(t_load,Res_load,'rx',t_load_lin,Res_load_lin,'y-')
        ax.set_title('')
        ax.set_xlabel('Time[s]')
        ax.set_ylabel('Resistance [Ohm]')
        ax.grid(True)


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
