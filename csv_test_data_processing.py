"""
FILE: csv_test_data_processing.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: Process electromechanical data measured from a stretched
conductive elastomer in real time.

TODO:
1) Make functions for script
2)
"""
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize

def diff_data(x,t):
    """
    DESCR: Differentiate 'x' w.r.t. 't' to get output list 'dx/dt'
    IN_PARAMS: x data, t data
    NOTES:
    TODO:
    """
    dx_dt = []
    dx_dt.append((x[1]-x[0])/(t[1]-t[0])) # Slightly erroneous first value to maintain len(x) == len(dx_dt)
    for i in range(1,len(x)-1):
        dx_dt.append(((x[i]-x[i-1])/(t[i]-t[i-1])+(x[i+1]-x[i])/(t[i+1]-t[i]))/2)
    dx_dt.append((x[len(x)-1]-x[len(x)-2])/(t[len(x)-1]-t[len(x)-2])) # Slightly erroneous final value to maintain len(x) == len(dx_dt)
    return dx_dt

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
        print(i)
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

    # Differentiate displacement w.r.t. time
    V_tot = diff_data(P_tot,tP_tot)
    A_tot = diff_data(V_tot,tP_tot)

    # Plot measurements over time
    fig1, axs1 = plt.subplots(3, 1, constrained_layout=True)

    ax = axs1[0]
    ax.plot(tR_tot, R_tot,'r-')
    ax.set_title('')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs1[1]
    ax.plot(tP_tot, Strain_tot,'r-')
    ax.set_title('')
    ax.set_ylabel('Strain')
    ax.grid(True)

    ax = axs1[2]
    ax.plot(tF_tot, Stress_tot,'r-')
    ax.set_title('F')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Stress [Pa]')
    ax.grid(True)

    # ax = axs1[3]
    # ax.plot(tP_tot, V_tot,'r-')
    # ax.set_title('F')
    # ax.set_xlabel('Time [s]')
    # ax.set_ylabel('Velocity [mm/s]')
    # ax.grid(True)
    #
    # ax = axs1[4]
    # ax.plot(tP_tot, A_tot,'r-')
    # ax.set_title('F')
    # ax.set_xlabel('Time [s]')
    # ax.set_ylabel('Acceleration [mm/s]')
    # ax.grid(True)
    #
    # # Overlap Res and strain plots in time
    # fig2, axs2 = plt.subplots()
    #
    # axs2.plot(tR_tot, R_tot,'r-')
    # axs2.set_ylabel('Resistance [Ohm]')
    # axs2.grid(True)
    #
    # axs2a = axs2.twinx()
    # axs2a.set_navigate(False)
    # axs2a.plot(tP_tot, Strain_tot,'b-')
    # axs2a.set_ylabel('Strain')
    # axs2a.grid(True)

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
    R_load = R_tot[0:111]
    # R_load = np.concatenate((R_tot[0:111],R_tot[221:288],R_tot[355:408]))
    # R_load = np.concatenate((R_load,R_tot[461:508]))
    Strain_load = Strain[0:111]

    # Strain_load = np.concatenate((Strain[0:111],Strain[221:288],Strain[355:408]))
    # Strain_load = np.concatenate((Strain_load,Strain[461:508]))
    Stress_load = Stress_tot[0:111]

    # R_unload = np.concatenate(R_tot[111:221],R_tot[288:355],R_tot[408:461],R_tot[508:555])
    # Strain_unload = np.concatenate(Strain[111:221],Strain[288:355],Strain[408:461],Strain[508:555])

    # Curve fitting code (curve_fit func using non-lin lstsqr)
    def f(x, a, b, c):
        return a * np.exp(b * x) + c
    popt, pcov = optimize.curve_fit(f, Strain_load, R_load)

    # Calc R squared value here...
    SS_res = 0
    # SS_tot = 0
    for i in range(len(R_load)):
        SS_res = SS_res + (R_load[i]-f(Strain_load[i], *popt))**2

    print("Formula:%.2f * exp(%.2f*(Strain)) + %.2f" % (popt[0],popt[1],popt[2]))
    print("SS_res = %.2f" % (SS_res))

    # def f(x, b, c, t):
    #     return b**(t*x)+c

    # def f(x, a, b, c):
    #     return a*np.exp(b*x)+c
    #
    # def residual(p, x, y):
    #     return y - f(x, *p)
    #
    # p0 = [1.0,1.0]
    #
    # # popt, pcov = optimize.leastsq(residual, p0, args=(Strain_load, R_load))
    # popt, pcov = optimize.curve_fit(f, Strain_load, R_load)
    # # Determine the R_square value between 0 and 1. 1 is a strong correlation
    # # R_sqr_load = 1 - pcov/(R_load.size*R_load.var())
    #
    # # print("R = %.4f*(Strain)^6 + %.4f; R_sqr = %.4f" % (popt[0], popt[1], R_sqr_load))
    #
    Strain_load_lin = np.linspace(min(Strain_load),max(Strain_load) , 20)
    Res_lin = f(Strain_load_lin,*popt)

    plt.figure()
    ax3 = plt.plot(Strain_load,R_load,'x',Strain_load_lin,Res_lin,'-')
    plt.xlabel('Strain')
    plt.ylabel('Res [Ohm]')



    # 3D surface plot
    fig3d = plt.figure()
    ax3d = fig3d.add_subplot(111, projection='3d')
    ax3d.plot_trisurf(np.log(tP_tot), Stress_tot, Strain_tot)
    #
    # fig3da = plt.figure()
    # ax3da = fig3da.add_subplot(111, projection='3d')
    # ax3da.scatter(V_tot, Strain_tot, R_tot)

    # Use linear least squares to find Young's modulus -> Stress = Y * Strain + offset_error
    A = np.vstack([Strain_tot,np.ones(len(Strain_tot))]).T
    model = np.linalg.lstsq(A, Stress_tot, rcond=None)
    Y, offset_error= model[0]
    resid = model[1]
    # Determine the R_square value between 0 and 1. 1 is a strong correlation
    R_sqr = 1 - resid/(Stress_tot.size*Stress_tot.var())
    print("Y = %.4f, offset_error = %.4f, R_sqr = %.4f" % (Y, offset_error, R_sqr))

    Strain_lin = np.linspace(min(Strain_tot),max(Strain_tot) , 5)
    Stress_lin = Y * Strain_lin + offset_error

    plt.figure()
    ax3 = plt.plot(Strain_tot,Stress_tot,'x',Strain_lin,Stress_lin,'-')
    plt.xlabel('Strain')
    plt.ylabel('Stress [Pa]')

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
