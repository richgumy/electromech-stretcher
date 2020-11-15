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
from matplotlib.ticker import MaxNLocator
import numpy as np

def main(input_filename):
    R = np.array([])
    tR = np.array([])
    P = np.array([])
    tP = np.array([])
    F = np.array([])
    tF = np.array([])

    # Write to csv
    with open(input_filename + '.csv', 'r', newline='') as csvfile:
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
    ax.plot(tR_tot, R_tot,'ro')#,tR, R,'x')
    ax.set_title('')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs1[1]
    ax.plot(tP_tot, Strain_tot,'ro')#, tP, Strain, 'x')
    ax.set_title('')
    ax.set_ylabel('Strain')
    ax.grid(True)

    ax = axs1[2]
    ax.plot(tF_tot, Stress_tot,'ro')#, tF, Stress, 'x')
    ax.set_title('F')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Stress [Pa]')
    ax.grid(True)

    # Plot Res vs XX measurements
    fig2, axs2 = plt.subplots(2, 1, constrained_layout=True)

    ax = axs2[0]
    ax.plot(Strain_tot, R_tot,'ro')
    ax.set_title('')
    ax.set_xlabel('Strain')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    ax = axs2[1]
    ax.plot(Stress_tot, R_tot,'ro')
    ax.set_title('')
    ax.set_xlabel('Stress[Pa]')
    ax.set_ylabel('Resistance [Ohm]')
    ax.grid(True)

    # print(Stress)
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
	input_filename = "1st_test_num10"
    # input2
	import sys
	if len(sys.argv)>1: input_filename   = (sys.argv[1])
	# if len(sys.argv)>2: input2   =int(sys.argv[2])

	main(input_filename)
