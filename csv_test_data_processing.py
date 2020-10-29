"""
FILE: csv_test_data_processing.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: Process electromechanical data measured from a stretched
conductive elastomer in real time.

TODO:
1) Make functions
"""
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import numpy as np

R = np.array([])
tR = np.array([])
P = np.array([])
tP = np.array([])
F = np.array([])
tF = np.array([])

# Write to csv
with open('testy2.csv', 'r', newline='') as csvfile:
    data = csv.reader(csvfile, delimiter=',')
    line_count = 0
    for row in data:
        if line_count == 0:
            print(f'Column names are {", ".join(row)}')
            line_count = 1
        else:
            R = np.append(R,float(row[0]))
            tR = np.append(tR,float(row[1]))
            P = np.append(P,float(row[2]))
            tP = np.append(tP,float(row[3]))
            F = np.append(F,float(row[4]))
            tF = np.append(tF,float(row[5]))
            print(row[0],row[1],row[2],row[3],row[4],row[5])

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

# Plot measurements over time
fig1, axs1 = plt.subplots(3, 1, constrained_layout=True)

ax = axs1[0]
ax.plot(tR_tot, R_tot,'ro',tR, R,'x')
ax.set_title('R')
ax.grid(True)

ax = axs1[1]
ax.plot(tP_tot, P_tot,'ro', tP, P, 'x')
ax.set_title('P')
ax.grid(True)

ax = axs1[2]
ax.plot(tF_tot, F_tot,'ro', tF, F, 'x')
ax.set_title('F')
ax.grid(True)

# Plot stress over strain
Stress = F_tot/(10e-3*4e-3) # force/cross-sectional-area
Strain = -P_tot/30 # dx/x
print(Stress)
# Use linear least squares to find Young's modulus -> Stress = Y * Strain + offset_error
A = np.vstack([Strain,nAp.ones(len(Strain))]).T
model = np.linalg.lstsq(A, Stress, rcond=None)
Y, offset_error= model[0]
resid = model[1]
# Determine the R_square value between 0 and 1. 1 is a strong correlation
R_sqr = 1 - resid/(Stress.size*Stress.var())
print("Y = %.4f, offset_error = %.4f, R_sqr = %.4f" % (Y, offset_error, R_sqr))

Strain_lin = np.linspace(min(Strain),max(Strain) , 5)
Stress_lin = Y * Strain_lin + offset_error


plt.figure()
ax2 = plt.plot(Strain,Stress,'x',Strain_lin,Stress_lin,'-')

plt.show()
