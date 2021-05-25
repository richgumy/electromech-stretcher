"""
Just absolute garbage.py
"""
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

def sinusoid_maker(amplitude_strain, freq, offset_strain):
    """
    DESCR: Setup grbl parameters to apply a quantised sinusoidal strain waveform
    IN_PARAMS: amplitude_strain of sinusoid in mm, frequency in Hz, offset_strain in mm (default = amplitude_strain/2)
    OUTPUT: velocity in mm/s
    NOTES:    Requires serial library and Grbl
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
              Cannot produce compresive strain only tensile
    """
    delta_t = 0.01 # hard code the quantisation time
    T = 1 / freq # period
    n = int(T / delta_t) # samples per wavelength
    t = np.linspace(0,T,n+1) # wavelength time array
    x_prev = 0 # previous x used for differentiation
    x = np.zeros(n+1) # output step array
    dx_dt = np.zeros(n+1) # output velocity array
    for i in range(n+1):
        x[i] = amplitude_strain*np.sin(2*np.pi*freq*t[i]) + offset_strain
        dx_dt[i] = (x[i] - x_prev)/(delta_t)
        if i == 1:
            dx_dt[i-1] = dx_dt[i]
        x_prev = x[i]
    return t, x, dx_dt

t, strain_profile, velocity_profile = sinusoid_maker(2,1,2)

# plotting figures by creating aexs object
# using subplots() function
fig, ax = plt.subplots(figsize = (10, 5))
# plt.title('Example of Two Y labels')

# using the twinx() for creating another
# axes object for secondry y-Axis
ax2 = ax.twinx()
ax.plot(t, strain_profile, color = 'g')
ax2.plot(t, velocity_profile, color = 'b')

# giving labels to the axises
ax.set_xlabel('t')
ax.set_ylabel('strain', color = 'g')

# secondary y-axis label
ax2.set_ylabel('velocity', color = 'b')

# defining display layout
plt.tight_layout()

# show plot
plt.show()

# def MAF(x,N):
#     """
#     DESCR: Moving average filter for a set of data averaging +/- N/2 each point
#     IN_PARAMS: x data, sample period to average over
#     NOTES: Requires numpy library. Use an odd numbered N size window.
#     TODO:
#     """
#     x_padded = np.pad(x, (N//2, N-1-N//2), mode='edge')
#     x_smooth = np.convolve(x_padded, np.ones((N,))/N, mode='valid')
#
#     return x_smooth
#
# data = [1,2,3,4,5,6,7,8]
#
# new_data = MAF(data,3)
#
# print(data)
# print(new_data)

# def f(x, a0, a1, tau_1): # generalised Kelvin SLS relaxation model for n = 3
#     return a0 + a1 * np.exp(-x/tau_1)
#
# t = np.linspace(0,30000,1000)
# y = f(t,170,1320,17400)
#
# plt.plot(t,y)
#
# plt.show()

#
# step_profile = [1,2,3]
# filename = "test this juan"
#
# log_file = open("log.txt","a")
# log_file.write(str(datetime.date(datetime.now()))+'\n')
# log_file.write(str(datetime.time(datetime.now()))+'\n')
# log_file.write('filename='+filename+'\n')
# log_file.write('step profile='+str(step_profile)+'\n')
# log_file.write('velocity profile='+str(velocity_profile)+'\n')
# log_file.write('diff_min(convergence checker)='+str(diff_min)+'\n')
# log_file.write('iter_max(convergence timeout)='+str(iter_max)+'\n')
# log_file.close()


# def main(input_filename):
#     if (input_filename[-4:] != '.csv'): input_filename = input_filename + '.csv' # Assume file is .csv if not explicitly stated in main parameter input
#     print(input_filename)
#     Rout = np.array([])
#     Rin = np.array([])
#     tR = np.array([])
#     P = np.array([])
#     tP = np.array([])
#     F = np.array([])
#     tF = np.array([])
#
#     with open(input_filename, 'r', newline='') as csvfile:
#         data = csv.reader(csvfile, delimiter=',')
#         line_count = 0
#         for row in data:
#             if line_count == 0:
#                 # print(f'Column names are {", ".join(row)}')
#                 line_count = 1
#             else:
#                 R = row[0].strip('][').split(', ')
#                 Rout = np.append(Rout,float(R[0]))
#                 Rin = np.append(Rin,float(R[1]))
#                 tR = np.append(tR,float(row[1]))
#                 P = np.append(P,float(row[2]))
#                 tP = np.append(tP,float(row[3]))
#                 F = np.append(F,float(row[4]))
#                 tF = np.append(tF,float(row[5]))
#
#     plt.plot(tR, Rout)
#     # plt.plot(Rin, tR)
#     plt.show()
#
# # The real main driver
# if __name__ == "__main__":
# 	# Default input parameters
# 	input_filename = ""
# 	if len(sys.argv)>1: input_filename   = (sys.argv[1])
# 	# if len(sys.argv)>2: input2   =int(sys.argv[2])
#
# 	main(input_filename)
