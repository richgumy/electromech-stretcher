"""
Just absolute garbage.py
"""
import sys
import csv
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import random

def sinusoid_maker(amplitude, freq, offset, t_tot, delta_t):
    """
    DESCR: Setup grbl parameters to apply a quantised sinusoidal waveform
    IN_PARAMS: amplitude of sinusoid, frequency in Hz, offset_strain in (default = amplitude/2)
    OUTPUT: time, x and dx_dt arrays
    NOTES:    Requires serial library and Grbl
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
              Cannot produce negative values for x
    """
    n = int(t_tot / delta_t) # samples per wavelength
    t = np.linspace(0,t_tot,n+1) # wavelength time array
    x_prev = 0 # previous x used for differentiation
    x = np.zeros(n+1) # output step array
    dx_dt = np.zeros(n+1) # output velocity array
    for i in range(n+1):
        x[i] = amplitude*np.sin(2*np.pi*freq*t[i]) + offset
        dx_dt[i] = (x[i] - x_prev)/(delta_t)
        if i == 1:
            dx_dt[i-1] = dx_dt[i]
        x_prev = x[i]
    return t, x, dx_dt


def rand_multisine_generator(num_sines, freq_range, amp_range, t_tot, max_amp):
    """
    DESCR: Generate a quantised multisinusoidal strain waveform
    IN_PARAMS: number of sinusoids, range of sinusoid amplitudes [mm] in form [min,max],
    frequency range of sinusoids [Hz] in form [min,max], total duration of multisine
    OUTPUT: time, strain and velocity profiles. Amplitude and frequency arrays
    NOTES:    Requires libraries: serial and Grbl and random
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
              Cannot produce compresive strain only tensile
              Can use the amplitude and frequency arrays to confirm that FFT on the data makes sense
    TODO:
     1. Make all generated start from zero
        -> move extremely slowly at a constant speed to the start point of waveform strain
            OR time shift all sinusoids apart from one
            OR shift to start point and wait for resistance to relax until starting experiment
    2. Currently has an unlimited size for the output arrays, which can be too large for arduino
        GRBL to handle. GRBL can handle a max. of approx 50 g-code commands. Yet to confirm...
    """
    ampl_array = np.zeros(num_sines)
    freq_array = np.zeros(num_sines)

    # Minimize the amount of data points required (fs >= 10*max(freq_range))
    Ts = 1 / (10 * max(freq_range))

    # makes a 't' array of correct size
    t, test_sine, test_v_sine = sinusoid_maker(1, 1, 1, t_tot, Ts)

    sine_array = np.zeros([num_sines, len(t)])
    sine_velocity_array = np.zeros([num_sines, len(t)])

    # generate a set number of 'random' sinusoids
    for i in range(num_sines):
        ampl_array[i] = round(random.uniform(amp_range[0],amp_range[1]))
        freq_array[i] = round(random.uniform(freq_range[0],freq_range[1]),3)
        t, sine_array[i], sine_velocity_array[i] = sinusoid_maker(ampl_array[i], freq_array[i], ampl_array[i], t_tot, Ts)

    # sum all sines and ensure min value is set to zero.
    sine_array_tot = sum(sine_array) - min(sum(sine_array))
    sine_velocity_tot = sum(sine_velocity_array)

    # scale down sine array such that the max. value is equal to max_amp
    scale_factor = max(sine_array_tot) / max_amp
    sine_array_tot = sine_array_tot / scale_factor
    sine_velocity_tot = sine_velocity_tot / scale_factor

    print("Characteristics of sinusoids:")
    print("Amplitudes:", ampl_array)
    print("Frequencies:", freq_array)
    print("Scale_factor:", scale_factor)
    print("Max. speed:",max(sine_velocity_tot))
    print("Max. strain:",max(sine_array_tot))
    print("Size of array:",len(sine_array_tot))
    input("Press Ctrl+C if these values are not ok, else press Enter")

    return t, sine_array_tot, sine_velocity_tot, ampl_array, freq_array


# # Combine all strain and velocity profiles to create a multi-sine waveform
# strain_profile_tot = strain_profile0 + strain_profile1 + strain_profile2
# velocity_profile_tot = velocity_profile0 + velocity_profile1 + velocity_profile2

t, strain_profile_tot, velocity_profile_tot, ampl_array, freq_array = rand_multisine_generator(6, [0.2,5], [0.5,6], 5, 12)

strain_profile_tot = -strain_profile_tot
# plotting figures by creating aexs object
# using subplots() function
fig, ax = plt.subplots(figsize = (10, 5))
# plt.title('Example of Two Y labels')

# using the twinx() for creating another
# axes object for secondry y-Axis
ax2 = ax.twinx()

ax.plot(t, strain_profile_tot, color = 'g')
ax2.plot(t, velocity_profile_tot, color = 'b')

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
