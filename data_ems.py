"""
FILE: data_ems.py
AUTHOR: R Ellingham
DATE MODIFIED: June 2021
PROGRAM DESC: Library for a few data processing functions for the electromech stretcher
project
DATE CREATED: June 2021

TODO:
"""

import csv
import numpy as np
import random
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from scipy import signal

## Preprocessing functions

def write_PosResForce_to_CSV(filename,resistance_outer, resistance_inner, time_R, displacement, time_d, force, time_f):
    """
    DESCR: Loads 6 columns of data into a filename.csv file
    IN_PARAMS: resistance(outer electrodes), resistance (inner electrodes), displacement, force, and their timestamps
    OUTPUT: N/A
    NOTES:  Requires all inputs to be present to operate
    TODO: Make more generalised function for different data logging
    """
    with open(filename+'.csv', 'w', newline='') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        for i in range(len(time_R)):
            data.writerow([resistance_outer[i],resistance_inner[i], time_R[i], displacement[i], time_d[i],
                force[i], time_f[i]])
    return 0

def stringify_list(in_list, dimension=1):
    """
    DESCR: Turn all elements from a 1,2,or 3 dimension list into a list of strings
    IN_PARAMS: list, dimension of list
    OUTPUT: new list of strings
    NOTES:  Requires all evenly dimensioned list (e.g. couldn't have >>stringify_list([[1,2],1], 2))
    TODO: Automatically detect dimension of input list
    """
    new_list1 = []
    for i1 in range(len(in_list)):
        if dimension > 1:
            new_list2 = []
            for i2 in range(len(in_list[i1])):
                if dimension > 2:
                    new_list3 = []
                    for i3 in range(len(in_list[i1][i2])):
                        new_list3.append(str(in_list[i1][i2][i3]))
                    new_list2.append(new_list3)
                else:
                    new_list2.append(str(in_list[i1][i2]))
            new_list1.append(new_list2)
        else:
            new_list1.append(str(in_list[i1]))
    return new_list1

def average_list(in_list):
    """
    DESCR: Finds average of a list of numbers
    IN_PARAMS: list of numbers (int, float, ...)
    OUTPUT: Average value
    """
    sum = 0
    for item in in_list:
        sum = sum + float(item)
    return sum/len(in_list)

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
    DESCR: Geenerate a quantised multisinusoidal strain waveform
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
        GRBL to handle. GRBL can handle a max. of approx 50?? g-code commands. Yet to confirm...
    3. Output velocity is offset!! Is should be centered around zero!!
    """
    ampl_array = np.zeros(num_sines)
    freq_array = np.zeros(num_sines)

    # Minimize the amount of data points required (fs >= 10*max(freq_range))
    Ts = 1 / (20 * max(freq_range))

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

    return t, sine_array_tot, sine_velocity_tot, ampl_array, freq_array, Ts


## Postprocessing functions

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
    DESCR: Fixes any randomly seen zero strain values
    IN_PARAMS: Raw position array
    NOTES:
    TODO:
    """
    new_strain_arr = strain_arr
    for i in range(1,len(strain_arr)):
        if abs(strain_arr[i] - strain_arr[i-1]) > 1:
            new_strain_arr[i] = strain_arr[i-1]
    return new_strain_arr

def interpolate_RFS_data(Ro,Ri,tR,P,tP,F,tF):
    """
    DESCR: Interpolates measurement data that has an aperiodic sampling period
    IN_PARAMS: Resistance outer, resistance inner, resistance time, position, 
    position time, force, force time.
    OUTPUT: Interpolated numpy arrays for each measurement to get a constant 
    sampling period
    NOTES: Interpolation will introduce error (haven't noticed any significant 
                                                error as of yet)
    TODO:
    """
    # Interpolate all data
    Fintp_P = np.interp(tP,tF,F)
    Fintp_R = np.interp(tR,tF,F)
    
    Pintp_F = np.interp(tF,tP,P)
    Pintp_R = np.interp(tR,tP,P)
    
    # timing data for resistance for inital testing is averaged (i.e. R and Ri measurements assumed to be taken at same time)
    Rointp_P = np.interp(tP,tR,Ro)
    Rointp_F = np.interp(tF,tR,Ro)
    
    Riintp_P = np.interp(tP,tR,Ri)
    Riintp_F = np.interp(tF,tR,Ri)
    
    # Generate new interpolated data for the '_tot' arrays
    # interpolation order is [force,pos,res]
    F_tot = np.append(F,[Fintp_P,Fintp_R])
    tF_tot = np.append(tF,[tP,tR])
    
    P_tot = np.append(Pintp_F,[P,Pintp_R])
    tP_tot = np.append(tF,[tP,tR])
    
    Ro_tot = np.append(Rointp_F,[Rointp_P,Ro])
    tR_tot = np.append(tF,[tP,tR])
    
    Ri_tot = np.append(Riintp_F,[Riintp_P,Ri])
    
    ## Sort all _tot data into time order
    # Make new arrays for sorted data
    F_ = np.zeros(len(tF_tot))
    P_ = np.zeros(len(tF_tot))
    Ro_ = np.zeros(len(tF_tot))
    Ri_ = np.zeros(len(tF_tot))
    t_ = np.zeros(len(tF_tot))
    
    ind = np.argsort(tF_tot)
    for i in range(len(tF_tot)):
        F_[i] = F_tot[ind[i]]
        P_[i] = P_tot[ind[i]]
        Ro_[i] = Ro_tot[ind[i]]
        Ri_[i] = Ri_tot[ind[i]]
        t_[i] = tF_tot[ind[i]]
    
    # Interpolate again to get time in a linear form
    t_lin = np.linspace(min(t_),max(t_),int(len(t_)/3)) # divide by 3 to get original size of data
    F_t = np.interp(t_lin,t_,F_)
    P_t = abs(np.interp(t_lin,t_,P_))
    Ro_t = abs(np.interp(t_lin,t_,Ro_))
    Ri_t = abs(np.interp(t_lin,t_,Ri_))
    
    return t_lin, F_t, P_t, Ro_t, Ri_t

def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes: 
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows

def notch_filtering(wav, fs, w0, Q):
    """ Apply a notch (band-stop) filter to the signal.
    
    Args:
        wav: Waveform.
        fs: Sampling frequency of the waveform.
        w0: The frequency to filter. See scipy.signal.iirnotch.
        Q: See scipy.signal.iirnotch.
        
    Returns:
        wav: Notch filtered waveform.
    """
    b, a = signal.iirnotch(2 * w0/fs, Q)
    wav = signal.lfilter(b, a, wav)
    return wav


def write_processed_data(input_filename, Ro_t, Ri_t, resistivity_t, Strain_t, Strain_true_t, Stress_pois_t, Stress_eng_t, Stress_true_t, t_lin):
    """
    DESCR: Write interpolated processed data to CSV
    IN_PARAMS: Ro_t, Ri_t, resistivity_t, Strain_t, Strain_true_t, Stress_pois_t, Stress_eng_t, Stress_true_t, t_lin
    OUTPUT: None
    NOTES:
    TODO:
    """
    with open('INTERP_'+input_filename, 'w', newline='') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        for i in range(len(t_lin)):
            data.writerow([Ro_t[i],Ri_t[i], resistivity_t[i], Strain_t[i], Strain_true_t[i], Stress_pois_t[i], Stress_eng_t[i], Stress_true_t[i], t_lin[i]])
    return 0