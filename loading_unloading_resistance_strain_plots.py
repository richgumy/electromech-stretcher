"""
FILE: data_relaxation_fitting_spyderv2.py
AUTHOR: R Ellingham
DATE: Dec 2020
PROGRAM DESC: Fit a relaxation model to stress relaxation data.

TODO:
1)
2)
"""

import csv
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize
from scipy import signal

import matplotlib.lines as mlines
import matplotlib.patches as mpatches

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




font = {'family' : 'normal',
        'size'   : 16}

matplotlib.rc('font', **font)

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
    TODO: Make function work for saw tooth data (i.e. if sign != diff(strain)/abs(diff(strain)): set index; sign = -sign)
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
    IN_PARAMS: Raw posistion array
    NOTES:
    TODO:
    """
    new_strain_arr = strain_arr
    for i in range(1,len(strain_arr)):
        if abs(strain_arr[i] - strain_arr[i-1]) > 1:
            new_strain_arr[i] = strain_arr[i-1]
    return new_strain_arr

def notch_filtering(wav, fs, w0, Q):
    """ Apply a notch (band-stop) filter to the signal.
    
    Args:
        wav: Waveform.
        fs: Sampling frequency of the waveform.
        w0: The frequency to filter. See scipy.signal.iirnotch.
        Q: See scipy.signal.iirnotch.
        
    Returns:
        wav: Filtered waveform.
    """
    b, a = signal.iirnotch(2 * w0/fs, Q)
    wav = signal.lfilter(b, a, wav)
    return wav

### MAIN CODE STARTS HERE ###

# Test specimen dimensions:
spec_length = 40e-3
spec_width = 10e-3
spec_thickness = 4e-3

# input_filename = "1_CB0_v1_0.3Strain.csv"
# input_filename = "2_7-5_E4pin_20mm_v13_0.1_0.2_0.3_strain.csv"
# input_filename = "1_10_E4pin_20mm_v9_0.3Strain.csv"
# input_filename = "2_7-5_4Epin_rand_sawtooth.csv"
# input_filename = "2_7-5_E4pin_20mm_v19_0.3Strain.csv"
# input_filename = "2_7-5_4Epin_20mm_quasistatic_v1.csv"
input_filename = "2_7-5_4Epin_20mm_10prestrain3.csv"

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
            P[i] = float(row[3])
            tP[i] = float(row[4])
            F[i] = float(row[5])
            tF[i] = float(row[6])
            print("%.2f" % (i/len(data)))
            i = i + 1
            # print(row[0],row[1],row[2],row[3],row[4],row[5])
            
P = pos_outlier_corrector(P)

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
Ri_t = np.interp(t_lin,t_,Ri_)

# Use a MAF on the resistance data if it is an AC measurement
Ri_t = abs(MAF(Ri_t,8))

# Calc strain from displacement
Strain_t = abs(P_t/(spec_length*1000)) # dx/x

# Calc stress from force and changing strain
poisson_ratio = 0.29 # Poisson's ratio. Found experimentally using #2_7.5%dragonskin10NV specimen
d_A = ((spec_width*spec_thickness)*(-Strain_t*poisson_ratio+1)*(-Strain_t*poisson_ratio+1))
Stress_t = F_t/d_A # force/cross-sectional-area
resistivity_t = Ri_t/d_A

# Determine strain gauge factor dRes/dStrain
Gauge_factor = ((Ri_t - min(Ri_t))/Ri_t)
Gauge_factor_t = np.zeros(len(Ri_t))
for i in range(len(Gauge_factor)):
    if Strain_t[i] == 0:
        Gauge_factor_t[i] = 0
    else:
        Gauge_factor_t[i]  = Gauge_factor[i] / Strain_t[i]
        
# Filter osicllation seen in strain ramp part of data with frequency of approx. 1.85Hz
fs = 1/(t_lin[1]-t_lin[0])
f0 = 1.85 # notch freq
Q = f0 / 0.5 # notch freq / bandwidth
Strain_t_fil = notch_filtering(Strain_t,fs,f0,Q)

# Plot measurements over time
fig1, axs1 = plt.subplots(4, 1, constrained_layout=True)

ax = axs1[0]
ax.plot(t_lin, R_t,'r-')
ax.set_title('')
ax.set_ylabel('Resistance outer [Ohm]')
ax.grid(True)

ax = axs1[1]
ax.plot(t_lin, Ri_t,'b-')#),t_lin, Ri_t_fil,'b-')
ax.set_title('')
ax.set_ylabel('Resistance inner [Ohm]')
ax.set_xlabel('Time [s]')
ax.grid(True)

ax = axs1[2]
ax.plot(t_lin, Stress_t,'r-')
ax.set_title('')
ax.set_ylabel('Stress [Pa]')
ax.grid(True)

ax = axs1[3]
ax.plot(t_lin, Strain_t,'r-', t_lin, Strain_t_fil, 'b')
ax.set_title('')
ax.set_ylabel('Strain')
ax.grid(True)

## Plot relaxation and fit curve
# Split data into piece-wise data
# strain_splits = split_ramp_data(Strain_t)
# for i in range(len(strain_splits)):
#     print(t_lin[strain_splits[i]])

strain_splits = []
# for i in range(9):
#     strain_splits.append(i*160+4)
#     strain_splits.append(i*160+4)

# strain_splits = [0,0,161,161,322,322,403,403,483,483,536,536,590,590]
time_splits = [0, 24.2, 42.31, 51, 60.66, 66.9, 73.16, 77.9, 82.66, 100.77, 118.88, 128.2, 137.2, 143.45, 149.7, 154.95, 159.7, 177.81, 195.9, 205.2, 214.4, 220.61, 226.86, 231.11, 236.2] # add more
for i in range(len(time_splits)):
    strain_splits.append(time_splits[i]/(t_lin[1]-t_lin[0]))
    strain_splits.append(time_splits[i]/(t_lin[1]-t_lin[0]))

for i in range(len(strain_splits)):
    strain_splits[i] = strain_splits[i].astype(int)

# take chunks of stress and resistance relaxing values from index i1 to i2
i1 = []
i2 = []
for i in range(int(len(strain_splits)/2)):
    i1.append(int(2*i + 1))
    i2.append(int(2*i + 2))
    
print(i1,i2)

## Plotting and fitting the loading and unloading of the different specimens

def lin_func(x,m,c):
    return m*x + c

consts_lin_fit = []
p_init = [1,1]

constsu_lin_fit = []
pu_init = [1,1]

start_offset = 0 # index offsets
end_offset = 0

plot_counter = 0

try:
    fig3, ax3 = plt.subplots(figsize=(16,8))
    for i in range(0,int(len(strain_splits)/4)):
        Strain_load = Strain_t[int(strain_splits[4*i+1])-start_offset:int(strain_splits[4*i+2])+end_offset]
        Strain_unload = Strain_t[int(strain_splits[4*i+3])-start_offset:int(strain_splits[4*i+4])+end_offset]

        Res_load = Ri_t[int(strain_splits[4*i+1])-start_offset:int(strain_splits[4*i+2])+end_offset]
        # Res_load_min = np.min(Res_load)
        # Res_load = Res_load - Res_load_min
        Res_unload = Ri_t[int(strain_splits[4*i+3])-start_offset:int(strain_splits[4*i+4])+end_offset]
        # Res_unload_min = np.min(Res_load)
        # Res_unload = Res_load - Res_load_min
        
        # Res-strain loading fitting
        poptS, pcovS = optimize.curve_fit(lin_func, Strain_load, Res_load, p0=p_init, maxfev=50000)
        p_init = poptS
        Strain_load_lin = np.linspace(min(Strain_load),max(Strain_load) , 100)
        Res_load_lin = lin_func(Strain_load_lin,*poptS)
        consts_lin_fit.append(poptS) # Store fitted parameters of each relaxation
        
        # Res-strain unloading fitting
        poptSu, pcovSu = optimize.curve_fit(lin_func, Strain_unload, Res_unload, p0=pu_init, maxfev=50000)
        pu_init = poptSu
        Strain_unload_lin = np.linspace(min(Strain_unload),max(Strain_unload) , 100)
        Res_unload_lin = lin_func(Strain_load_lin,*poptSu)
        constsu_lin_fit.append(poptSu) # Store fitted parameters of each relaxation
        
        ## Plot res-strain loading/unloading of specimen
        plot_counter = plot_counter + 1
        if plot_counter == 5:
            plot_counter = 1
        plot_colour = {1:'#ff0000',2:'#00ff00',3:'#0000ff',4:'0.3',5:'#ff4444',6:'#44ff44',7:'#4444ff',8:'0.7'}           
            
        # fig3, ax3 = plt.subplots(figsize=(16,8))
        # Loading
        line_load, = ax3.plot(100*Strain_load,Res_load,color=plot_colour[plot_counter],marker='',ls='-')
        # ax3.plot(Strain_load_lin,Res_load_lin,color='y',ls='-')
        # Unloading
        line_unload, = ax3.plot(100*Strain_unload,Res_unload,color=plot_colour[4+plot_counter],marker='',ls='-')
        # ax3.plot(Strain_unload_lin,Res_unload_lin,color='g',ls='-')
        # ax3.legend(["Loading", "Unloading"])
        
        # line_load, = ax3.plot(100*Strain_t,Ri_t,color='r',marker='',ls='-')
        ax3.set_ylabel('Resistance [Ohm]')
        ax3.set_xlabel('Strain [%]')
        add_arrow_to_line2D(ax3, line_load, arrow_locs=np.linspace(0., 1., 200),
                        arrowstyle='->')
        add_arrow_to_line2D(ax3, line_unload, arrow_locs=np.linspace(0., 1., 200),
                        arrowstyle='->')

finally:

    plt.show()

## Plot colour gradient that matches dStrain_dt
# for i in range(len(t_lin)-1):
#     hex_mid = str(hex(int(dS_dt_colour[i])))[2:]
#     if len(hex_mid) == 1:
#         hex_mid = "0" + hex_mid
#     hex_clr = "#ff" + hex_mid + "00"
#     print(clr)
#     line_load, = ax3.plot(Strain_t[i],Ri_t[i],color=hex_clr,marker='.')