"""
FILE: stretch_n_measure_smu.py
AUTHOR: R Ellingham
DATE MODIFIED: Apr 2021
PROGRAM DESC: Gather data from a stretched conductive elastomer in real time
DATE CREATED: Oct 2020
using serial communication. Writing the data to a CSV file ready for analysis.

Parameters measured:
-> Resistance, Strain, Force, and their respective time stamps (for each individual measurement)

TODO:
1) Change scheduling system so that force data can be recorded at a much faster rate (10kHz+)
    a) Then we can use this rapid sample data to see any vibration reasonsance
    b) Parallel threads could be handy?
    c) Synchronise all of the 'clocks' of each measurement device
    d) Obtain chunks of data at a time from each measurement device
2) Change AC measurement to a true AC instead of this hacky toggling shit we got going on
3) Push a bunch of the functions in the code into a couple of libraries (linactuator, smuread, loadcell)
4) Make code again from scratch with new scheduling system

"""

import random
import sys, traceback
import csv
import matplotlib.pyplot as plt
import numpy as np
import serial
import serial.tools.list_ports
import time
from datetime import datetime
import re
import pyvisa
import nidaqmx
from nidaqmx.constants import LineGrouping
from nidaqmx.constants import Edge
from nidaqmx.constants import AcquisitionType
from k2600 import K2600 # see k2600.py for usage
import k2600
# from cmeter350610 import Cmeter350610

MAX_LOADCELL_FORCE = 6.5 # Maximum loadcell force in Newtons

### Linear actuator functions:
def write_g(serial_handle, msg):
    """
    DESCR: Sends encoded message to serial_handle device
    IN_PARAMS: serial handle if from serial.Serial(), g-code message to send
    OUTPUT: Serial write status
    NOTES:    Requires serial library
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """
    return serial_handle.write((msg+'\n').encode())

def readchar_g(serial_handle):
        """
        DESCR: Reads decoded character from serial_handle device
        IN_PARAMS: serial handle if from serial.Serial(), g-code message to send
        OUTPUT: Serial read message character
        NOTES:    Requires serial library
                  Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
        """
        return serial_handle.read().decode()

def readbuf_g(serial_handle):
        """
        DESCR: Reads whole decoded message from from serial_handle device
        IN_PARAMS: serial handle if from serial.Serial(), g-code message to send
        OUTPUT: Serial read message list, each item is a newline
        NOTES:    Requires serial library
                  Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
        """
        message_lines = []
        message = serial_handle.readline().decode()
        while message:
            message_lines.append(message)
            message = serial_handle.readline().decode()
        return message_lines

def readline_g(serial_handle):
        """
        DESCR: Reads decoded g-code string until newline character from serial_handle device
        IN_PARAMS: serial handle if from serial.Serial(), g-code message to send
        OUTPUT: Serial read message string
        NOTES:    Requires serial library
                  Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
        """
        return serial_handle.readline().decode()

def init_motion_params(serial_handle, max_accel="50", units="mm", steps_p_mm="250", motion="lin_interp", dist_mode="rel"):
    """
    DESCR: Sets linear actuator limits and params
    IN_PARAMS: max. acceleration, position units, steps per 1mm displacement, linear motion control mode", relative or absolute displacement
    OUTPUT: N/A
    NOTES:  Requires serial library
            Requires time library
            Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """
    write_g(serial_handle, "\r\n\r\n") # Wake up grbl
    time.sleep(2)   # Wait for grbl to initialize
    serial_handle.flushInput()
    #set acceleration
    write_g(serial_handle,"$122="+max_accel)
    #set units
    if units=="mm":
        write_g(serial_handle,"G21")
    elif units=="inch":
        write_g(serial_handle,"G20")
    else:
        print("invalid unit")
    #set steps per mm
    write_g(serial_handle,"$102="+steps_p_mm)
    #set motion mode
    if motion=="lin_interp":
        write_g(serial_handle,"G1") #acceleration and speed limited
    else:
        write_g(serial_handle,"G0") #no limit on acceleration or speed
    #set distance mode
    if dist_mode=="rel":
        write_g(serial_handle,"G91")
    elif dist_mode=="abs":
        write_g(serial_handle,"G90")
    else:
        print("invalid distance mode")
    return 0

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

def linear_travel(serial_handle, lin_speed="60", displacement="1"):
    """
    DESCR: Sends linear motion step to lin actuator
    IN_PARAMS: linear speed in mm/min, displacement in mm.
    OUTPUT: N/A
    NOTES:    Requires serial library and Grbl
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """
    serial_handle.flushInput()
    #set speed and displacement
    write_g(serial_handle,"G21G1"+"F"+str(lin_speed)+"Z"+str(displacement))
    return 0

def read_pos(serial_handle):
    """
    DESCR: Reads current linear position of lin actuator
    IN_PARAMS: N/A
    OUTPUT: Absolute position in mm (string)
    NOTES:  Linear actuator is open loop so position is only estimated with motor steps
            Requires regex (re) library
            Requires serial library
            Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """

    write_g(serial_handle,"?")
    raw_status = readline_g(serial_handle)
    current_position = raw_status.strip('<>\n\r').split(',')[-1] # string manipulation of raw status data
    if re.search('[a-zA-Z]', current_position):
        current_position = 0.0
    else:
        try:
            current_position = float(current_position)
        except ValueError:
            print("could not convert string('%s') to float -> current_pos set to zero" % current_position)
            current_pos = 0.0
    serial_handle.flushInput()
    return current_position

def list_serial_devices():
    """
    DESCR: Creates list of the names of all currrently connected serial comport devices
    IN_PARAMS: N/A
    OUTPUT: List of available comports
    NOTES: Requires serial.tools.list_ports library
    """
    print("Fetching list of devices...")
    avail_devs = serial.tools.list_ports.comports()
    devs_list = []
    for i in range(len(avail_devs)):
        devs_list.append(avail_devs[i].device)
    return devs_list



### 4/2 wire ohmmeter data acquisition functions:
def init_smu_ohmmeter_params(smu_handle,outer_I,outer_Vmax,num_wire=2,nplc=1):
    """
    DESCR: Initialise parameters for the use of the smu as a 2 or 4 wire ohmmeter.
    Using smua as the outer connection (or 2 wire). 2 wire is default.
    IN_PARAMS: smu_handle, Current source current value, current source max
    voltage, number of wire measurement, number of power line cycles per measurment
    OUTPUT: N/A
    NOTES:  Requires pyvisa and k2600 library
    TODO: create a 2wire version
    """
    if num_wire == 2:
        smu_handle.display.screen(smu_handle,smu_handle.display.SMUA) # display both SMUA on screen
        # smu_handle.display.smua.measure.func(smu_handle,smu_handle.display.MEASURE_OHMS) # display ohms measurement

        smu_handle.smua.source.func(smu_handle,smu_handle.smua.OUTPUT_DCAMPS) # set output mode to dc amps (current source)

        smu_handle.smua.source.leveli(smu_handle,outer_I) # set current source value

        smu_handle.smua.source.limitv(smu_handle,outer_Vmax) # set voltage limit (if too low this may restrict the current sourced)

        smu_handle.smua.measure.nplc(smu_handle,nplc) # set number of powerline cycle to integrate over for volts measurement to 1

        smu_handle.smua.source.output(smu_handle,1) # turn on channel a

    elif num_wire == 4:
        smu_handle.display.screen(smu_handle,smu_handle.display.SMUA_SMUB) # display both SMUs on screen
        # smu_handle.display.smua.measure.func(smu_handle,smu_handle.display.MEASURE_DCVOLTS) # display volt measurement on smu
        # smu_handle.display.smub.measure.func(smu_handle,smu_handle.display.MEASURE_DCVOLTS)

        smu_handle.smua.source.func(smu_handle,smu_handle.smua.OUTPUT_DCAMPS) # set output mode to dc amps (current source)
        smu_handle.smub.source.func(smu_handle,smu_handle.smub.OUTPUT_DCAMPS)

        smu_handle.smua.source.leveli(smu_handle,outer_I) # set current source values
        smu_handle.smub.source.leveli(smu_handle,0)

        smu_handle.smua.source.limitv(smu_handle,outer_Vmax) # set voltage limit (if too low this may restrict the current sourced)
        smu_handle.smub.source.limitv(smu_handle,outer_Vmax)

        smu_handle.smua.measure.nplc(smu_handle,nplc) # set number of powerline cycle to integrate over for volts measurement to 1
        smu_handle.smub.measure.nplc(smu_handle,nplc)

        smu_handle.smua.source.output(smu_handle,1) # turn on channel a
        smu_handle.smub.source.output(smu_handle,1) # turn on channel b

    return 0

def init_smu_AC_pulse(smu_handle):
    """
    DESCR: Initialises parameters for an AC pulse train of +/- Isrc
    IN_PARAMS: smu_handle
    OUTPUT: Current resistance reading, average time since start time, time taken for reading
    NOTES:  Requires pyvisa,time and k2600 libraries
    TODO: 1) Add AC measurement mode functionality +ve and -ve pulses
    2) Output timing for v and i measurements to determine if they are approx. 1PLC ea(+message send time)?
    """
    I_src = float(smu_handle.smua.source.leveli(smu_handle)) # get source current value for pulse train
    smu_handle.smua.makebuffer(smu_handle,buf_size,"Vbuffer") # MAKE FUNCTION FOR THIS IN K2600
    smu_handle.smua.trigger.source.limitV(smu_handle,sweepSourceLimitVval)
    smua.trigger.source.listv({-I_src, I_src}) # MAKE FUNCTION FOR THIS IN K2600
    smua.trigger.count = INFINITE # MAKE FUNCTION FOR THIS IN K2600
    smuX.trigger.source.action = smuX.ENABLE# MAKE FUNCTION FOR THIS IN K2600
    #smuX.trigger.arm.set() # MAKE FUNCTION FOR THIS IN K2600 The SMU automatically clears all the event detectors when the smuX.trigger.initiate() function is executed. This function should be called after the sweep is initiated.

def read_smu_res(smu_handle, num_wire=2, mode="NORMAL"):
    """
    DESCR: Gives a resistance reading
    IN_PARAMS: Start time, timeout
    OUTPUT: Current resistance reading, average time since start time, time taken for reading
    NOTES:  Requires pyvisa,time and k2600 libraries
    TODO: 1) Add AC measurement mode functionality +ve and -ve pulses
    2) Output timing for v and i measurements to determine if they are approx. 1PLC ea(+message send time)?
    """
    if mode == "AC": # alternates the direction of the current source each read
        # print("AC mode")
        I_src = -1 * float(smu_handle.smua.source.leveli(smu_handle)) # toggle source current
        # print("I:",I_src)
        smu_handle.smua.source.leveli(smu_handle,I_src) # set current source value
        # print("Isrc toggled")

    outer_res = 0
    inner_res = 0
    t_s = time.time()
    if num_wire == 2:
        ia = float(smu_handle.smua.measure.i(smu_handle))
        va = float(smu_handle.smua.measure.v(smu_handle))
        outer_res = va/ia
        # outer_res = float(smu_handle.smua.measure.r(smu_handle))
    elif num_wire == 4:
        ia = float(smu_handle.smua.measure.i(smu_handle))
        va = float(smu_handle.smua.measure.v(smu_handle))
        vb = float(smu_handle.smub.measure.v(smu_handle))
        outer_res = va/ia
        inner_res = vb/ia
    t_f = time.time()
    t_d = t_f - t_s
    # print("td:",t_d)
    current_res = [outer_res, inner_res]

    return [current_res, t_d]

### Load cell data acquisition functions:
def init_loadcell_params(loadcell_handle):
    """
    DESCR: Initialise parameters for 500g tal221 loadcell
    IN_PARAMS: N/A
    OUTPUT: N/A
    NOTES:  Requires nidaqmx & nidaqmx.constants libraries
    TODO: Add input params
    """
    loadcell_handle.ai_channels.add_ai_force_bridge_two_point_lin_chan('cDAQ1Mod3/ai3',
    name_to_assign_to_channel='500gLoadcell',
    min_val=-5.0, max_val=5.0,
    units=nidaqmx.constants.ForceUnits.NEWTONS,
    voltage_excit_val=6,
    bridge_config=nidaqmx.constants.BridgeConfiguration.FULL_BRIDGE,
    nominal_bridge_resistance=1000.0,
    first_electrical_val=0.0079, second_electrical_val=-0.683,
    electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
    first_physical_val=0.0, second_physical_val=4.905,
    physical_units=nidaqmx.constants.BridgePhysicalUnits.NEWTONS)

    loadcell_handle.timing.cfg_samp_clk_timing(10000, samps_per_chan=10)
    return 0



## Data processing fucntions:
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
    sum = 0
    for item in in_list:
        sum = sum + float(item)
    return sum/len(in_list)

## Zero force calibration function
def auto_zero_cal(loadcell_handle, serial_handle, tolerance):
    """
    DESCR: Rough PI controller to get stress in test rig to zero
    IN_PARAMS: loadcell_handle, serial_handle, zero tolerance [N]
    NOTES: Not well tested.
    TODO: Fix instability - maybe due to 'relative/absolute' travel of linear actuator?
    """
    error = tolerance*2 # start error bigger than tolerance
    error_prev = 0
    buf_size = 40
    f_buf = [0]*buf_size
    speed = 180
    Kp = 4 # control P gain constant
    Ki = 0.5 # control I gain constant
    while abs(error) > tolerance:
        for i in range(buf_size):
            raw_data = loadcell_handle.read(2) # read 1 data point
        force_av = sum(raw_data)/len(raw_data)
        print('force:%.5fN' % force_av)
        error = Kp*force_av + Ki*(error - error_prev)
        print('error:%.5f' % error)
        if error > 2 : error = 2
        if error < -2 : error = -2
        linear_travel(serial_handle, str(speed), str(error)) # (s, "speed" , "dist")
        wait_time = abs(float(error)/float(speed))
        print('t=%.5f' % wait_time)
        prev_error = error
        time.sleep(wait_time) # add a bit of settling time for stress relaxation
    time.sleep(2)
    for i in range(buf_size):
        raw_data = loadcell_handle.read(2) # read 1 data point
    force = sum(raw_data)/len(raw_data)
    print("Final force = %.5fN" % force)
    time.sleep(2)

def main():
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    # print(avail_devs)
    # device_index = int(input("Which linear actuator device from the list? (eg. index 0 or 1 or 2 or ...):"))
    # s = serial.Serial(avail_devs[device_index],115200,timeout=2) # Connect to port. GRBL operates at 115200 baud
    s = serial.Serial("COM8",115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    init_motion_params(s,dist_mode="abs") # Init Grbl and set travel displacement mode to absolute

    # set csv file name for data capture
    filename = input("File name? !Caution will overwrite files without warning!\n (e.g sample number+CB %+electrode type+distance between electrodes+repetitions of test=\n=samp1_CB7-5_Epin_20mm_v2):")

    ## Setup SMU connection
    rm = pyvisa.ResourceManager()
    available_devs = rm.list_resources()
    # print(available_devs)
    # device_index = input("Which device address from the list? (eg. index 0 or 1 or 2 or ...):")
    # ohmmeter = K2600(available_devs[int(device_index)])
    ohmmeter = K2600(available_devs[4])
    meas_wires = 4 # is it a 2 or 4 wire resistance measurement?
    I_src = 10e-6 # constant current source value
    V_max = 20 # max current source value
    meas_mode = "AC"
    init_smu_ohmmeter_params(ohmmeter,I_src,V_max,num_wire=meas_wires)

    ## Setup Hioki 3506-10 C-meter device connection
    # cmeter = Cmeter350610(freq=1E3)

    ## Setup loadcell connection
    loadcell = nidaqmx.Task()
    init_loadcell_params(loadcell)

    # Initialise lists for storing measurement data
    pos_data = []
    time_data_pos = []

    res_data_o = []
    res_data_i = []
    time_data_res = []
    avg_time_data_res = []

    # cap_data = []
    # esr_data = []
    # time_data_cap = []
    # avg_time_data_cap = []

    force_data = []
    time_data_force = []

    # set new zero
    write_g(s,"G90")
    write_g(s,"G10 P0 L20 X0 Y0 Z0")

    ####################################
    ### Set velocity profile params: ###
    ####################################
    ## Total time
    t_tot = 600
    f_range = [0.001,0.2]
    strain_range = [0.5,6]
    num_sines = 8
    max_strain = 14

    t, strain_profile, velocity_profile, ampl_array, freq_array, delta_t = rand_multisine_generator(num_sines, f_range, strain_range, t_tot, max_strain)

    strain_profile = -strain_profile
    velocity_profile = abs(velocity_profile)*60 # only have positive velocities and scale to mm/minute

    # ensure there are no zero velocities in velocity_profile
    for velocity in velocity_profile:
        if velocity < 2:
            velocity = 2
    # ensure first strain point is reached at a constant slow speed for each experiment
    speed2start = 10 # mm/s
    velocity_profile[0] = speed2start

    t_tot = t_tot + (((max_strain/2)/speed2start)*60) # add time to get to initial starting strain point of multisine
    repeats = 1

    # plot strain and strain velocity profiles
    fig, ax = plt.subplots(figsize = (10, 5))
    ax2 = ax.twinx()
    ax.plot(t, -strain_profile, color = 'g')
    ax2.plot(t, velocity_profile, color = 'b')
    ax.set_xlabel('t')
    ax.set_ylabel('strain', color = 'g')
    ax2.set_ylabel('velocity', color = 'b')
    plt.show()

    relax_delay = 60 # amount of time(s) to record the resistive and stress relaxation
    ####################################
    ### End velocity profile params: ###
    ####################################

    ###
    # Begin measurement loop
    ###
    start_time = time.time() # ref time reset for automated test
    try:
        print("Reading data...")
        step_counter = 0

        current_pos = 0 # init for while loop condition

        block_counter = 0
        # send first block of steps
        for step_indx in range(20*block_counter,20*(block_counter+1)): # 20 seems to be the grbl buffer size
            strain = strain_profile[step_indx]
            velocity = velocity_profile[step_indx]
            linear_travel(s, velocity, strain)
            print("step set ", step_indx)
        block_counter = block_counter + 1

        start_loop_time = time.time()
        current_time = start_loop_time + 0.1
        loop_time = 0.0
        pos_err_tol = 0.01
        pos_err = 10 # init pos_err as an appropriately high value

        print("time to start:",(strain_profile[0]/speed2start)*60)

        while (loop_time < t_tot):
            # print("current_time: ",current_time)
            # print("block_counter: ",block_counter)
            # print("delta_t? ",((current_time - (abs(strain_profile[0]/speed2start)*60)) - 20*delta_t*(block_counter)))
            # print("delta_t:", delta_t)
            # # Send one lot of 20 steps every ~19 steps taken. Multiplied by 10 to round the time to the
            # #   nearest 100 millisecond. Current time minus the time taken to get to the start point.
            # #   Threshold for starting next block is when the 19th step is being executed
            if abs((current_time - (abs(strain_profile[0]/speed2start)*60)) - 19*delta_t*(block_counter)) < delta_t:
                blk_s = 20*block_counter
                blk_f = 20*(block_counter+1)

                if (blk_f > len(strain_profile)): # stop out of bounds indexing
                    blk_f = len(strain_profile)
                if (blk_s > len(strain_profile)): # stop out of bounds indexing
                    blk_s = len(strain_profile)

                for step_indx in range(blk_s,blk_f): # 20 seems to be the grbl buffer size
                    strain = strain_profile[step_indx]
                    velocity = velocity_profile[step_indx]
                    linear_travel(s, velocity, strain)
                    print("step set ", step_indx)
                block_counter = block_counter + 1

            loop_time = time.time() - start_loop_time
            print("loop_time-",loop_time)

            # Read resistance
            current_res, t_d = read_smu_res(ohmmeter,num_wire=meas_wires,mode=meas_mode)
            r_stop_time = time.time() - start_time
            t_avg = r_stop_time - t_d/2
            res_data_o.append(current_res[0])
            res_data_i.append(current_res[1])
            avg_time_data_res.append(t_avg)
            time_data_res.append(t_d)

            # Read position
            current_pos = read_pos(s)
            current_time = time.time() - start_time
            pos_data.append(-current_pos)
            time_data_pos.append(current_time)

            print("current_res:",current_res)
            print("current pos:",current_pos)

            # # Read capacitance
            # t_s_cap = time.time()
            # current_cap, current_esr = cmeter.get_measurement()
            # cap_data.append(current_cap)
            # esr_data.append(current_esr)
            # t_f_cap = time.time()
            # t_avg_cap = t_s_cap - start_time + (t_f_cap-t_s_cap)/2
            # avg_time_data_cap.append(t_avg_cap)
            # time_data_cap.append(t_f_cap-t_s_cap)

            # Read force
            t_s_force = time.time()
            raw_f_data = loadcell.read(2) # read 2 data points from buffer
            force_avg = average_list(raw_f_data)
            if (force_avg > MAX_LOADCELL_FORCE):
                raise NameError('Maximum force of {}N for loadcell exceeded'.format(MAX_LOADCELL_FORCE))
            t_f_force = time.time()
            t_d_force = t_f_force - t_s_force
            # print(t_d_force)
            force_data.append(force_avg)
            current_time = time.time() - start_time
            time_data_force.append(current_time)
            # print(force_avg)

            # print("Step complete ", step_counter)

        log_file = open("log.txt","a")
        log_file.write(str(datetime.date(datetime.now()))+'\n')
        log_file.write(str(datetime.time(datetime.now()))+'\n')
        log_file.write('filename='+filename+'\n')
        log_file.write('strain profile='+str(strain_profile)+'\n')
        log_file.write('velocity profile='+str(velocity_profile)+'\n')
        log_file.write('Frequencies:'+str(freq_array)+' Amplitudes:'+str(ampl_array))
        log_file.write('MEASUREMENT:\n num wires='+str(meas_wires)+', Isrc='+str(I_src)+', Vmax='+str(V_max)+', Type='+str(meas_mode)+'\n')
        log_file.write('total time='+str(t_tot)+'\n')
        log_file.write('repeats='+str(repeats)+'\n')
        log_file.write(' \n')
        log_file.close()

        # Write data to CSV file
        write_PosResForce_to_CSV(filename,res_data_o, res_data_i, avg_time_data_res, pos_data, time_data_pos,
            force_data, time_data_force)
        # write_PosResForce_to_CSV(filename,cap_data, esr_data, avg_time_data_cap, pos_data, time_data_pos,
        #     force_data, time_data_force)

        # Plot all data
        plot_q = input("Plot all data?")
        if (plot_q == 'y'):
            fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)

            ax1.plot(avg_time_data_res, res_data_o)
            ax1.set(ylabel='Resistance_outer[Ohm]')

            ax2.plot(avg_time_data_res, res_data_i)
            ax2.set(ylabel='Resistance_inner[Ohm]')

            # ax1.plot(avg_time_data_cap, esr_data)
            # ax1.set(ylabel='ESR [Ohms]')
            #
            # ax2.plot(avg_time_data_cap, cap_data)
            # ax2.set(ylabel='Capacitance[Farads]')

            ax3.plot(time_data_pos, pos_data)
            ax3.set(ylabel='Position[mm]')

            ax4.plot(time_data_force, force_data)
            ax4.set(xlabel='Time[s]', ylabel='Force[N]')

            plt.show()
            input("press enter!")
            fig.savefig(filename)

    except Exception: # if test sequence stops abruptly log it and save data so far.
        print(traceback.format_exc())
        log_file = open("log.txt","a")
        log_file.write(str(datetime.date(datetime.now()))+'\n')
        log_file.write(str(datetime.time(datetime.now()))+'\n')
        log_file.write("ERROR occurred mid sequence! @ time ="+str(time.time() - start_time)+'\n')
        log_file.write('filename='+filename+'\n')
        log_file.write('strain profile='+str(strain_profile)+'\n')
        log_file.write('velocity profile='+str(velocity_profile)+'\n')
        log_file.write('Frequencies:'+str(freq_array)+' Amplitudes:'+str(ampl_array))
        log_file.write('MEASUREMENT:\n num wires='+str(meas_wires)+', Isrc='+str(I_src)+', Vmax='+str(V_max)+', Type='+str(meas_mode)+'\n')
        log_file.write('total time='+str(t_tot)+'\n')
        log_file.write('repeats='+str(repeats)+'\n')
        log_file.write(' \n')
        log_file.close()

    finally: # regardless of program success complete the following
        write_PosResForce_to_CSV(filename,res_data_o, res_data_i, avg_time_data_res, pos_data, time_data_pos,
            force_data, time_data_force)
        # write_PosResForce_to_CSV(filename,cap_data, esr_data, avg_time_data_cap, pos_data, time_data_pos,
        #     force_data, time_data_force)

        init_motion_params(s) # Init Grbl and set travel displacement mode to relative
        # zero the specimen
        auto_zero_cal(loadcell, s, 0.005)
        print("Disconnecting serial and pyvisa connections...")
        # Close smu pyvisa connection
        ohmmeter.disconnect()
        # Close serial port
        s.close()
        print("done")

    return 0

# The real main driver
if __name__ == "__main__":
	# Default input parameters
	# input_filename = "1st_test_num10" main(input_filename, input2)
    # input2 = ??
	# import sys
	# if len(sys.argv)>1: input_filename   = (sys.argv[1])
	# if len(sys.argv)>2: input2   =int(sys.argv[2])

	main()
