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


# import random
import sys, traceback
# import csv
import matplotlib.pyplot as plt
# import numpy as np
import serial
# import serial.tools.list_ports
import time
from datetime import datetime
import re
import pyvisa
import nidaqmx

# Homebrew libraries:
import data_ems
import linact_grbl
from k2600 import K2600 # see k2600.py for usage
import k2600
# from cmeter350610 import Cmeter350610 # see cmeter350610.py for usage

MAX_LOADCELL_FORCE = 6.5 # Maximum loadcell force in Newtons

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
    TODO: 
        1) Add AC measurement mode functionality +ve and -ve pulses
        2) Output timing for v and i measurements to determine if they are approx. 1PLC ea(+message send time)?
    
    NOT COMPLETE
    
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

##############################################MAIN##############################################

def main():
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    # print(avail_devs)
    # device_index = int(input("Which linear actuator device from the list? (eg. index 0 or 1 or 2 or ...):"))
    # s = serial.Serial(avail_devs[device_index],115200,timeout=2) # Connect to port. GRBL operates at 115200 baud
    s = serial.Serial("COM8",115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    linact_grbl.init_motion_params(s,dist_mode="abs") # Init Grbl and set travel displacement mode to absolute

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
    linact_grbl.write_g(s,"G90")
    linact_grbl.write_g(s,"G10 P0 L20 X0 Y0 Z0")

    ####################################
    ### Set velocity profile params: ###
    ####################################
    ## Total time
    t_tot = 600
    f_range = [0.001,0.2]
    strain_range = [0.5,6]
    num_sines = 8
    max_strain = 14

    t, strain_profile, velocity_profile, ampl_array, freq_array, delta_t = data_ems.rand_multisine_generator(num_sines, f_range, strain_range, t_tot, max_strain)

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
            linact_grbl.linear_travel(s, velocity, strain)
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
                    linact_grbl.linear_travel(s, velocity, strain)
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
            current_pos = linact_grbl.read_pos(s)
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
            force_avg = data_ems.average_list(raw_f_data)
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
        data_ems.write_PosResForce_to_CSV(filename,res_data_o, res_data_i, avg_time_data_res, pos_data, time_data_pos,
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
        data_ems.write_PosResForce_to_CSV(filename,res_data_o, res_data_i, avg_time_data_res, pos_data, time_data_pos,
            force_data, time_data_force)
        # write_PosResForce_to_CSV(filename,cap_data, esr_data, avg_time_data_cap, pos_data, time_data_pos,
        #     force_data, time_data_force)

        linact_grbl.init_motion_params(s) # Init Grbl and set travel displacement mode to relative
        # zero the specimen
        linact_grbl.auto_zero_cal(loadcell, s, 0.005)
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
