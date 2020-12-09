"""
FILE: main.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: Gather data from a stretched conductive elastomer in real time
using serial communication. Writing the data to a CSV file ready for analysis

Parameters measured:
-> Resistance, Strain, Force, and Time (for each individual measurement)

TODO:
1) Complete data processing functions
2) Minimise all instructions while data being recorded (use post-processing as
much as possible)
3) Make object oriented GUI so other people can use it?
4) rename prog
"""

import csv
import matplotlib.pyplot as plt
import serial
import serial.tools.list_ports
import time
import re
import pyvisa
import nidaqmx
from nidaqmx.constants import LineGrouping
from nidaqmx.constants import Edge
from nidaqmx.constants import AcquisitionType

MAX_LOADCELL_FORCE = 100 # Maximum loadcell force in Newtons

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
    current_position = raw_status.strip('<>\n\r').split(',')[-1] #string manipulation
    if re.search('[a-zA-Z]', current_position):
        current_position = 0
    else:
        # print("Position data handler:",current_position)
        current_position = float(current_position)
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



### Ohmmeter data acquisition functions:
def init_ohmmeter_params(ohmmeter_handle):
    """
    DESCR: Initialise parameters for ohmmeter
    IN_PARAMS: N/A
    OUTPUT: N/A
    NOTES:  Requires pyvisa library
    TODO: Add input params
    """
    ohmmeter_handle.baud_rate = 115200
    ohmmeter_handle.timeout = 5000

    ohmmeter_handle.write(":CONF:RES AUTO,1000,(@101)") #set range and resolution of resistance measurements
    ohmmeter_handle.write(":ROUT:CHAN:DEL 0,(@101)") #set delay between scans
    ohmmeter_handle.write(":SENS:RES:NPLC 1,(@101)") #set PLC(i.e. 50Hz power line cycle for NZ) to 0.02 for 400us integration time... or higher for better noise suppression
    ohmmeter_handle.write("ROUT:SCAN (@101)") # Set scan channel
    ohmmeter_handle.write(":TRIG:COUNT 1") # Set number of measurements taken each scan

    return 0

def read_ohmmeter(ohmmeter_handle,start_time, timeout=10000):
    """
    DESCR: Read ohmmeter resistance, will timeout if nothing recieved within timeout period (default 10s)
    IN_PARAMS: Start time, timeout
    OUTPUT: Current resistance reading, average time since start time, time taken for reading
    NOTES:  Requgires pyvisa library
    TODO: Add input params
    """
    t_s = time.time()
    ohmmeter_handle.write(":INIT") # Begin scan
    time.sleep(0.002) # Wait for scan to finish
    current_res = ohmmeter_handle.query("R?")
    time_out_R = 0;
    while len(current_res) <= 4: # If scan not finished keep requesting until result appears
        current_res = ohmmeter_handle.query("R?")
        time_out_R = time_out_R + 1
        if time_out_R > 10000:
            print("Ohmmeter timeout error.")
            break
    current_res = float(current_res.strip()[5:-4])*10**float(current_res.strip()[-1:]) # Waste of processing time
    t_f = time.time()
    t_avg = (t_f + t_s)/2-start_time # Record time of measurement halfway between when measurement was requested and when it was received, assuming continuous integration
    t_d = t_f - t_s # How long did it take to receive the message from time of request
    return [current_res, t_avg, t_d]

### Load cell data acquisition functions:
def init_loadcell_params(loadcell_handle):
    """
    DESCR: Initialise parameters for loadcell
    IN_PARAMS: N/A
    OUTPUT: N/A
    NOTES:  Requires nidaqmx & nidaqmx.constants libraries
    TODO: Add input params
    """
    #adding ABRITRARY linear table of values for load cell reading, got lazy and will scale later
    loadcell_handle.ai_channels.add_ai_force_bridge_table_chan("cDAQ2Mod3/ai3",
        name_to_assign_to_channel="10kgLoadcell",
        min_val=-5.0e-3, max_val=5.0e-3,
        bridge_config=nidaqmx.constants.BridgeConfiguration.FULL_BRIDGE,
        voltage_excit_source= nidaqmx.constants.ExcitationSource.INTERNAL,
        voltage_excit_val=2.5, nominal_bridge_resistance=1000.0,
        electrical_vals= [0, 1, 2, 3, 4],
        electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
        physical_vals=[-1, 0, 1, 2, 3],
        physical_units=nidaqmx.constants.BridgePhysicalUnits.NEWTONS,
        custom_scale_name=None)
    loadcell_handle.timing.cfg_samp_clk_timing(100, source="",
        active_edge=Edge.RISING, sample_mode=AcquisitionType.FINITE,
        samps_per_chan=5)
    return 0



## Data processing fucntions:
def write_PosResForce_to_CSV(filename,resistance, time_R, displacement, time_d, force, time_f):
    """
    DESCR: Loads 6 columns of data into a filename.csv file
    IN_PARAMS: resistance, displacement, force, and their timestamps
    OUTPUT: N/A
    NOTES:  Requires all inputs to be present to operate
    TODO: Make more generalised function for different data logging
    """
    with open(filename+'.csv', 'w', newline='') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        for i in range(len(time_R)):
            data.writerow([resistance[i], time_R[i], displacement[i], time_d[i],
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

def auto_zero_cal(loadcell_handle, serial_handle, tolerance):
    """
    DESCR: Basic proportional controller to calibrate zero force position for
    stretching device.
    IN_PARAMS: loadcell_handle, serial_handle, zero tolerance [N]
    NOTES: Not well tested.
    TODO:
    """
    error = tolerance*2 # start error bigger than tolerance
    buf_size = 100
    f_buf = [0]*buf_size
    speed = 100
    Kp = 1 # control gain constant
    while abs(error) > tolerance:
        for i in range(buf_size):
            raw_data = loadcell_handle.read(1) # read 1 data point
            f_buf[i] = 452.29*float(raw_data[0]) + 98.155 # magic nums to get force into newtons
            sleep(0.001) # Allow material to 'settle' and gain a better average
        force_av = sum(f_buf)/len(f_buf)
        print(force_av)
        error = Kp*force_av
        if error > 2 : error = 2
        if error < -2 : error = -2
        linear_travel(serial_handle, str(speed), str(error)) # (s, "speed" , "dist")
        print(abs(error)/speed)
        time.sleep((abs(error)/speed)) # add a bit of settling time for stress relaxation
    time.sleep(2)
    raw_data = loadcell_handle.read(1) # read 1 data point
    force = 452.29*float(raw_data[0]) + 98.155 # magic nums to get force into newtons
    print("Final force ="+str(force))
    time.sleep(2)

def main():
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    # print(avail_devs)
    # device_index = int(input("Which linear actuator device from the list? (eg. index 0 or 1 or 2 or ...):"))
    # s = serial.Serial(avail_devs[device_index],115200,timeout=2) # Connect to port. GRBL operates at 115200 baud
    s = serial.Serial("COM4",115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    init_motion_params(s) # Init Grbl

    ## Setup 34970a connection
    rm = pyvisa.ResourceManager()
    available_devs = rm.list_resources()
    # print(available_devs)
    # device_index = input("Which DAQ device from the list? (eg. index 0 or 1 or 2 or ...):")
    # ohmmeter = rm.open_resource(available_devs[int(device_index)])
    ohmmeter = rm.open_resource(available_devs[3]) # comment if port unknown
    print("Connecting to 34970a DAQ unit...")
    init_ohmmeter_params(ohmmeter)

    ## Setup loadcell connection
    loadcell = nidaqmx.Task()
    init_loadcell_params(loadcell)

    # Initialise lists for storing measurement data
    pos_data = []
    time_data_pos = []

    res_data = []
    time_data_res = []
    avg_time_data_res = []

    force_data = []
    time_data_force = []

    # compression_data = [["current_resistance", "measure_time", "time_taken", "thickness"]]
    # single_sided_electrode_data = [["current_resistance(single electrode)", "measure_time(single electrode)", "time_taken(single electrode"]]
    #
    # # Begin manual testing
    # start_time = time.time() # ref time for manual testing
    # test_specimen_num = input("Input test specimen num:")
    #
    # specimen_thickness = float(input("What is the approx. specimen thickness(in mm)?"))
    # compressed_thickness = input("Enter the compressed test thickness of the test specimen(mm):")
    #
    # compress_test_ans = input("Compression resistivity test?(y or n):")
    # if compress_test_ans == 'y':
    #     # Read compression data
    #     print("Electrode compression test")
    #     print("==========================")
    #     res_input = input("Press enter to read resistance for single-sided electrode test")
    #     for i in range(5):
    #         single_sided_electrode_data.append(read_ohmmeter(ohmmeter,0)) # function ouputs [current_resistance, measure_time, time_taken]
    #         print(" %.4f Ohms in %.4f s" % (single_sided_electrode_data[i+1][0], single_sided_electrode_data[i+1][2]))
    #
    #     num_compr_str = 4
    #     num_compr_res_sample = 5
    #     i_compr_ld = 0
    #     i_compr_uld = 0
    #     res_input = 0
    #
    #     # Loading clamps with compressive strain
    #     while ((res_input != "q") and (i_compr_ld != num_compr_str)):
    #         res_input = input("Press enter to read resistance or 'q' to quit to the next stage")
    #         for i in range(num_compr_res_sample):
    #             compress_n_strain = read_ohmmeter(ohmmeter,0)
    #             compress_n_strain.append(specimen_thickness)
    #             compression_data.append(compress_n_strain)
    #             print(compression_data[(i_compr_ld*num_compr_res_sample)+1+i])
    #         specimen_thickness = specimen_thickness - 0.5 # 0.5mm pitch for M3 bolt
    #         i_compr_ld = i_compr_ld + 1
    #         print(i_compr_ld)
    #         print("Turn each clamp bolt 1 CW revolution")
    #
    #     # Unloading clamps of compressive strain
    #     while ((res_input != "q") and (i_compr_uld != num_compr_str)):
    #         print("Now turn each clamp bolt 1 CCW revolution")
    #         res_input = input("Press enter to read resistance or 'q' to quit to the next stage")
    #         for i in range(num_compr_res_sample):
    #             compress_n_strain = read_ohmmeter(ohmmeter,0)
    #             compress_n_strain.append(specimen_thickness)
    #             compression_data.append(compress_n_strain)
    #             print(compression_data[(i_compr_ld*num_compr_res_sample+i_compr_uld*num_compr_res_sample)+1+i])
    #         specimen_thickness = specimen_thickness + 0.5 # 0.5mm pitch for M3 bolt
    #         i_compr_uld = i_compr_uld + 1
    #         print(i_compr_uld)
    #
    #
    #     format_compr_list = [["Specimen thickness:",specimen_thickness],["Test compressed thickness:",compressed_thickness]]
    #
    #     compression_filename = "compression_resistance_#%s" % (test_specimen_num)
    #     with open(compression_filename+'.csv', 'w', newline='') as csvfile:
    #         data = csv.writer(csvfile, delimiter=',')
    #         print(compression_data)
    #         for i in range(len(compression_data)):
    #             if i <= 1:
    #                 data.writerow(stringify_list(compression_data[i],1)+stringify_list(single_sided_electrode_data[i],1)+stringify_list(format_compr_list[0][i],1)+stringify_list(format_compr_list[1][i],1))
    #             else:
    #                 data.writerow(stringify_list(compression_data[i],1))

    # # TODO: make a parsing error handler for manual jog mode
    # jog_input = input("Enter jog input in mm. Enter 'a' for auto zero cal OR 'q' to manually set zero(start) position:")
    # if jog_input == 'a':
    #     auto_zero_cal(loadcell,s,0.1) # moves stretcher until measured force is zero'd or close enough (f < tol)
    # else:
    #     while (jog_input != "q"):
    #         linear_travel(s, "150", str(jog_input)) # (s, "speed" , "dist")
    #         if jog_input == 'f':
    #             raw_data = loadcell.read(1) # read 1 data point
    #             force = 452.29*float(raw_data[0]) + 98.155
    #             print(force)
    #         jog_input = input("Enter jog dist OR 'f' for force measurement >>")

    # set new zero
    write_g(s,"G90")
    write_g(s,"G10 P0 L20 X0 Y0 Z0")

    # Send desired motion g-code to grbl
    # set_travel_input = input("Set step profile (i.e.[-3,-6,-3,0] stretch material 3mm, 6mm, 3mm ... for strains of 10%, 20%, 10% ...\n>>")
    # set_speed_input = input("How zoomy shall we do the stretchy? [mm/min]: ")
    # linear_travel(s, set_speed_input, set_travel_input) # (s, "speed" , "dist")

    ####################################
    ### Set velocity profile params: ###
    ####################################
    step_profile = [-12,0] # travel 3mm, 6mm ... for strains of 10%, 20% ...
    velocity_profile = [140] # set travel speeds in mm/s
    lag_delay = 40
    ###
    # Begin measurement loop
    ###
    print("Reading data...")
    step_counter = 0
    start_time = time.time() # ref time reset for automated test procedure
    for velocity in velocity_profile:
        for step in step_profile:
            linear_travel(s, velocity, step)
            print("Linear motion set! %dmm @ %dmm/s" % (step,velocity))
            # print("%s" % readbuf_g(s))
            current_pos = 0 # init for while loop condition
            lag_start = 0 # to capture data from just after the strain has stopped
            lag = 0
            while ((float(current_pos) != float(step)) or (lag < lag_delay)):
                if (lag_start == 0) and (float(current_pos) >= float(step)):
                    lag_start = time.time()
                if float(current_pos) >= float(step):
                    lag = time.time() - lag_start

                # Read position
                current_pos = read_pos(s)
                current_time = time.time() - start_time
                pos_data.append(current_pos)
                time_data_pos.append(current_time)

                # Read resistance
                current_res, t_avg, t_d = read_ohmmeter(ohmmeter, start_time)
                res_data.append(current_res)
                avg_time_data_res.append(t_avg)
                time_data_res.append(t_d)

                # Read force
                raw_data = loadcell.read(1) # read 1 data point
                force = 452.29*float(raw_data[0]) + 98.155 # scale data to Newtons (waste of processing time)
                if (force > MAX_LOADCELL_FORCE):
                    raise NameError('Maximum force of {}N for loadcell exceeded'.format(MAX_LOADCELL_FORCE))
                force_data.append(force)
                current_time = time.time() - start_time
                time_data_force.append(current_time)
                print(current_pos)
            step_counter = step_counter + 1
            print("Step complete", step_counter)


    # Write data to CSV file
    filename = input("CSV file name? !Caution will overwrite files without warning!: ")
    write_PosResForce_to_CSV(filename,res_data, avg_time_data_res, pos_data, time_data_pos,
        force_data, time_data_force)

    # Data processing
    # - Best line of fit and R-square values
    # - FFT to see vibrations
    # - Elastic modulus
    # -


    # Plot all data
    plot_q = input("Plot all data?")
    if (plot_q == 'y'):
        fig, (ax1, ax2, ax3) = plt.subplots(3)

        ax1.plot(avg_time_data_res, res_data)
        ax1.set(ylabel='Resistance[Ohm]')

        ax2.plot(time_data_pos, pos_data)
        ax2.set(ylabel='Position[mm]')

        ax3.plot(time_data_force, force_data)
        ax3.set(xlabel='Time[s]', ylabel='Force[N]')

        plt.show()
        plt.savefig(filename)


    # Wait here until grbl is finished to close serial port and file.
    input("Press <Enter> to exit and stop serial Grbl communications.")

    # Close serial port
    s.close()

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
