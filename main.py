"""
FILE: g-main.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: Gather data from a stretched conductive elastomer in real time using
serial communication. Writing the data to a CSV file ready for analysis

Parameters measured:
-> Resistance, Strain, Force, and Time (for each individual measurement)

TODO:
1) Make functions for all desired serial requests/commands
2) Ensure all significant time delays are accounted-for/read
3)
"""
import csv
import matplotlib.pyplot as plt
import serial
import serial.tools.list_ports
import time
import re
import pyvisa

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
    NOTES:    Requires serial library
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """
    serial_handle.flushInput()
    #set speed and displacement
    write_g(serial_handle,"G21G1"+"Z"+displacement+"F"+lin_speed)
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
    ohmmeter.baud_rate = 115200
    ohmmeter.timeout = 5000

    ohmmeter.write(":CONF:RES 10000000,1000,(@101)") #set range and resolution of resistance measurements
    ohmmeter.write(":ROUT:CHAN:DEL 0,(@101)") #set delay between scans
    ohmmeter.write(":SENS:RES:NPLC 1,(@101)") #set PLC(i.e. 50Hz power line cycle for NZ) to 0.02 for 400us integration time... or higher for better noise suppression
    ohmmeter.write("ROUT:SCAN (@101)") # Set scan channel
    ohmmeter.write(":TRIG:COUNT 1") # Set number of measurements taken each scan

    return 0



### Load cell data acquisition functions:




## Data processing fucntions:
def write_all_to_CSV(resistance, time_R, displacement, time_d, force, time_f):
    with open('data.csv', 'w', newline='') as csvfile:
        data = csv.writer(csvfile, delimiter=',')
        for i in range(len(time_R)):
            data.writerow([resistance[i], time_R[i], displacement[i], time_d[i], force[i], time_f[i])
    return 0

def main():
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    print(avail_devs)
    device_index = int(input("Which linear actuator device from the list? (eg. index 0 or 1 or 2 or ...):"))
    s = serial.Serial(avail_devs[device_index],115200,timeout=2) # Connect to port. GRBL operates at 115200 baud
    # s = serial.Serial("COM4",115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    init_motion_params(s) # Init Grbl

    ## Setup 34970a connection
    rm = pyvisa.ResourceManager()
    available_devs = rm.list_resources()
    print(available_devs)
    device_index = input("Which DAQ device from the list? (eg. index 0 or 1 or 2 or ...):")
    ohmmeter = rm.open_resource(available_devs[int(device_index)])
    # ohmmeter = rm.open_resource(available_devs[2]) # comment if port unknown
    print("Connecting to 34970a DAQ unit...")
    init_ohmmeter_params(ohmmeter)


    print("Sending motion command...")
    # Send desired motion g-code to grbl
    linear_travel(s, "150", "10")

    print("Reading position data...")
    pos_data = []
    time_data = []
    start_time = time.time()
    for i in range(100):
        # Read position
        current_pos = read_pos(s)
        current_time = time.time() - start_time
        pos_data.append(current_pos)
        time_data.append(current_time)

        # Read resistance
        
        # Read force

        time.sleep(0.01)
    print(pos_data)
    plt.plot(time_data, pos_data,'ro')
    plt.yticks(range(1,int(max(pos_data)+max(pos_data)/10),int(max(pos_data)/10)))
    plt.show()
    # Wait here until grbl is finished to close serial port and file.
    input("Press <Enter> to exit and stop serial communications.")

    # Close serial port
    s.close()

    return 0

main()
