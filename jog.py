"""
FILE: jog.py
AUTHOR: R Ellingham
DATE: Nov 2020
PROGRAM DESC: Program used to jog a linear actuator Grbl device back and forth
using serial g-code commands
"""

import serial
import serial.tools.list_ports
import time


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
    NOTES:    Requires serial library and Grbl device attached
              Requires serial_handle e.g. serial_handle = serial.Serial("COM4",115200)
    """
    serial_handle.flushInput()
    #set speed and displacement
    write_g(serial_handle,"G21G1"+"Z"+str(displacement)+"F"+str(lin_speed))
    return


def main():
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    print(avail_devs)
    device_index = int(input("Which linear actuator device from the list? (eg. index 0 or 1 or 2 or ...):"))
    s = serial.Serial(avail_devs[device_index],115200,timeout=2) # Connect to port. GRBL operates at 115200 baud
    # s = serial.Serial("COM8",115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    init_motion_params(s) # Init Grbl
    jog_x = input("Enter jog displacement in mm, x:")
    jog_v = input("Enter jog speed in mm/s, v:")
    while(jog_x != 'q'):
        linear_travel(s, str(jog_v), str(jog_x)) # (s, "speed" , "dist")
        jog_x = input("x:")
        jog_v = input("v:")
    s.close()

    return 0

main()
