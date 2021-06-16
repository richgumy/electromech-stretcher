"""
FILE: linact_grbl.py
AUTHOR: R Ellingham
DATE MODIFIED: June 2021
PROGRAM DESC: Library for a the linear actuator functions for the electromech stretcher
project
NOTES: Requires application code to have py serial library and make a serial instance/s
DATE CREATED: June 2021

TODO:
"""

import serial
import time
import re


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