"""
FILE: g-code_sender.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: A script to send g-code commands to an arduino device with GRBL
loaded via a USB serial connection. Specific task for this application are to
send velocity profiles and gather real time data for displacement derviatives.

TODO:
1) Make functions for all desired serial requests/commands
2) Ensure any time delays are accounted for
"""
import csv
import matplotlib.pyplot as plt
import serial
import serial.tools.list_ports
import time



# def DummyFunctionTemplate(IN_PARAMS)
# """
# DESCR: What the function does, very briefly how it works
# IN_PARAMS: Function input parameters
# OUTPUT: Function output value/s
# NOTES: Any issues/limitations with the functionality or use of the function
# EXAMPLES: A couple of examples to show how a range of params can affect output
# """
#   start = function_code_here
#   while (code != complete):
#       do_another_func() # OooOOo this function is doing...
#   return OUTPUT

# Set g-code
g_code_stream = ['$$','?']

def list_devices():
    """
    DESCR: Creates list of the names of all currrently connected serial comport devices
    IN_PARAMS: N/A
    OUTPUT: List of available comports
    NOTES: Requires 'import serial.tools.list_ports'
    """
    avail_devs = serial.tools.list_ports.comports()
    devs_list = []
    for i in range(len(avail_devs)):
        devs_list.append(avail_devs[i].device)
    return devs_list

avail_devs = list_devices()
print(avail_devs)
device_index = int(input("Which device from the list (eg. index 0 or 1 or 2 or ...):"))
s = serial.Serial(avail_devs[device_index],115200) # Connect to port. GRBL operates at 115200 baud

# Wake up grbl
s.write("\r\n\r\n")
time.sleep(2)   # Wait for grbl to initialize
s.flushInput()  # Flush startup text in serial input

# Stream g-code to grbl
for i in len(g_code_stream):
    msg = g_code_stream[i]
    print ('Sending: ' + msg)
    s.write(msg.encode()) # encode msg then send g-code to grbl
    grbl_out = s.readline() # Wait for grbl response with carriage return
    print (' : ' + grbl_out.strip())
    time.sleep(1)

# Wait here until grbl is finished to close serial port and file.
raw_input("  Press <Enter> to exit and disable grbl.")

# Close serial port
s.close()
