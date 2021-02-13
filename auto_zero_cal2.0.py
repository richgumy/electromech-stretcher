"""
FILE: auto_zero_cal.py
AUTHOR: R Ellingham
DATE: Oct 2020
PROGRAM DESC: Basic proportional controller to calibrate zero force/stress/strains
position for the stretch-n-yeet5000 device.

Parameters measured:
-> Strain and Force

TODO:
1)
"""
print("importing some thicc libraries...")

from stretch_n_measure_smu import *

print("Thicc libraries loaded.")

def auto_zero_crtl(loadcell_handle, serial_handle, tolerance):
    """
    DESCR: Rough PI controller to get stress in test rig to zero
    IN_PARAMS: loadcell_handle, serial_handle, zero tolerance [N]
    NOTES: Not well tested.
    TODO:
    """
    error = tolerance*2 # start error bigger than tolerance
    error_prev = 0
    buf_size = 40
    f_buf = [0]*buf_size
    speed = 180
    Kp = 2 # control P gain constant
    Ki = 0.2 # control I gain constant
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

def main(tolerance): # Tolerance = How close to zero force/stress is close enough?
    ## Setup grbl serial coms:
    avail_devs = list_serial_devices()
    print(avail_devs)
    dev = input("What is the index of the desired COM port? ")
    s = serial.Serial(avail_devs[int(dev)],115200,timeout=2) #comment this and uncomment above for interactive choice of com port
    print("Connecting to grbl device...")
    init_motion_params(s) # Init Grbl

    ## Setup loadcell connection
    loadcell = nidaqmx.Task()
    init_loadcell_params(loadcell)

    auto_zero_crtl(loadcell, s, tolerance)

if __name__ == "__main__":
	# Default input parameters
	tolerance = 0.005
	import sys
	if len(sys.argv)>1: tolerance   = (sys.argv[1])
	# if len(sys.argv)>2: input2   =int(sys.argv[2])

	main(tolerance)
