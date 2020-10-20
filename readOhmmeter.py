import csv
import matplotlib.pyplot as plt
import pyvisa
import serial
import time


# Setup 34970a connection
rm = pyvisa.ResourceManager()
available_devs = rm.list_resources()
print(available_devs)
device_index = input("Which device from the list (eg. index 0 or 1 or 2 or ...):")
inst = rm.open_resource(available_devs[int(device_index)])
# inst = rm.open_resource(available_devs[2]) # comment if port unknown
inst.baud_rate = 115200
inst.timeout = 5000
# inst.read_termination = '\r'
inst.write(":CONF:RES 10000000,1000,(@101)") 
# print("Resistance range set to 10MOhm, 1kOhm resolution")
# print(inst.query(":CONF? (@101)"))

inst.write(":ROUT:CHAN:DEL 0,(@101)")
# print("delay set 0")
# print(inst.query(":ROUT:CHAN:DEL? (@101)"))
# print(inst.query(":CONF? (@101)"))
# inst.read_termination = '\r'
inst.write(":SENS:RES:NPLC 1,(@101)") #set PLC(i.e. 50Hz power line cycle for NZ) to 0.02 for 400us integration time... or higher for better noise suppression
print("PLC set to %.2f" % float(inst.query(":SENS:RES:NPLC? (@101)")))


inst.write("ROUT:SCAN (@101)") # Set scan channel
inst.write(":TRIG:COUNT 1") # Set number of measurements taken each scan

buf_size = 100
res_vals = []
avg_sample_times = [] #
sample_times = [] # time taken for sample
time_difs = []
current_res = '0' # init start resistance to zero string
i = 0

mode = input("Press 'm' for manual mode any key for auto mode: ")
if mode == 'm':
    while(i < buf_size):
        ready = input("Press 'y' to continue to take next sample: ")
        start_time = time.time()
        if ready == 'y':
            t_s = time.time()
            res_vals.append(float((inst.query(":MEAS:RES? (@101)")).strip()[1:-4]))
            t_f = time.time()
            # t_dif = (t_f + t_s)/2-start_time
            sample_times.append(t_s)
            i = i + 1
        else:
            break
        ready = ''
else:
        start_time = time.time() # Record start time to make measurements relative not absolute
        while(i < buf_size):
            t_s = time.time()
            inst.write(":INIT") # Begin scan
            time.sleep(0.002) # Wait for scan to finish
            while len(current_res) <= 4: # If scan not finished keep questioning until result appears
                current_res = inst.query("R?")
            # print(current_res)
            res_vals.append(current_res)
            t_f = time.time()
            t_avg = (t_f + t_s)/2-start_time # Record time of measurement halfway between when measurement was requested and when it was received
            avg_sample_times.append(t_avg)
            t_d = t_f - t_s # How long did it take to receive the message from time of request
            sample_times.append(t_d)
            i = i + 1
print("Delay\nRange and res\nPLC..")
print(inst.query(":ROUT:CHAN:DEL? (@101)"))
print(inst.query(":CONF? (@101)"))
print(inst.query(":SENS:RES:NPLC? (@101)"))

with open('data.csv', 'w', newline='') as csvfile:
    data = csv.writer(csvfile, delimiter=',')
    for i in range(len(sample_times)):
        res_vals[i] = float(res_vals[i].strip()[5:-4])*10**float(res_vals[i].strip()[-1:])
        data.writerow([sample_times[i],res_vals[i]])

print(res_vals)
print(sample_times)
print(avg_sample_times)
plt.plot(avg_sample_times, res_vals,'ro')
plt.xlabel('Time[s]')
plt.ylabel('Resistance[Ohms]')
plt.show()
