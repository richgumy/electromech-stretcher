import csv
import matplotlib.pyplot as plt
import pyvisa
import serial
import time

# Setup Serial arduino connection
ardSerial = serial.Serial()

# Setup 34970a connection
rm = pyvisa.ResourceManager()
available_devs = rm.list_resources()
print(available_devs)
device_index = input("Which device from the list (eg. 0 or 1 or 2 or...):")
inst = rm.open_resource(available_devs[int(device_index)])
inst.timeout = 5000
inst.read_termination = '\r'

buf_size = 10
res_vals = []
sample_times = []
time_difs = []
i = 0

start_time = time.time()
mode = input("Press 'm' for manual mode any key for auto mode: ")
if mode == 'm':
    while(i < buf_size):
        ready = input("Press 'y' to continue take next sample: ")
        if ready == 'y':
            t_s = time.time()
            res_vals.append(float((inst.query(":MEASure:RESistance? (@101)")).strip()[1:-4]))
            t_f = time.time()
            t_dif = (t_f + t_s)/2-start_time
            sample_times.append(t_dif)
            i = i + 1
        else:
            break
        ready = ''
else:
        while(i < buf_size):
            t_s = time.time()
            current_res = inst.query(":MEASure:RESistance? (@101)")
            # print(current_res[-1:])
            res_vals.append(float(current_res.strip()[1:-4])*10**float(current_res[-1:]))
            t_f = time.time()
            t_dif = (t_f + t_s)/2-start_time
            sample_times.append(t_dif)
            i = i + 1

with open('data.csv', 'w', newline='') as csvfile:
    data = csv.writer(csvfile, delimiter=',')
    for i in range(len(sample_times)):
        data.writerow([sample_times[i],res_vals[i]])

print(res_vals)
print(sample_times)
plt.plot(sample_times, res_vals,'ro')
plt.xlabel('Time[s]')
plt.ylabel('Resistance[Ohms]')
plt.show()
