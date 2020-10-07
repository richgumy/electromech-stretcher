import csv
import pyvisa
import time
import matplotlib.pyplot as plt

rm = pyvisa.ResourceManager()
available_devs = rm.list_resources()
inst = rm.open_resource(available_devs[0])
# print(inst.query("*IDN?"))
buf_size = 10
res_vals = []
sample_times = []
time_difs = []
i = 0

start_time = time.time()

while(i < buf_size):
    t_s = time.time()
    res_vals.append(float((inst.query(":MEASure:RESistance? (@101)")).strip()[1:-4]))
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
plt.ylabel('Resistance[M Ohms]')
plt.show()
