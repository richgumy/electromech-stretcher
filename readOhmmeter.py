import pyvisa
import time
start_time = time.time()

rm = pyvisa.ResourceManager()
available_devs = rm.list_resources()
inst = rm.open_resource(available_devs[0])
# print(inst.query("*IDN?"))
buf_size = 5
res_vals = []
sample_times = []
i = 0

while(i < buf_size):
    t_s = time.time()
    res_vals.append(inst.query(":MEASure:RESistance? (@101)"))
    t_f = time.time()
    t_dif = t_f - t_s
    sample_times.append(t_dif)
    i = i + 1

# sum_res = 0
# sum_time = 0
# for i in range(len(sample_times)):
#     sum_res = sum_res + int(res_vals[i])
#     sum_time = sum_time + int(sample_times[i])
#
# average_res = sum_res/len(res_vals)
# average_time = sum_time/len(sample_times)

print(res_vals)
print(sample_times)
# print(average_res)
# print(average_time)
