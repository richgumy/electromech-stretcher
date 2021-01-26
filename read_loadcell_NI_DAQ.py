
import nidaqmx
from nidaqmx.constants import LineGrouping
from nidaqmx.constants import Edge
from nidaqmx.constants import AcquisitionType
import matplotlib.pyplot as plt
import time
import numpy as np

#keep the plat from closing between data sets
# plt.ion()

i = 0
force_data = []
time_data = []

#calibration factors found through NI MAX: y = 0.456x - 0.364
# electrical_vals= [-2.99, 7.7, 17.9, 23.7, 28]
# physical_vals=[-0.5 ,0, 0.476, 0.75, 0.949]
# Rated output = 1 +/- 0.15 mV/V
# (1mV/V) / (10kg)??

#code portion to read from the NI instrument
task = nidaqmx.Task()
    #adding linear table of values for load cell reading: loadcell SN 666133
# task.ai_channels.add_ai_force_bridge_table_chan("cDAQ2Mod3/ai3",
# name_to_assign_to_channel="10kgLoadcell",
#     min_val=-5.0e-3, max_val=5.0e-3,
#     bridge_config=nidaqmx.constants.BridgeConfiguration.FULL_BRIDGE,
#     voltage_excit_source= nidaqmx.constants.ExcitationSource.INTERNAL,
#     voltage_excit_val=2.5, nominal_bridge_resistance=1000.0,
#     electrical_vals= [0, 1, 2, 3, 4],
#     electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
#     physical_vals=[-1, 0, 1, 2, 3],
# physical_units=nidaqmx.constants.BridgePhysicalUnits.NEWTONS,
#     custom_scale_name=None)

task.ai_channels.add_ai_force_bridge_two_point_lin_chan('cDAQ2Mod3/ai3',
name_to_assign_to_channel='10kgLoadcell',
min_val=-2.0, max_val=2.0,
units=nidaqmx.constants.ForceUnits.NEWTONS,
voltage_excit_val=10,
nominal_bridge_resistance=1000.0,
first_electrical_val=0.0483, second_electrical_val=0.121,
electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
first_physical_val=0.0, second_physical_val=0.007522,
physical_units=nidaqmx.constants.BridgePhysicalUnits.NEWTONS)
## 4.7N offset from zero i.e. F = v*x - 4.7
task.timing.cfg_samp_clk_timing(50000, source="",
active_edge=Edge.RISING, sample_mode=AcquisitionType.FINITE)
print("Loadcell configuration = successful!\n")
start_time = time.time()
raw_data = task.read(1000) # Caps out out 1000 samples for unkown reason?
# while i<100:
#     raw_data = task.read(500)
#     # print(raw_data)
#     avg_sum = 0
#     for sample in range(len(raw_data)):
#         avg_sum = avg_sum + raw_data[sample]
#     avg = avg_sum/len(raw_data)
#     force = 452.29*float(raw_data[0]) + 98.155 # cal 10/2020
#     # # cal 14/12/2020 -> 460*float(raw_data[0]) + 100.02
#     force_data.append(force)
#     time_i = time.time() - start_time
#     time_data.append(time_i)
#     # #plt.scatter(i,data[1],c='b')
#     # # plt.pause(0.01)
#     i=i+1
#     print(force)
    # time.sleep(0.1)
force_data = np.array(raw_data) - 4.7
time_tot = time.time() - start_time
time_lin = np.linspace(0,0.02,1000)
print("total time: %.3f" % time_tot)
plt.scatter(time_lin,force_data,c='r')
# plt.scatter(time_data,force_data,c='r')
plt.show()
