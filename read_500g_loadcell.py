
import nidaqmx
from nidaqmx.constants import LineGrouping
from nidaqmx.constants import Edge
from nidaqmx.constants import AcquisitionType
import matplotlib.pyplot as plt
import time
import numpy as np

i = 0
force_data = []
time_data = []


def average_list(in_list):
    sum = 0
    for item in in_list:
        sum = sum + float(item)
    return sum/len(in_list)


#code portion to read from the NI instrument
task = nidaqmx.Task()

task.ai_channels.add_ai_force_bridge_two_point_lin_chan('cDAQ2Mod3/ai3',
name_to_assign_to_channel='500gLoadcell',
min_val=-5.0, max_val=5.0,
units=nidaqmx.constants.ForceUnits.NEWTONS,
voltage_excit_val=6,
bridge_config=nidaqmx.constants.BridgeConfiguration.FULL_BRIDGE,
nominal_bridge_resistance=1000.0,
first_electrical_val=0.0079, second_electrical_val=-0.683,
electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
first_physical_val=0.0, second_physical_val=4.905,
physical_units=nidaqmx.constants.BridgePhysicalUnits.NEWTONS)
task.timing.cfg_samp_clk_timing(1000,samps_per_chan=10)
# task.ai_channel(ai_lowpass_enable=True, ai_lowpass_cutoff_freq=100)
print("Loadcell configuration = successful!\n")
start_time = time.time()
while True:
    raw_data = task.read(10) # Caps out out 1000 samples for unkown reason?
    force_avg = average_list(raw_data)
    print(force_avg)
    time.sleep(1)

force_data = np.array(raw_data)
time_tot = time.time() - start_time
time_lin = np.linspace(0,0.1,1000)
print("total time: %.3f" % time_tot)
plt.scatter(time_lin,force_data,c='r')
# plt.scatter(time_data,force_data,c='r')
plt.show()
