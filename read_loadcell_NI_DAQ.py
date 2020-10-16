
import nidaqmx
from nidaqmx.constants import LineGrouping
from nidaqmx.constants import Edge
from nidaqmx.constants import AcquisitionType
import matplotlib.pyplot as plt

#keep the plat from closing between data sets
plt.ion()

i = 0

#calibration factors found through NI MAX: y = 0.456x - 0.364

#code portion to read from the NI instrument
with nidaqmx.Task() as task:
    #adding linear table of values for load cell reading: loadcell SN 666133
    task.ai_channels.add_ai_force_bridge_table_chan("cDAQ2Mod3/ai3",
 name_to_assign_to_channel="10kgLoadcell",
        min_val=-5.0e-3, max_val=5.0e-3,
        bridge_config=nidaqmx.constants.BridgeConfiguration.FULL_BRIDGE,
        voltage_excit_source= nidaqmx.constants.ExcitationSource.INTERNAL,
        voltage_excit_val=2.5, nominal_bridge_resistance=1000.0,
        electrical_vals= [-2.99, 7.7, 17.9, 23.7, 28],
        electrical_units=nidaqmx.constants.BridgeElectricalUnits.M_VOLTS_PER_VOLT,
        physical_vals=[-0.5 ,0, 0.476, 0.75, 0.949],
physical_units=nidaqmx.constants.BridgePhysicalUnits.KILOGRAM_FORCE,
        custom_scale_name=None)
    task.timing.cfg_samp_clk_timing(100, source="",
    active_edge=Edge.RISING, sample_mode=AcquisitionType.FINITE, samps_per_chan=30)
    print("Loadcell configuration = successful!\n")
    while i<50:
        data = task.read(30)
        print("Data received!")
        plt.scatter(i,data[0],c='r')
        #plt.scatter(i,data[1],c='b')
        plt.pause(0.01)
        plt.show()
        i=i+1
        print(data)
plt.show()
