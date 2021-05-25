"""
FILE: c_meter_350610.py
AUTHOR: R Ellingham
DATE MODIFIED: May 2021
PROGRAM DESC: A library of commonly used visa commands for the Hioki 350610 C-meter
(capacitance) measurment device.

Uses pyvisa to form connection, and read & write commands

Some common commands:

*IDN? - Name of device ('*' char means universal across IEEE488.2 compliant devices)
*CLS - clear event register
*RST - Re-initialises device
*TRG - Perform sampling once
*(MORE register cmds available) - ...

:AVER - set number of measurements to average over (NR1)
:CAL:CAB - set the cable length in meters 0,1,2 m (NR1)
:RANGE - set range of expected capacitance measurement in farads (1kHz=9-24,1MHz=1-12) (NR1)
:FREQ - set measurement freq to 1E3 or 1E6 [kHz] (NR3)
:FREQ:SHIFT - shifts the measurement frequency from -2% to +2% (e.g. 1=>1%=>1.01MHz, only works for 1MHz freq) (NR1)
:LEV? - set the measurement signal voltage value to 1 or 0.5 volts (NR2)
:MEAS? - what is the current capacitance value? outputs in format: <Measurement Status (NR1)>,C <Measurement Value (NR3)>, D
(Q) <Measurement Value (NR2)>,<Panel Load Number (NR1)> (e.g. '2,  0.0549E-09, 999999,0\n')


*
note: input values are formatted as NRf:
NR1 = integer (+1,-23, 34)
NR2 = fixed point data (+1.23, -23.45, 34.56)
NR3 = float exp representation (+1.0E-2, -2.3E+4)
*


TODO:
- Make a class
"""

import pyvisa
import time
import numpy as np

class Cmeter350610:
    def __init__(self, port_auto=1, freq=1E3, cct_equiv='PAR', **kwargs):
        """
        DESCR: Initialises Hioki 3506-10 c_meter object
        IN_PARAMS: automatically choose device from port or manually select (1 or 0),
        measurement frequency, cct_equiv (SER=series, PAR=parallel)
        OUTPUT: N/A
        NOTES: defaults set for all inputs
        TODO: Add initialisation error
        """
        # init class params
        self.port_auto = port_auto
        self.frequency = freq
        self.cct_equiv = cct_equiv
        self._connection_kwargs = kwargs
        self.rm = pyvisa.ResourceManager()

        self.connect(**kwargs)

        self.connection.write("*IDN?")
        print("Device connected: " + self.connection.read().rstrip())
        # Set device params
        self.connection.write(":CIRC "+str(cct_equiv))
        self.connection.write(":FREQ "+str(freq))
        print("frequency and equivalent circuit values set.")

    def connect(self, **kwargs) -> bool:
        """
        Connects to a Cmeter 3506-10 device
        @param kwargs: Keyword arguments for Visa connection (e.g. baud_rate/
        timeout)
        @return Whether connection has succeeded
        """
        try:
            available_devs = self.rm.list_resources()
            # create visa connection automatically...
            self.connection = self.rm.open_resource(available_devs[-1])
            print("connecting to default port..")
            # if self.auto:
            #     self.connection = self.rm.open_resource(available_devs[-1])
            #     print("connecting to default port..")
            # # ...or manually choose port
            # else:
            #     print(available_devs)
            #     index = input("Choose index of desired device:")
            #     self.connection = rm.open_resource(available_devs[index])
            self.connection.read_termination = "\n"
            self.connected = True
            # logger.debug("Connected to Cmeter")

        except ValueError:
            self.connection = None
            self.connected = False
            raise
        except ConnectionError:
            # logger.info(
            #     "Connection error. Please check that no other program is connected."
            # )
            self.connection = None
            self.connected = False
        except AttributeError:
            # logger.info("Invalid VISA address...")
            self.connection = None
            self.connected = False
        except Exception:
            # logger.info("Could not connect to Cmeter at given address")
            self.connection = None
            self.connected = False

        return self.connected

    def disconnect(self):
        """
        Disconnects Cmeter device if present
        @return Connection status
        """
        try:
            close_state = self.connection.close
            # logger.info("Disconnected successfully from Cmeter")
        except AttributeError:
            # logger.info("Can't disconnect because Cmeter is not connected you silly goose!")
            self.connection = None
            self.connected = False
        return self.connected

    def get_measurement(self):
        self.connection.write("MEAS?")
        c_reading_raw = self.connection.read().split(',')
        c_reading = float(c_reading_raw[1])
        d_reading = float(c_reading_raw[2])
        print(d_reading)
        ## TODO: if d_reading == 999999 (or whatever the error number is): ESR = 0
        ESR = d_reading * (2*np.pi*self.frequency*c_reading)
        return c_reading, ESR

# cmeter = init()
#
# while(1):
#     t_s = time.time()
#     capac = get_measurement(cmeter)
#     t_f = time.time()
#     print((t_f - t_s)*1000)
#     print(capac)
#     time.sleep(0.1)
