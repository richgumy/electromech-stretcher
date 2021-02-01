"""
FILE: k2600.py
AUTHOR: R Ellingham
DATE: Jan 2021
PROGRAM DESC: Driver for communicating withe the Keithley2634B SMU. The program
has hard coded specific Lua commands specified in the datasheet
but as a class of K2600. For example to call the Lua function
'serial.baud' you can now use K2600_instance.serial.baud in your code.

Another version of this driver is yet to be made where any function can be
assumed as an attribute of K2600. For example K2600.function is not explicitly
declared within the code but will be converted to a string, verfied as a valid
Lua command and then sent to the SMU via a pyvisa write or query command.

Uses pyvisa as to form connection.

TODO:
1 - Hard code listed functions
2 - Make a query response handler function
"""
import sys
import pyvisa
import numpy as np
import time

import logging
from typing import (
    IO,
    Optional,
    Any,
    Dict,
    Union,
    List,
    Tuple,
    Sequence,
    Iterable,
    Iterator,
)

# Setup logger
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

_ch = logging.StreamHandler()
_ch.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.addHandler(_ch)

def debug_enable(stream_output: Optional[IO] = sys.stderr, level: int = logging.DEBUG) -> None:
    """
    Displays verbose logger ouput in console
    """
    logger.setLevel(level)
    _ch.setStream(stream_output)
    _ch.setLevel(level)

class K2600:
    """
    Contains some useful Lua commands:(add more as required)
    """
    def __init__(
        self,
        visa_address: str,
        raise_keithley_errors: bool = False,
        **kwargs,
        ) -> None:
        self.visa_address = visa_address
        self._connection_kwargs = kwargs

        self.raise_keithley_errors = raise_keithley_errors

        self.rm = pyvisa.ResourceManager()

        self.connect(**kwargs)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}({self.visa_address})"

    def connect(self, **kwargs) -> bool:
        """
        Connects to Keithley device
        @param kwargs: Keyword arguments for Visa connection (e.g. baud_rate/
        timeout)
        @return Whether connection has succeeded
        """
        try:
            self.connection = self.rm.open_resource(self.visa_address,**kwargs)
            self.connection.read_termination = "\n"
            self.connected = True
            logger.debug("Connected to Keithley at %s.", self.visa_address)
        except ValueError:
            self.connection = None
            self.connected = False
            raise
        except ConnectionError:
            logger.info(
                "Connection error. Please check that no other program is connected."
            )
            self.connection = None
            self.connected = False
        except AttributeError:
            logger.info("Invalid VISA address %s.", self.visa_address)
            self.connection = None
            self.connected = False
        except Exception:
            logger.info("Could not connect to Keithley at %s.", self.visa_address)
            self.connection = None
            self.connected = False

        return self.connected

    def disconnect(self):
        """
        Disconnects Keithley device if present
        @return Connection status
        """
        try:
            close_state = self.connection.close
            logger.info("Disconnected successfully from Keithley at %s.", self.visa_address)
        except AttributeError:
            logger.info("Can't disconnect because %s is not connected you silly goose!", self.visa_address)
            self.connection = None
            self.connected = False
        return self.connected

    ### Lua functions as classes and functions ###
    class beeper:
        def enable(self, state):
            self.connection.write("beeper.enable('%d')" % state)
        def beep(self, length=1, freq=1000):
            self.connection.write("beeper.beep(%.2f,%.2f)" % (length, freq))

    def delay(self, time):
        self.connection.write("delay('%d')" % time)

    class display:
        def clear(self):
            self.connection.write("display.clear()")
        def settext(self, text: str):
            self.connection.write("display.settext('%s')" % text)
            print("display.settext('%s')" % text)
        class smua:
            class measure:
                func = 0

    def makegetter(self):
        self.connection.write("makegetter()")

    def makegetter(self):
        self.connection.write("makegetter()")

    def print(self, value):
        self.connection.write("print(value)")

    def PulseIMeasureV(self):
        self.connection.write("PulseIMeasureV()")

    class os:
        def time(self):
            self.connection.write("os.time()")

    class smua:
        def makebuffer(self):
            self.connection.write("smua.makebuffer()")
        class measure:
            def v(self):
                self.connection.write("smua.measure.v()")
            def i(self):
                self.connection.write("smua.measure.i()")
            def r(self):
                self.connection.write("smua.measure.r()")
            def r(self):
                self.connection.write("smua.measure.r()")
            class filter:
                def enable(self, value=None):
                    if value != None:
                        self.connection.write("smua.measure.filter.enable = %d", value)
                    else:
                        return = self.connection.query("smua.measure.filter.enable")
                def type(self, value=None):
                    if value != None:
                        self.connection.write("smua.measure.filter.type = %d", value)
                    else:
                        return = self.connection.query("smua.measure.filter.type")
        def reset(self):
            self.connection.write("smua.reset()")

        def savebuffer(self):
            self.connection.write("smua.savebuffer()")

        def sense(self, value=None):
            if value != None:
                self.connection.write("smua.sense = %d", value)
            else:
                return = self.connection.query("smua.sense")

        class source:
            def levelv(self, value=None):
                if value != None:
                    self.connection.write("smua.source.levelv = %d", value)
                else:
                    return = self.connection.query("smua.source.levelv")
            def leveli(self, value=None):
                if value != None:
                    self.connection.write("smua.source.leveli = %d", value)
                else:
                    return = self.connection.query("smua.source.leveli")
            def limitv(self, value=None):
                if value != None:
                    self.connection.write("smua.source.limitv = %d", value)
                else:
                    return = self.connection.query("smua.source.limitv")
            def limiti(self, value=None):
                if value != None:
                    self.connection.write("smua.source.limiti = %d", value)
                else:
                    return = self.connection.query("smua.source.limiti")
            def output(self, value=None):
                if value != None:
                    self.connection.write("smua.source.output = %d", value)
                else:
                    return = self.connection.query("smua.source.output")

    class smub:
        def makebuffer(self):
            self.connection.write("smub.makebuffer()")
        class measure:
            def v(self):
                self.connection.write("smub.measure.v()")
            def i(self):
                self.connection.write("smub.measure.i()")
            def r(self):
                self.connection.write("smub.measure.r()")
            def r(self):
                self.connection.write("smub.measure.r()")
            class filter:
                def enable(self, value=None):
                    if value != None:
                        self.connection.write("smub.measure.filter.enable = %d", value)
                    else:
                        return = self.connection.query("smub.measure.filter.enable")
                def type(self, value=None):
                    if value != None:
                        self.connection.write("smub.measure.filter.type = %d", value)
                    else:
                        return = self.connection.query("smub.measure.filter.type")
        def reset(self):
            self.connection.write("smub.reset()")

        def savebuffer(self):
            self.connection.write("smub.savebuffer()")

        def sense(self, value=None):
            if value != None:
                self.connection.write("smub.sense = %d", value)
            else:
                return = self.connection.query("smub.sense")

        class source:
            def levelv(self, value=None):
                if value != None:
                    self.connection.write("smub.source.levelv = %d", value)
                else:
                    return = self.connection.query("smub.source.levelv")
            def leveli(self, value=None):
                if value != None:
                    self.connection.write("smub.source.leveli = %d", value)
                else:
                    return = self.connection.query("smub.source.leveli")
            def limitv(self, value=None):
                if value != None:
                    self.connection.write("smub.source.limitv = %d", value)
                else:
                    return = self.connection.query("smub.source.limitv")
            def limiti(self, value=None):
                if value != None:
                    self.connection.write("smub.source.limiti = %d", value)
                else:
                    return = self.connection.query("smub.source.limiti")
            def output(self, value=None):
                if value != None:
                    self.connection.write("smub.source.output = %d", value)
                else:
                    return = self.connection.query("smub.source.output")

    class timer:
        def reset(self):
            self.connection.query("timer.reset()")
        class measure:
            def t(self):
                return self.connection.query("timer.measure.t()") # time since reset


debug_enable()

k = K2600('TCPIP0::169.254.0.1::inst0::INSTR')

# k.display.clear(k)
k.display.settext(k,"hello")

k.disconnect()
