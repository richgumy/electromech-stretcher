"""
FILE: keithley2634B.py
AUTHOR: R Ellingham
DATE: Jan 2021
PROGRAM DESC: Driver for communicating withe the Keithley2634B SMU. The program
uses the same Lua commands as specified in the datasheet
but as a class of keithley2634B. For example to call the Lua function
'serial.baud' you can now use Keithley2634B_instance.serial.baud in your code.

Uses pyvisa as to form connection.

TODO:
1 - Make a function which turns class.LuaClass.LuaFunc() into a string
'LuaClass.Luafunc()' to be sent using pyvisa.query or pyvisa.write functions
>> possible implementation use ErrorFlag to turn 'LuaClass.LuaFunc()' From
class.LuaClass.LuaFunc() into a string and then call function to write this to
the pyvisa connection
"""

# System imports
import sys
import numpy as np
import logging # For logging error/warning/other messages during runtime
import time
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
# External imports
import pyvisa

# Setup logger
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

_ch = logging.StreamHandler()
_ch.setFormatter(formatter)

logger = logging.getLogger(__name__)
logger.addHandler(_ch)

# Use log_to_screen in code to print out logger data
def log_to_screen(level: int = logging.DEBUG) -> None:
    log_to_stream(sys.stderr, level)

def log_to_stream(stream_output: Optional[IO], level: int = logging.DEBUG) -> None:
    logger.setLevel(level)
    _ch.setStream(stream_output)
    _ch.setLevel(level)

# Define class for Keithley2634b smu
class K2634b():
    """
    Main class for generating Lua commands
    @param visa_address: Visa address of the instrument.
    @param kwargs: Keyword arguments for Visa connection (e.g. baud_rate/
    timeout)
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

    def __getattr__(self, name: str):
        """Returns the attribute matching passed name."""
        # Get internal dict value matching name.
        value = self.__dict__.get(name)
        if not value: # if 'name' isn't found as a function of the class send it as a Lua function
            # Raise AttributeError if attribute value not found.
            self.connection.write('%s.beep(0.5,500)' % name)
            raise AttributeError()
        # Return attribute value.
        return value

    def multi_getattr(self, attr, default=None):
        return 0

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
        try:
            close_state = self.connection.close
            logger.info("Disconnected successfully from Keithley at %s.", self.visa_address)
        except AttributeError:
            logger.info("Can't disconnect because %s is not connected you silly goose!", self.visa_address)
            self.connection = None
            self.connected = False

log_to_screen() #debug mode
k = K2634b('TCPIP0::169.254.0.1::inst0::INSTR')
k.beeper.beep()

k.disconnect()
