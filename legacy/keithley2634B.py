e"""
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
def log_to_screen(stream_output: Optional[IO] = sys.stderr, level: int = logging.DEBUG) -> None:
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
        # try:
        #     # check if parsed function is a valid Lua commands
        #     # valid_cmd = check_lua_table(name) # TODO
        #     # logger.debug("%s is a valid function", name)
        #     print("not complete yet")
        # except AttributeError:
        # check if parsed function is a valid Lua commands
        # valid_cmd = check_lua_table(name) # TODO
        # logger.debug("%s is a valid function", name)
        logger.info("Attempting Lua cmd '%s'.", name)
        if not value: # if 'name' isn't found as a function of the class send it as a Lua function
            # Raise AttributeError if attribute value not found.
            # self.connection.write('%s.beep(0.5,500)' % name)
            # self.connection.write("display.clear()")
            # self.connection.write("display.settext('Running all the code')")
            # self.connection.write('smua.abort()')
            raise AttributeError('hotdawg')
        # Return attribute value.
        return value

    def multi_getattr(self, attr, default=None):

        return 0

    def starWars(self):
        """
        Plays Star Wars main theme song
        """
        song_notes = [('d1',0.33),('d1',0.33),('d1',0.33),('g1',2),('d2',2),('c2',0.33),('b2',0.33)
        ,('a2',0.33),('g2',2),('d2',1),('c2',0.33),('b2',0.33),('a2',0.33),('g2',2),('d2',1),('c2',0.33),
        ('b2',0.33),('c2',0.33),('a2',2),('d1',0.33),('d1',0.33),('d1',0.33),('g1',2),('d2',2),('c2',0.33),('b2',0.33)
        ,('a2',0.33),('g2',2),('d2',1),('c2',0.33),('b2',0.33),('a2',0.33),('g2',2),('d2',1),('c2',0.33),
        ('b2',0.33),('c2',0.33),('a2',2),('d1',0.75),('d1',0.25),('e1',1.5),('e1',0.5),('c2',0.5),('b2',0.5),
        ('a2',0.5),('g1',0.5),('g1',0.33),('a2',0.33),('b2',0.33)] # Good enough lol
        # [(NoteOctave,Length),..]
        for note in song_notes:
            self.connection.write('beeper.beep(%.2f,%.2f)' % (note[1]/2, self.getFrequency(note[0])))

    def getFrequency(self, note, A4=440):
        """
        Gets frequency for a given note and octave in format "NoteOctave"
        """
        notes = ['a', 'a#', 'b', 'c', 'c#', 'd', 'd#', 'e', 'f', 'f#', 'g', 'g#']

        octave = int(note[2]) if len(note) == 3 else int(note[1])

        keyNumber = notes.index(note[0:-1]);

        if (keyNumber < 3) :
            keyNumber = keyNumber + 12 + ((octave - 1) * 12) + 1;
        else:
            keyNumber = keyNumber + ((octave - 1) * 12) + 1;

        return A4 * 2** ((keyNumber- 49) / 12)

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

log_to_screen() # debug mode
k = K2634b('TCPIP0::169.254.0.1::inst0::INSTR')
print(k)

k.disconnect()
