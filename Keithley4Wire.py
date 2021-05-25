# Test 4-wire measurement using both channels of the Keithley SMU
from time import sleep
import visa

smu_resource = "TCPIP0::169.254.0.1::inst0::INSTR"

source_I = 10e-6  # Current to attempt to source
limit_V = 10  # Max voltage to apply doing so

rm = visa.ResourceManager("@py")

# Connect to SMU
smu = rm.open_resource(smu_resource)

# Write a list of commands to the SMU, saving typing smu.write() so much
def smu_write_list(command_list):
    for command in command_list:
        smu.write(command)

def list_to_newlines_string(command_list):
    new_string = ""
    for command in command_list:
        new_string = new_string+command+'\n'
    return new_string

# Avoids confusion about initial state and if errors are from
# this program run or the last one
smu_initialisation = ["reset()", "errorqueue.clear()"]
smu_write_list(smu_initialisation)

# Try and clear out any rubbish data from buffers somewhere
# This seems to happen if the SMU was previously vomiting data
smu.timeout = 2000
try:
    smu.read_raw()
except:
    pass
print("SMU: " + smu.query("*IDN?").rstrip())

# Looks like it is autorange by default
# print(smu.query("print(smua.source.autorangev)"))
# print(smu.query("print(smua.source.autorangei)"))

# Looks like it is 2 wire per channel by default
# print(smu.query("print(smua.sense==smua.SENSE_LOCAL)"))
# print(smu.query("print(smua.sense==smua.SENSE_REMOTE)"))

# Commands to set up Channel A which is sourcing current
# smu_setup_a = [
#     "smua.source.func=smua.OUTPUT_DCAMPS",
#     "smua.source.leveli={:.5e}".format(source_I),
#     "smua.source.limitv={:.5e}".format(limit_V),
#     "smua.measure.nplc=1",
#     "smua.source.output = smua.OUTPUT_ON",
# ]
# smu_setup_a = [
#     "smua.reset()"
#     "smua.source.rangev = 5"
#     "smua.source.rangei = 1"
#     "smua.source.levelv = 0"
#     "smua.measure.rangev = 5"
#     "smua.measure.rangei = 1"
#     "smua.measure.nplc = 0.01"
#     "smua.measure.autozero = smua.AUTOZERO_ONCE\n"
#     "smua.nvbuffer1.clear()"
#     "smua.nvbuffer1.appendmode = 1"
#     "smua.source.output = smua.OUTPUT_ON\n"
#     "ConfigPulseVMeasureI(smua, 0, 5, 1, 0.002, 0.2, 10, smua.nvbuffer1, 1)"
#     "f2, msg2 = InitiatePulseTest(1)"
# ]

smu_setup_a = [
    "period_timer = trigger.timer[1]"
    "pulse_timer = trigger.timer[2]"
    "smua.trigger.source.listi({-10e-6,10e-6})"
    "smua.trigger.source.action = smua.ENABLE\n"
    "smua.source.rangei = 0.1\n"
    "smua.trigger.measure.action = smua.DISABLE\n"
    "pulse_timer.delay = 0.0006\n"
    "pulse_timer.stimulus = period_timer.EVENT_ID\n"
    "pulse_timer.count = 1\n"
    "period_timer.delay = 0.005\n"
    "period_timer.count = 9\n"
    "period_timer.stimulus = smua.trigger.SWEEPING_EVENT_ID\n"
    "period_timer.passthrough = true\n"
    "smua.trigger.source.stimulus = period_timer.EVENT_ID\n"
    "smua.trigger.endpulse.action = smua.SOURCE_IDLE\n"
    "smua.trigger.endpulse.stimulus = pulse_timer.EVENT_ID\n"
    "smua.trigger.count = 1\n"
    "smua.trigger.arm.count = 10\n"
    "smua.source.output = smua.OUTPUT_ON\n"
    "smua.trigger.initiate()\n"
    "waitcomplete()\n"
    # "nvbuffer1.clear()"
    # "nvbuffer2.clear()"
]
# Commands to set up Channel B which is measuring voltage
# smu_setup_b = [
#     "smub.source.func=smub.OUTPUT_DCAMPS",
#     "smub.source.leveli=0",  # Voltmeter=source 0 current
#     "smub.source.limitv={:.5e}".format(limit_V),
#     "smub.measure.nplc=1",
#     "smub.source.output = smub.OUTPUT_ON",
# ]

# Function to make a measurement with SMU
# Note without trying two sourcing currents, or setting some
# tiny smub current (1/1e6 of smua?) and reading back if it flows,
# you can't tell if the smub leads are disconnected,
# this could be a good idea
def get_data():
    # outer_I = float(smu.query("print(smua.measure.i())"))
    # outer_V = float(smu.query("print(smua.measure.v())"))
    # inner_V = float(smu.query("print(smub.measure.v())"))
    # return {"outer_I": outer_I, "outer_V": outer_V, "inner_V": inner_V}
    return outer_I


# Make the measurements. In a try->finally block to attempt to turn off SMU
# if measurement is interrupted
try:
    smu_write_list(smu_setup_a)
    # smu_write_list(smu_setup_b)
    have_errors = smu.query("print(errorqueue.count)").rstrip()
    print(have_errors)
    if not have_errors.isnumeric():
        have_errors = 0
    else:
        have_errors = round(float(have_errors))
    if have_errors > 0:
        print("SMU pulse init Errors:")
        for error in range(have_errors):
            smu.write("code, message = errorqueue.next()")
            print(smu.query("print(code,message)"))
    else:
        print("No pulse init errors")
    sleep(1)
    while True:
        print(smu.query("print(nvbuffer1.readings)"))
        # print(get_data())
        sleep(2)

# Stop on CTRL-C
except KeyboardInterrupt:
    pass

# When program stops, try and turn off SMU outputs and print errors
finally:

    # Turn off SMU outputs
    smu.write("smua.source.output = smua.OUTPUT_OFF")
    smu.write("smub.source.output = smub.OUTPUT_OFF")

    # Print any SMU errors
    have_errors = round(float(smu.query("print(errorqueue.count)").rstrip()))
    if have_errors > 0:
        print("SMU Errors:")
        for error in range(have_errors):
            smu.write("code, message = errorqueue.next()")
            print(smu.query("print(code,message)"))
    else:
        print("No errors")
