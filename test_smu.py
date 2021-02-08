"""
FILE: test_smu.py
AUTHOR: R Ellingham
DATE MODIFIED: Feb 2021
"""

from k2600 import K2600
import k2600

import time

k = K2600('TCPIP0::169.254.0.1::inst0::INSTR')
k2600.debug_enable() # run verbose version program

def main(test):
    test = int(test)
    ## EXAMPLE 2-WIRE MEASUREMENT ##
    if test == 0:
        k.display.smua.measure.func(k,k.display.MEASURE_OHMS) # display ohms measurement

        k.smua.sense(k,k.smua.SENSE_LOCAL) # choose resistance measurement config`

        k.smua.source.leveli(k,10e-6) # set current source value

        k.smua.source.limitv(k,20) # set voltage limit (if too low this may restrict the current sourced)

        k.smua.source.output(k,1) # turn on channel a

        for sample in range(200):
            # take 1000 measurements of resistance
            print(k.smua.measure.r(k))

        k.disconnect() # disconnects smus and shuts off both a&b output channels

    if test == 1:
        ## EXAMPLE 4-WIRE MEASUREMENT ##

        k.display.screen(k,k.display.SMUA_SMUB) # display both smu a and b
        # k.display.smua.measure.func(k,k.display.MEASURE_DCVOLTS) # display volts
        # k.display.smub.measure.func(k,k.display.MEASURE_DCVOLTS) # display volts

        k.smua.source.func(k,k.smua.OUTPUT_DCAMPS)
        k.smub.source.func(k,k.smub.OUTPUT_DCAMPS)

        k.smua.source.leveli(k,10e-6) # set current source values
        k.smub.source.leveli(k,0)

        k.smua.source.limitv(k,10) # set voltage limit (if too low this may restrict the current sourced)
        k.smub.source.limitv(k,10)

        k.smua.measure.nplc(k,1) # set number of powerline cycle to integrate over for volts measurement to 1
        k.smub.measure.nplc(k,1)

        k.smua.source.output(k,1) # turn on channel b
        k.smub.source.output(k,1) # turn on channel a

        for sample in range(10):
            # take x measurements of resistance
            ia = k.smua.measure.i(k)
            va = k.smua.measure.v(k)
            vb = k.smub.measure.v(k)
            print(ia,"A(a)\t",end='')
            print(va,"V(a)\t",end='')
            print(vb,"V(b)\t",end='')
            print(float(va)/float(ia),"ohms(a)\t",end='')
            print(float(vb)/float(ia),"ohms(b)")
            time.sleep(0.5)

        k.disconnect() # disconnects smus and shuts off both a&b output channels

    return 0

if __name__ == "__main__":
	test_num = 0 # Default input parameter
	import sys
	if len(sys.argv)>1: test_num = (sys.argv[1])
	main(test_num)
