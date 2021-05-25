# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 12:13:11 2021

@author: rel80
"""

import matplotlib.pyplot as plt

# 2_7-5_4Epin_20mm capacitance (@10kHz) strain data
ci_arr = [53,32,24,20]
co_arr = [32,25,22,19]
strain_arr = [0,10,20,30]

plt.plot(strain_arr,ci_arr,'r')
plt.plot(strain_arr,co_arr,'b')
plt.xlabel("Strain [%]")
plt.ylabel("Capacitance [pF]")
plt.legend(["Inner electrodes","Outer electrodes"])
plt.show()