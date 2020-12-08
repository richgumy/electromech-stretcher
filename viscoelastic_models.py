## import libraries
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize

# Step function for general use (only approximate,no 0- and 0+)
step_input = lambda x,steps : np.append(np.zeros(1),np.linspace(x,x,steps-1))
# Stress relaxation function for Maxwell model
maxwell_stress_lax = lambda e0,E,mu,t : E*e0*np.exp(-(E/mu)*t)
# Creep function
maxwell_strain_creep = lambda s0,E,mu,t : s0*((1/E)*t+1/mu)
# Dirac delta function
# dt = 2a -> where a is determines the width of the dirac delta function
dirac = lambda t,dt : 2/(dt*np.sqrt(np.pi))*np.exp(-(2*t/dt)**2)
# Relaxation function
kv_stress_lax = lambda e0,E,mu,t : e0*(E + mu*dirac(t,t[1]))
# Creep function
kv_strain_creep = lambda s0,E,mu,t : (s0/E)*(1-np.exp(-(E/mu)*t))

# Creep function
kv_strain_super_creep = lambda s0,E,mu,tau,t : (s0/E)*(1-np.exp(-(E/mu)*(tau-t))) # change this to match above equation
# Abitrary parameter values
s0 = [0.1,0.2,0.3]
tau = [20,40,60]
E = 1
mu = 2
steps = 100
t = np.linspace(0,10,steps)
for i in range(len(s0)):
    plt.plot(kv_strain_super_creep(s0[i],E,mu,tau[i],t),'b-')
    plt.plot(step_input(s0[i],steps),'r-')
plt.xlabel('time')
plt.ylabel('r=Stress, b=Strain')
