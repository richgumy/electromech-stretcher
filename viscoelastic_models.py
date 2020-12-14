## import libraries
import csv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib import cm # colour map
import numpy as np
from scipy import optimize

# Step function for general use (only approximate,no 0- and 0+)
# t is where the step starts in time (within the range of steps)
step_input = lambda a,steps,tau, : np.append(np.zeros(t),np.linspace(a,a,(steps-tau)))
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
kv_strain_super_creep = lambda s0,E,mu,t,tau : (s0/E)*(1-np.exp(-(E/mu)*(t-tau)))
# Abitrary parameter values
s0 = [0.1,0.2,0.3]
tau = [2,4,6]
E = 1
mu = 2
steps = 100
t = np.linspace(0,10,steps)
print(t)
for i in range(len(s0)):
    plt.plot(t,kv_strain_super_creep(s0[i],E,mu,t,tau[i]),'b-')
    plt.plot(t,step_input(s0[i],steps,tau[i]),'r-')
    # plt.plot(t,s0[i]*heaviside(t-tau[i], 0.5))
plt.xlabel('time')
plt.ylabel('r=Stress, b=Strain')

plt.plot()
plt.show()
