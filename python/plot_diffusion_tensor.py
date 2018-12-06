#########################################################################
#       This script plots the diffusion tensor of the                   #
#       magnetic fields generated by generate_samples.exe               #
#                                                                       #
#       Written by Odd-Einar C. Nervik                                  #
#########################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

d_i = np.loadtxt('data/eigenvalues.dat') 

i = len(d_i)

t = []

for k in range(0, i):
    t.append(10**k)                     # t now holds the time-vector to plot the values against.


fig1 = plt.figure(1)

plt.title("Diffusion tensor")
plt.xlabel("Time in years")
plt.ylabel(r'$ \mathrm{d}_i$ $[\mathrm{pc}^2 / \mathrm{y}]$')

for k in range(0, 3):
    plt.loglog(t, d_i[:,k], linestyle='-', marker='x')

plt.show()