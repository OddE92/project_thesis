#########################################################################
#       This script plots the linear approximation of the               #
#       particle trajectory generated by linear_approx.exe              #
#                                                                       #
#       Written by Odd-Einar C. Nervik                                  #
#########################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

trajectory = np.loadtxt('cpp_runge-kutta/RK4_approx_stepsizectrl.dat')                               # Reads the data from the RK4 approx
#trajectory2 = np.loadtxt('cpp_runge-kutta/RK4_approx_t.dat')


fig1 = plt.figure(1)
ax = fig1.gca(projection='3d')

ax.plot3D(trajectory[:,0], trajectory[:,1], trajectory[:,2])

ax.set_title('Trajectory RK4 step size control')
ax.set_xlabel('x [pc]')
ax.set_ylabel('y [pc]')
ax.set_zlabel('z [pc]')

#fig2 = plt.figure(2)

#plt.plot(trajectory[:,0], trajectory[:,1])

#fig2 = plt.figure(2)
#ax = fig2.gca(projection='3d')

#ax.plot3D(trajectory2[:,0], trajectory2[:,1], trajectory2[:,2])

#ax.set_title('Trajectory RK4')
#ax.set_xlabel('x')
#ax.set_ylabel('y')
#ax.set_zlabel('z')

plt.show()