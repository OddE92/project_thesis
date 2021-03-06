#########################################################################
#       This script plots the linear approximation of the               #
#       particle trajectory generated by linear_approx.exe              #
#                                                                       #
#       Written by OddE92                                               #
#########################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

trajectory = np.loadtxt('cpp_trajectory_linear_approximation/linear_approx_traj.dat')               # Reads the data from the generated file.



fig = plt.figure(3)
ax = fig.gca(projection='3d')

ax.plot3D(trajectory[:,0], trajectory[:,1], trajectory[:,2])

plt.xlabel('Bx')
plt.ylabel('By')
ax.set_zlabel('Bz')

#plt.show()

fig2 = plt.figure(4)

plt.plot(trajectory[:,0], trajectory[:,1])

plt.show()
