#########################################################################
#       This script plots a heat map of the particle distribution       #
#       in the (x,y)-plane at recorded times.                           #
#                                                                       #
#       Written by OddE92                                               #
#########################################################################

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize

rank1 = np.loadtxt('data/r_rank1.dat')
rank2 = np.loadtxt('data/r_rank2.dat')
rank3 = np.loadtxt('data/r_rank3.dat')
rank4 = np.loadtxt('data/r_rank4.dat')
rank5 = np.loadtxt('data/r_rank5.dat')


print(len(rank1))
time1 = []
time2 = []
time3 = []

numPairs = 46
a = 18
b = 36
c = 45

for i in range(0, len(rank1), numPairs):
    time1.append( [ rank1[i+a, 0], rank1[i+a, 1] ] )
    #time1.append( [ rank2[i+a, 0], rank2[i+a, 1] ] )
    #time1.append( [ rank3[i+a, 0], rank3[i+a, 1] ] )
    #time1.append( [ rank4[i+a, 0], rank4[i+a, 1] ] )
    #time1.append( [ rank5[i+a, 0], rank5[i+a, 1] ] )
    time2.append( [ rank1[i+b, 0], rank1[i+b, 1] ] )
    #time2.append( [ rank2[i+b, 0], rank2[i+b, 1] ] )
    #time2.append( [ rank3[i+b, 0], rank3[i+b, 1] ] )
    #time2.append( [ rank4[i+b, 0], rank4[i+b, 1] ] )
    #time2.append( [ rank5[i+b, 0], rank5[i+b, 1] ] )
    time3.append( [ rank1[i+c, 0], rank1[i+c, 1] ] )
    #time3.append( [ rank2[i+c, 0], rank2[i+c, 1] ] )
    #time3.append( [ rank3[i+c, 0], rank3[i+c, 1] ] )
    #time3.append( [ rank4[i+c, 0], rank4[i+c, 1] ] )
    #time3.append( [ rank5[i+c, 0], rank5[i+c, 1] ] ) 

Pytime1 = np.array(time1)
Pytime2 = np.array(time2)
Pytime3 = np.array(time3)


f, (ax1, ax2, ax3) = plt.subplots(3)
f.suptitle('Diffusion after 10^2, 10^4 and 10^5 years')

#figure1 = plt.figure(1)
#ax1 = plt.gca()
ax1.set_xlim([-600, 600])
ax1.set_ylim([-600, 600])

#ax1.title("100 years")

ax1.scatter(Pytime1[:,1], Pytime1[:,0])

#figure2 = plt.figure(2)

#ax2 = plt.gca()
ax2.set_xlim([-600, 600])
ax2.set_ylim([-600, 600])

#ax2.title("1000 years")

ax2.scatter(Pytime2[:,1], Pytime2[:,0])



#figure3 = plt.figure(3)
#ax3 = plt.gca()
ax3.set_xlim([-600, 600])
ax3.set_ylim([-600, 600])

#ax3.title("10'000 years")

ax3.scatter(Pytime3[:,1], Pytime3[:,0])


plt.show()
