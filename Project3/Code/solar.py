#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt

print "Enter the number of grid points"
n=raw_input()
print "Enter the name of the Figure"
Fig=raw_input()
#plt.plotfile('Solar_VV.dat', cols=(1,2), label='Mercury', delimiter=' ')
#plt.plotfile('Solar_VV.dat', cols=(5,6),label='Venus', newfig=False, delimiter=' ')
#plt.plotfile('Solar_VV.dat', cols=(9,10),label='Earth', newfig=False, delimiter=' ')
#plt.plotfile('Solar_VV.dat', cols=(13,14),label='Mars', newfig=False, delimiter=' ')
plt.plotfile('Solar_VV.dat', cols=(17,18),label='Jupiter', newfig=False, delimiter=' ')
plt.plotfile('Solar_VV.dat', cols=(21,22),label='Saturn', newfig=False, delimiter=' ')
plt.plotfile('Solar_VV.dat', cols=(25,26),label='Uranus', newfig=False, delimiter=' ')
plt.plotfile('Solar_VV.dat', cols=(29,30),label='Neptune', newfig=False, delimiter=' ')
plt.plotfile('Solar_VV.dat', cols=(33,34),label='Pluto', newfig=False, delimiter=' ')
plt.title('y vs x  '+n)
plt.xlabel('x(au)')
plt.ylabel('y(au)')
plt.legend()
plt.grid()
plt.savefig(Fig)
plt.show()
