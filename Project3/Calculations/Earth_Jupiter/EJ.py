#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt

print "Enter the number of grid points"
n=raw_input()
print "Enter the name of the Figure"
Fig=raw_input()
plt.plotfile('Ea_Ju_VVE3.dat', cols=(1,2), label='Earth', delimiter=' ')
plt.plotfile('Ea_Ju_VVE3.dat', cols=(5,6),label='Jupiter', newfig=False, delimiter=' ')
plt.title('y vs x  '+n)
plt.xlabel('x(au)')
plt.ylabel('y(au)')
plt.legend()
plt.axis([-10,10,-10,10])
plt.grid()
plt.savefig(Fig)
plt.show()
