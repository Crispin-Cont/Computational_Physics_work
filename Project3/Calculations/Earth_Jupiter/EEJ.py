#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt

with open(sys.argv[1],'r') as File, open(sys.argv[2],'r') as File2:
     print "Enter the number of grid points"
     n=raw_input()
     print "Enter the name of the Figure"
     Fig=raw_input()
     plt.plotfile(File, cols=(0,4), label='RK4', delimiter=' ')
     plt.plotfile(File2, cols=(0,4),label='Verlet', newfig=False, delimiter=' ')
     plt.title('Total Energy vs Time  '+n)
     plt.xlabel('Time(yr)')
     plt.ylabel(r'$E_{tot} (Au/yr)^2$')
     plt.legend()
     plt.grid()
     plt.savefig(Fig)
     plt.show()
