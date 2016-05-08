#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt

with open(sys.argv[1],'r') as File, open(sys.argv[2],'r') as File2:
     print "Enter the number of grid points"
     n=raw_input()
     print "Enter the name of the Figure"
     Fig=raw_input()
     plt.plotfile(File, cols=(0,3), label='RK4', delimiter=' ')
     plt.plotfile(File2, cols=(0,3),label='Verlet', newfig=False, delimiter=' ')
     plt.title('Verlet vs RK4  '+n)
     plt.xlabel('Time(years)')
     plt.ylabel(r'$E_{tot}$')
     plt.legend()
     plt.grid()
     plt.savefig(Fig)
     plt.show()
