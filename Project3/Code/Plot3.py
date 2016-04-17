#!/usr/bin/env python
import os
import sys
import matplotlib.pyplot as plt

File=str(sys.argv[1])

plt.plotfile('File', cols=(0,1), label='test', delimiter=' ')
plt.title('TEst ')
plt.xlabel('Test ')
plt.ylabel('Test')
plt.legend()
plt.grid()
plt.show()
