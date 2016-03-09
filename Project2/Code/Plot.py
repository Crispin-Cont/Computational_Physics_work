#!/usr/bin/env python
import matplotlib.pyplot as plt

plt.plotfile('I_5.dat', cols=(0,1),label=r'$|u_0(\rho)|^2$', delimiter=' ')
plt.plotfile('I_5.dat', cols=(0,2),label=r'$|u_1(\rho)|^2$', newfig=False, delimiter=' ')
plt.plotfile('I_5.dat', cols=(0,3),label=r'$|u_2(\rho)|^2$',newfig=False, delimiter=' ')
plt.title(r'Relative Radial Wavefunction for $\omega_r$=0.5')
plt.xlabel(r'Radial Coordinate $ \rho$')
plt.ylabel(r'Radial Wavefunction: $|u(\rho)|^2$')
plt.legend()
plt.grid()
plt.show()
