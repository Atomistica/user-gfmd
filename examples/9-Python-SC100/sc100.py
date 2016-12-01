#! /usr/bin/env python3
#
# This is an example of how to access the stiffness kernels from Python.

import numpy as np
import matplotlib.pyplot as plt

import gfmd

###

E = 5/2 # Young's modulus
nu = 1/2 # Poisson number

###

stiffness_1 = gfmd.SurfaceStiffness(1, 1, "sc100", "height 1")
stiffness_10 = gfmd.SurfaceStiffness(1, 1, "sc100", "height 10")
stiffness_100 = gfmd.SurfaceStiffness(1, 1, "sc100", "height 100")
stiffness_infty = gfmd.SurfaceStiffness(1, 1, "sc100", "height 100000")

special_points = [ ('$\Gamma$', (0., 0.)), ('$M$', (0.5, 0.5)), ('$X$', (0.5, 0.)), ('$\Gamma$', (0., 0.)) ]
colors = ['r', 'g', 'b', 'k']
labels = ['$h=a_0$', '$10 a_0$', '$100 a_0$', '$\infty$']

npts = 100
z = 0.0
for (sym1, (x1, y1)), (sym2, (x2, y2)) in zip(special_points[:-1], special_points[1:]):
    print(sym1, sym2)
    x = 2*np.pi*np.linspace(x1, x2, npts)
    y = 2*np.pi*np.linspace(y1, y2, npts)
    for color, label, stiffness in zip(colors, labels,
                                       [stiffness_1, stiffness_10, stiffness_100, stiffness_infty]):
        phi = stiffness(x, y)
        phizz = np.real(phi[:, 2, 2])
        plt.plot(np.linspace(z, z+1.0, npts), phizz/E, color+'-', label=label)
    z += 1.0
    labels = [None]*len(colors)

plt.xticks(np.arange(len(special_points)), [name for name, pts in special_points])
plt.gca().xaxis.grid(True)
plt.ylabel('Stiffness, $a_0 \Phi_{zz}/E$')
plt.xlabel('Wavevector, $a_0 q$')
plt.legend(loc='best')
plt.tight_layout()
plt.savefig('sc100.pdf')
plt.show()
