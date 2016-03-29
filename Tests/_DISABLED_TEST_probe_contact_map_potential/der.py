#! /usr/bin/env python

import sys
import numpy as np

###

def dx(x, y):
    """ First derivative using finite differences on an arbitrary grid.
        O(dx**2)
    """

    x        = np.asarray(x)
    y        = np.asarray(y)

    dx1      = x[2:] - x[1:-1]
    dx2      = x[1:-1] - x[:-2]

    dx1_sq   = dx1*dx1
    dx2_sq   = dx2*dx2

    r        = y.copy()
    # O(dx) on the boundary
    r[0]     = ( y[1]  - y[0]  )/( x[1]  - x[0]  )
    r[-1]    = ( y[-1] - y[-2] )/( x[-1] - x[-2] )
    # O(dx**2) in the middle
    r[1:-1]  = ( ( y[2:] - y[1:-1] )/dx1_sq - ( y[0:-2] - y[1:-1] )/dx2_sq )/(1.0/dx1+1.0/dx2)

    return r

###

x, y = np.loadtxt('pot.out', usecols=[1,2], unpack=True)
dy = dx(x, y)

np.savetxt('pot_nd.out', np.transpose([x,y,-dy]))
