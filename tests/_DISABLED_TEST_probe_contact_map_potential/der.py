# ======================================================================
# USER-GFMD - Elastic half-space methods for LAMMPS
# https://github.com/Atomistica/user-gfmd
#
# Copyright (2011-2016,2021)
#    Lars Pastewka <lars.pastewka@imtek.uni-freiburg>,
#    Tristan A. Sharp and others.
# See the AUTHORS file in the top-level USER-GFMD directory.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ======================================================================
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
