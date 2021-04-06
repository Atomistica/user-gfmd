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

import glob
from math import pi, sqrt
import sys

import numpy as np

###

R = 100.0
#G = 1.0
#nu = 0.3
#E = 2*G*(1+nu)/(1-nu*nu)
E = 1.39 # Contact modulus

###

fns = glob.glob('gfmd.*.r.f2.out')
fns.remove('gfmd.0.r.f2.out')
if len(fns) > 1:
    raise RuntimeError('More than one GFMD output found. Not sure which one to use.')

f_xy = np.loadtxt(fns[0])
nx, ny = f_xy.shape

###

r0 = 3.0
rbins = [ r0 ]
r2 = r0
while r2 < nx/4:
    r2 = sqrt(r2*r2+r0*r0)
    rbins += [ r2 ]

###

x = np.arange(nx)+0.5
x = np.where(x > nx/2, x-nx, x)
y = np.arange(ny)+0.5
y = np.where(y > ny/2, y-ny, y)

r_xy = np.sqrt( (x**2).reshape(-1,1) + (y**2).reshape(1,-1) )

### Pressure as a function of distance

N = np.sum(f_xy)
a = R*(3./4*( N/(E*R**2) ))**(1./3)
p0 = 3*N/(2*pi*a*a)

### Compute residual

fa_xy = np.where(r_xy<a, p0*np.sqrt(1-(r_xy/a)**2), np.zeros_like(r_xy))
res = np.sum( (f_xy - fa_xy)**2 )

if res > 1e-2:
    raise RuntimeError('Residual outside bounds: res = %f' % res)

###

if len(sys.argv) == 2 and sys.argv[1] == '--dump':
    print('Residual: ', res)
    print('Dumping f.out...')
    np.savetxt('f.out', np.transpose([r_xy.reshape(-1), f_xy.reshape(-1), fa_xy.reshape(-1)]))
