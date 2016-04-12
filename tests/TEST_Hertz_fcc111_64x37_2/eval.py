# ======================================================================
# USER-GFMD - Green's function molecular dynamics for LAMMPS
# https://github.com/Atomistica/user-gfmd
#
# Copyright (2011-2016) Lars Pastewka <lars.pastewka@kit.edu>, Tristan A. Sharp
# and others. See the AUTHORS file in the top-level USER-GFMD directory.
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

R = 100
E = 1.54

###

# FCC111 unit cell contains two atoms. We load forces on both.
fns = glob.glob('gfmd.*.r.f2.out')
fns.remove('gfmd.0.r.f2.out')
if len(fns) > 1:
    raise RuntimeError('More than one GFMD output found. Not sure which one to use.')

f1_xy = np.loadtxt(fns[0])
nx, ny = f1_xy.shape

fns = glob.glob('gfmd.*.r.f5.out')
fns.remove('gfmd.0.r.f5.out')
if len(fns) > 1:
    raise RuntimeError('More than one GFMD output found. Not sure which one to use.')

f2_xy = np.loadtxt(fns[0])

###

dx = sqrt(3.0) # Cell size in y-direction
A0 = 1.0*dx/2 # two atoms per unit cell

###

sx = dx*nx
sy = ny

x1 = dx*np.arange(nx)
x1 = np.where(x1 > sx/2, x1-sx, x1)
y1 = np.arange(ny)
y1 = np.where(y1 > sy/2, y1-sy, y1)

r1_xy = np.sqrt( (x1**2).reshape(-1,1) + (y1**2).reshape(1,-1) )

# Second atom is shifted
x2 = dx*(np.arange(nx)+0.5)
x2 = np.where(x2 > sx/2, x2-sx, x2)
y2 = np.arange(ny)+0.5
y2 = np.where(y2 > sy/2, y2-sy, y2)
r2_xy = np.sqrt( (x2**2).reshape(-1,1) + (y2**2).reshape(1,-1) )

### Pressure as a function of distance

r_xy = np.append(r1_xy.reshape(-1), r2_xy.reshape(-1))
f_xy = np.append(f1_xy.reshape(-1), f2_xy.reshape(-1))

### Fit contact radius and pressure

N = np.sum(f_xy)
a = R*(3./4*( N/(E*R**2) ))**(1./3)
p0 = 3*N/(2*pi*a*a)

### Compute residual

# pressure
p_xy = f_xy/A0
pa_xy = np.where(r_xy<a, p0*np.sqrt(1-(r_xy/a)**2), np.zeros_like(r_xy))
res = np.sum( (p_xy - pa_xy)**2 )

###

if len(sys.argv) == 2 and sys.argv[1] == '--dump':
     print 'Residual: ', res
     print 'Dumping p.out...'
     np.savetxt('p.out', np.transpose([r_xy.reshape(-1), p_xy.reshape(-1), pa_xy.reshape(-1)]))

# Residual is high. Not sure if something is wrong with FCC 111. Residual on
# FCC 100 is much lower.
if res > 0.002:
    raise RuntimeError('Residual outside bounds: res = %f' % res)
