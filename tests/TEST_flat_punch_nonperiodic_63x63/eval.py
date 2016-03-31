#! /usr/bin/env python

import glob
from math import pi, sqrt
import sys

import numpy as np

###

R = 100.0
G = 1.0
nu = 0.3
E = 2*G*(1+nu)/(1-nu*nu)

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

x = np.arange(nx)-float(nx)/2+0.5
#x = np.where(x > nx/2, x-nx, x)
y = np.arange(ny)-float(ny)/2+0.5
#y = np.where(y > ny/2, y-ny, y)

r_xy = np.sqrt( (x**2).reshape(-1,1) + (y**2).reshape(1,-1) )

### Fit contact radius and pressure

a = 13.7
N = np.min(f_xy[r_xy<a/2])

### Compute residual

fa_xy = np.where(r_xy<a, N/np.sqrt(1.0-(r_xy/a)**2), np.zeros_like(r_xy))
res = np.sum( (f_xy - fa_xy)**2 )

###

if len(sys.argv) == 2 and sys.argv[1] == '--dump':
    print 'Residual: ', res
    print 'Dumping f.out...'
    np.savetxt('f.out', np.transpose([r_xy.reshape(-1), f_xy.reshape(-1), fa_xy.reshape(-1)]))

if res > 1e-2:
    raise RuntimeError('Residual outside bounds: res = %f' % res)
