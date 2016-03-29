#! /usr/bin/env python

from math import sqrt

import numpy as np

f = open('spec_pts_fcc100x2.dat', 'r')
f.readline()
spec_pts = np.loadtxt(f)
f.close()

#R = np.array([ [ 1.0/sqrt(2.0), -1.0/sqrt(2.0) ],
#               [ 1.0/sqrt(2.0),  1.0/sqrt(2.0) ] ])

R = np.array([ [ 0.0, -1.0 ],
               [ 1.0,  0.0 ] ])


rot_spec_pts = [ ]
for p in spec_pts:
    rot_spec_pts += [ np.dot(R, p) ]

f = open('spec_pts_fcc100_90.dat', 'w')
f.write('%i\n' % len(rot_spec_pts))
np.savetxt(f, rot_spec_pts)
f.close()
