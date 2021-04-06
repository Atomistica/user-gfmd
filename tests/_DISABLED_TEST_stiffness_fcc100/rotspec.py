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
