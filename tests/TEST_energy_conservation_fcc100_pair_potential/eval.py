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

import sys

import numpy as np

sys.path += [ '..' ]
import lammps_log

###

l = lammps_log.log('log.lammps')
e = l.get('TotEng')

de = np.max(np.abs(e-np.mean(e)))

if de/np.mean(e) > 1e-4:
    raise RuntimeError('Fluctuations in total energy exceed threshold')

