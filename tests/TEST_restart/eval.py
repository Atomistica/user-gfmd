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

import numpy as np

###

def forward_to(f, token):
    l = f.readline()
    if l:
        l = l.split()
    while l is not None and (len(l) < 1 or l[0] != token):
        l = f.readline()
        if l:
            l = l.split()

def read_col(f, col, token):
    x = []
    l = ['']
    while l is not None and l[0] != token:
        l = f.readline()
        if l:
            l = l.split()
        if l is not None and l[0] != token:
            x += [float(l[col])]   
    return x

f = open('log.lammps')
forward_to(f, 'Step')
en1 = read_col(f, 2, 'Loop')

forward_to(f, 'Step')
en2 = read_col(f, 2, 'Loop')

assert abs(en1[-1]-en2[0]) < 1e-6