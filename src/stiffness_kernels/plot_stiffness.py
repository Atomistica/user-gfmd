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
from matplotlib.pyplot import figure, setp, rc, show

rc("font", family="Arial", size=12)
rc("xtick", labelsize=8)
rc("ytick", labelsize=8)

cm       = 0.39370  # 1/inch
ylabel_x = -0.13

###

fig = figure(figsize=[7*cm,7*cm])

###

d, x, y, phixxr, phixxi, phiyyr, phiyyi, phizzr, phizzi, phixyr, phixyi, phiyzr, phiyzi, phixzr, phixzi = np.loadtxt('stiffness.out', unpack=True)

specsym, specd, specx, specy = np.loadtxt('specpoints.out', unpack=True,
                                          dtype=str)
specd = np.array(specd, dtype=float)

print specsym
print specd

###

p1 = fig.add_subplot(111, xticks=specd, xticklabels=specsym)
p1.set_ylabel(r'$\Phi$')
p1.set_xlabel(r'$q$')

p1.plot(d,    phixxr, label=r'$\Phi_{xx}$')
p1.plot(d,    phiyyr, label=r'$\Phi_{yy}$')
p1.plot(d,    phizzr, label=r'$\Phi_{zz}$')
p1.plot(d, 10*phixyr, label=r'$10 \Phi_{xy}$')
p1.plot(d,  4*phiyzi, label=r'$-4 i \Phi_{yz}$')
p1.plot(d, -4*phixzi, label=r'$4 i \Phi_{xz}$')

p1.set_xlim(specd[0], specd[-1])
p1.set_ylim(-1.5, 5)
p1.get_xaxis().grid(True)

###

p1.legend(frameon=False, prop={'size':8})

###

fig.subplots_adjust(left=0.15, right=0.98, top=0.98, bottom=0.12, hspace=0.15)

fig.savefig('plot_stiffness.pdf')
