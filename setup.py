# ======================================================================
# Atomistica - Interatomic potential library and molecular dynamics code
# https://github.com/Atomistica/atomistica
#
# Copyright (2005-2015) Lars Pastewka <lars.pastewka@kit.edu> and others
# See the AUTHORS file in the top-level Atomistica directory.
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

from __future__ import absolute_import, print_function

import glob
import os
import sys
import re

import numpy as np
from numpy.distutils.core import setup, Extension

import versioneer

###

inc_dirs = [ np.get_include(),
             'src/LAMMPS_STUB',
             'src/main',
             'src/mathutils',
             'src/stiffness_kernels',
             'src/surfaces'
             ]

###

setup(
    name = 'gfmd',
    version = versioneer.get_version(),
    cmdclass = versioneer.get_cmdclass(),
    description = 'GFMD - Green\' function molecular dynamics',
    maintainer = 'Lars Pastewka',
    maintainer_email = 'lars.pastewka@kit.edu',
    url = 'https://github.com/Atomistica/user-gfmd',
    download_url = 'https://github.com/Atomistica/user-gfmd/tarball/'+versioneer.get_version(),
    license = 'GPLv3',
    packages = [
        'gfmd'
        ],
    package_dir = {
        'gfmd': 'src/python/gfmd'
        },
    ext_modules = [
        Extension(
            '_gfmd',
            ['src/LAMMPS_STUB/memory.cpp',
             'src/main/crystal_surface.cpp',
             'src/main/surface_stiffness.cpp',
             'src/mathutils/table2d.cpp',
             'src/mathutils/linearalgebra.cpp',
             'src/stiffness_kernels/debug_stiffness.cpp',
             #'src/stiffness_kernels/ft_stiffness.cpp',
             'src/stiffness_kernels/isotropic_stiffness.cpp',
             'src/stiffness_kernels/sc100_stiffness.cpp',
             'src/surfaces/dia100_surface.cpp',
             'src/surfaces/dia111_surface.cpp',
             'src/surfaces/fcc100_surface.cpp',
             'src/surfaces/fcc111_surface.cpp',
             'src/surfaces/sc100_surface.cpp',
             'src/python/gfmdmodule.cpp',
             'src/python/py_bicubic.cpp',
             'src/python/py_surface_stiffness.cpp'
            ],
            include_dirs = inc_dirs
            )
        ],
    )
