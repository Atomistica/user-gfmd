#! /usr/bin/env python

import sys

import numpy as np

sys.path += [ '../../../tools/python' ]
import pizza.log

###

l = pizza.log.log('log.lammps')
e = l.get('TotEng')

de = np.max(np.abs(e-np.mean(e)))

if de/np.mean(e) > 1e-4:
    raise RuntimeError('Fluctuations in total energy exceed threshold')

