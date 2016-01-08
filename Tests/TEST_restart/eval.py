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