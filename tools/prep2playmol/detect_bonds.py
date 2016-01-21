#!/usr/bin/python
import sys
import string
import math

natoms = int(sys.stdin.readline())
title = sys.stdin.readline()
atom = []

x = []
y = []
z = []
for i in range(natoms):
    line = string.split(sys.stdin.readline())
    atom.append(line[0])
    x.append(float(line[1]))
    y.append(float(line[2]))
    z.append(float(line[3]))

cutoff = 1.7
for i in range(natoms-1):
    for j in range(i+1,natoms):
        r = math.sqrt(math.pow(x[i] - x[j],2) + math.pow(y[i] - y[j],2) + math.pow(z[i] - z[j],2))
        if (r <= cutoff): print atom[i], atom[j]

