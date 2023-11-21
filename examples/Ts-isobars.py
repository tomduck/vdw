#! /usr/bin/env python3

import argparse

import numpy
from matplotlib import pyplot

import vdw

parser = argparse.ArgumentParser()
parser.add_argument('-o', dest='path')
path = parser.parse_args().path

## Calculations ##

# Reduced and non-dimensional variables
vr = numpy.linspace(0.36, 20, 1000)
pr = numpy.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2])
s = [vdw.s(vdw.Tr(p, vr), vr) for p in pr]
Tr = [vdw.Tr(p, vr) if p >= 1 else
      vdw.maxwell(vdw.Tr(p, vr), vr, p, line='isobar') for p in pr]

# Mixed-phase boundary
sl, sg, Trsat = [], [], []
for prsat in numpy.linspace(0.3, 1, 1000):
    if prsat < 1:
        ds = vdw.ds(pr=prsat)
        Trsat.append(vdw.Trsat(ds))
        vrl, vrg = vdw.vrsat(ds)
        sl.append(vdw.s(vdw.Tr(prsat, vrl), vrl))
        sg.append(vdw.s(vdw.Tr(prsat, vrg), vrg))

## Plotting ##

fig = pyplot.figure(figsize=[4, 2.8])
fig.set_tight_layout(True)

for s_, T, p in zip(s, Tr, pr):
    linewidth = 2 if p==1 else 1
    pyplot.plot(s_, T, 'k-', linewidth=linewidth)

pyplot.plot([0],[1],'ko')

pyplot.plot(sl, Trsat, 'k--')
pyplot.plot(sg, Trsat, 'k--')

pyplot.xlim(-2, 2)
pyplot.ylim(0.8, 1.2)

pyplot.xlabel('Nondimensional Entropy')
pyplot.ylabel('Reduced Temperature')

if path:
    pyplot.savefig(path, transparent=True)
else:
    pyplot.show()
