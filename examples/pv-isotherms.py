#! /usr/bin/env python3

"""Plots van der Waals isotherms using the Maxwell construction."""

import argparse

import numpy
from matplotlib import pyplot

import vdw


# Parse the args
parser = argparse.ArgumentParser()
parser.add_argument('-o', dest='path')
path = parser.parse_args().path

## Calculations ##

# Reduced variables
vr = numpy.linspace(0.36,20,1000)
Tr = numpy.array([0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4])
pr = [vdw.maxwell(T, vr, vdw.pr(T, vr)) for T in Tr]

# Mixed-phase boundary
vrl, vrg, prsat = [], [], []
for Trsat in numpy.linspace(0.6, 1, 1000):
    if Trsat < 1:
        ds = vdw.ds(Tr=Trsat)
        prsat.append(vdw.prsat(ds))
        vrsat = vdw.vrsat(ds)
        vrl.append(vrsat[0])
        vrg.append(vrsat[1])

## Plotting ##

fig = pyplot.figure(figsize=[4, 2.8])
fig.set_tight_layout(True)

for Tr_, pr_ in zip(Tr, pr):
    linewidth = 2 if Tr_==1 else 1
    pyplot.loglog(vr, pr_, 'k-', linewidth=linewidth, label='%.2f'%Tr_)

pyplot.plot([1],[1],'ko')

pyplot.loglog(vrl, prsat, 'k--')
pyplot.loglog(vrg, prsat, 'k--')

pyplot.xlim(0.3, 20)
pyplot.ylim(0.1, 10)

pyplot.xlabel('Reduced volume')
pyplot.ylabel('Reduced pressure')

ax = pyplot.gca()
ax.set_xticks([1,10])
ax.set_xticklabels(['1','10'])
ax.set_yticks([0.1,1,10])
ax.set_yticklabels(['0.1','1','10'])

if path:
    pyplot.savefig(path, transparent=True)
else:
    pyplot.show()
