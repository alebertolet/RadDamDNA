#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/28/22 7:59 PM

@author: alejandrobertolet
"""
import numpy as np

from RadDamDNA.SDD import *
import random

from scipy import interpolate

nboot = 20

maxdose = 8.0
Dose = np.linspace(0, maxdose, 100)
DSB = np.zeros(Dose.shape)
alldsbs = np.zeros([len(Dose), nboot])
SSB = np.zeros(Dose.shape)
BD = np.zeros(Dose.shape)
ntries = np.zeros(Dose.shape)
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/proton/sims/1MeV.txt/'
for i in range(nboot):
    damage = DamageToDNA()
    neworder = random.sample(range(225), 225)
    for j in neworder:
        path = basepath + str(j) + '/'
        damage.readSDDAndDose(path)
    damage.populateDamages(stopAtDose=maxdose)
    damage.computeStrandBreaks()
    dose, dsb = damage.getDoseResponseCurve(plot=False, q='dsb')
    dose, ssb = damage.getDoseResponseCurve(plot=False, q='ssb')
    dose, bd = damage.getDoseResponseCurve(plot=False, q='bd')
    dose = np.array(dose)
    dsb = np.array(dsb)
    ssb = np.array(dsb)
    bd = np.array(dsb)
    dose = np.insert(dose, 0, 0)
    dsb = np.insert(dsb, 0, 0)
    ssb = np.insert(ssb, 0, 0)
    bd = np.insert(bd, 0, 0)
    fdsb = interpolate.interp1d(dose, dsb)
    fssb = interpolate.interp1d(dose, ssb)
    fbd = interpolate.interp1d(dose, bd)
    for id, d in enumerate(Dose):
        if d >= 0 and d < np.max(dose):
            DSB[id] += fdsb(d)
            alldsbs[id, i] = fdsb(d)
            SSB[id] += fssb(d)
            BD[id] += fbd(d)
            ntries[id] += 1

vardsb = np.zeros(Dose.shape)
for j in range(alldsbs.shape[0]):
    vardsb[j] = np.var(alldsbs[j,:])

stddsb = np.sqrt(vardsb)

DSB = DSB / ntries
SSB = DSB / ntries
BD = BD / ntries

fig = plt.figure()
fig.set_size_inches((4, 4))
ax = fig.add_subplot(111)
ax.plot(Dose, DSB, 'k-')
ax.fill_between(Dose, DSB-stddsb, DSB+stddsb)
ax.set_xlim(0, 1.01*maxdose)
ax.set_ylim(0, None)
ax.set_xlabel('Dose (Gy)')
ax.set_ylabel('DSB')
ax.grid()
fig.tight_layout()
plt.show()