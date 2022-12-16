#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 12/16/22 4:40 PM

@author: alejandrobertolet
"""
from RadDamDNA.damage import *
import numpy as np
import random

fig = plt.figure()
fig.set_size_inches((12, 12))
ax = fig.add_subplot(111)

base = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/alpha/sims/'
nboot = 20
maxdose = 1.0

energies = ['8MeV', '6MeV', '4MeV', '2MeV', '1MeV']
energies.reverse()

E = np.array([1, 2, 4, 6, 8])
N = np.array(E.shape)
Nw = np.array(E.shape)
stdN = np.array(E.shape)
stdNw = np.array(E.shape)

for ie, e in enumerate(energies):
    basepath = base + e + '.txt/'
    nfiles = len(os.listdir(basepath))
    #Section to get what directories actually contains both dose and SDD. Disregard others!
    listOfAvailableDirs = []
    for j in range(nfiles):
        newpath = basepath + str(j) + '/'
        if str(j) in os.listdir(basepath):
            files = os.listdir(newpath)
            if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
                if os.path.getsize(newpath + 'DNADamage_sdd.txt') > 0: #only those with actual data
                    listOfAvailableDirs.append(j)

    meann = np.array([])
    meanwn = np.array([])
    # Bootstrapping
    for i in range(nboot):
        damage = DamageToDNA()
        # Getting new random sequence along the available directories each time
        neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
        for j in neworder:
            path = basepath + str(j) + '/'
            damage.readSDDAndDose(path, version='1.0')

        damage.populateDamages(stopAtDose=maxdose)
        damage.computeStrandBreaks()
        damage.getDamageClusterDistribution(plot=False)
        meann = np.append(meann, damage.meanN)
        meanwn = np.append(meanwn, damage.meanWeightedN)

    N[ie] = np.mean(meann)
    Nw[ie] = np.mean(meanwn)
    stdN[ie] = np.std(meann)
    stdNw[ie] = np.std(meanwn)



ax.plot(E, N, '-', label='Mean N', color='red')
ax.fill_between(E, N-stdN, N+stdN, color='red', alpha=0.15)
ax.plot(E, Nw, '-', label='Weighted mean N', color='blue')
ax.fill_between(E, N-stdNw, N+stdNw, color='blue', alpha=0.15)

ax.set_xlabel('Energy (MeV)')
ax.set_ylabel('Number of damages per site')
fig.tight_layout()
plt.show()