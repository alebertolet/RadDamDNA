#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/31/22 9:37 AM

@author: alejandrobertolet
"""

from RadDamDNA.medras import *
from RadDamDNA.damage import *
import random

maxdose = 2.0
base = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/'
particles = ['proton']
energies = ['10MeV']
for ip, e in enumerate(energies):
    basepath = base + particles[ip] + '/sims/' + e + '.txt/'
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
    damage = DamageToDNA()
    # Getting new random sequence along the available directories each time
    neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
    for j in neworder:
        path = basepath + str(j) + '/'
        damage.readSDDAndDose(path, defectiveChromosomeNumber=True)
    damage.populateDamages(stopAtDose=maxdose)
    damage.computeStrandBreaks()
    damage.printDamageCount()
    damage.writeSDD('/Users/ai925/source/workspace/repair/proton_10MeV_' + str(maxdose) + 'Gy_0.sdd')

    repair = MedrasRepair(damage=damage)
    repair.setVideoForOneEvent(recalculateDamages=False)
    #repair.repairSimulation(recalculateDamages=False)