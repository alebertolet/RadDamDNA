#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from RadDamDNA.damage import *
import random

maxdose = 2.0
damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/alpha/sims/8MeV.txt/'
nfiles = len(os.listdir(basepath))
# Section to get what directories actually contains both dose and SDD. Disregard others!
listOfAvailableDirs = []
for j in range(nfiles):
    newpath = basepath + str(j) + '/'
    if str(j) in os.listdir(basepath):
        files = os.listdir(newpath)
        if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
            if os.path.getsize(newpath + 'DNADamage_sdd.txt') > 0:  # only those with actual data
                listOfAvailableDirs.append(j)
neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
for i, e in enumerate(neworder):
    path = basepath + str(e) + '/'
    damage.readSDDAndDose(path, version='1.0')

damage.populateDamages(getVideo=False, stopAtDose=maxdose)
damage.computeStrandBreaks()
damage.printDamageCount()
#damage.getDoseResponseCurve(q='SSB')
#damage.getDoseResponseCurve(q='DSB')
damage.getDamageClusterDistribution(savefile='alpha8MeV2Gy.png')
#damage.getDoseResponseCurve(q='BD')
#damage.produce3DImage()
#damage.produce2DImages()

#sddfile = 'test.sdd'
#damage.writeSDD(sddfile)