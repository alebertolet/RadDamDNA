#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from RadDamDNA.damage import *
import random

damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/alpha/sims/1MeV.txt/'
neworder = random.sample(range(225), 150)
neworder = [0, 1]
times = [0, 0]
for i, e in enumerate(neworder):
    path = basepath + str(e) + '/'
    damage.readSDDAndDose(path, defectiveChromosomeNumber=True, particleTime=times[i])

damage.populateDamages(getVideo=False)
damage.computeStrandBreaks()
damage.printDamageCount()
#damage.getDoseResponseCurve(q='SSB')
#damage.getDoseResponseCurve(q='DSB')
damage.getDoseResponseCurve(q='BD')
#damage.produce3DImage()
#damage.produce2DImages()

sddfile = 'test.sdd'
damage.writeSDD(sddfile)