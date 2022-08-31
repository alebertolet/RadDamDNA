#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from RadDamDNA.SDD import *
import random

damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/alpha/sims/1MeV.txt/'
neworder = random.sample(range(225), 150)
neworder = [0, 1]
for i in neworder:
    path = basepath + str(i) + '/'
    damage.readSDDAndDose(path, defectiveChromosomeNumber=True)

damage.populateDamages(getVideo=False)
damage.computeStrandBreaks()
damage.printDamageCount()
damage.getDoseResponseCurve('SSB')
damage.getDoseResponseCurve('DSB')
damage.getDoseResponseCurve('BD')
#damage.produce3DImage()
#damage.produce2DImages()

sddfile = 'test.sdd'
damage.writeSDD(sddfile)
