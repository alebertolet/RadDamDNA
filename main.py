#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from RadDamDNA.SDD import *
import random

damage = DamageToDNA()
basepath = '/Users/alebertolet/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/proton/sims/1MeV.txt/'
#neworder = random.sample(range(250), 2)
neworder = [0, 1]
for i in neworder:
    path = basepath + str(i) + '/'
    damage.readSDDAndDose(path)

damage.populateDamages(getVideo=False)
damage.computeStrandBreaks()
damage.printDamageCount()
damage.plotDoseResponseCurve()
#damage.produce3DImage()
#damage.produce2DImages()

sddfile = 'test.sdd'
damage.writeSDD(sddfile)
