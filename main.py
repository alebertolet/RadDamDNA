#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from RadDamDNA.SDD import *
import random

damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/proton/sims/1MeV.txt/'
neworder = random.sample(range(250), 50)
for i in neworder:
    path = basepath + str(i) + '/'
    try:
        damage.readSDDAndDose(path)
    except:
        pass

damage.populateDamages()
damage.computeStrandBreaks()
damage.printDamageCount()
damage.plotDoseResponseCurve()
damage.produce3DImage()
damage.produce2DImages()