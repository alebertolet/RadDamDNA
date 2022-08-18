#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 11:25 AM

@author: alejandrobertolet
"""

from SDD import *

damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/proton/sims/0.8MeV.txt/'
for i in range(0, 100):
    path = basepath + str(i) + '/'
    try:
        damage.readSDDAndDose(path)
    except:
        pass

damage.populateDamages()
damage.computeStrandBreaks()
damage.plotDoseResponseCurve()