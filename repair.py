#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/31/22 9:37 AM

@author: alejandrobertolet
"""

from RadDamDNA.medras import *

path = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/alpha/sims/1MeV.txt/0/'
repair = MedrasRepair()
repair.repairSimulation(path)