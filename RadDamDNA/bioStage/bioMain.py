#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 4:56 PM

@author: alejandrobertolet
"""
from RadDamDNA.damage import *
from RadDamDNA.bioStage.running import Simulator
import random

####################
#### SCRIPT TO READ DAMAGE DATA FROM TOPAS-NBIO
maxdose = 2.0
damage = DamageToDNA()
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/'
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
##############

# Time options is a list with initial, final times and number of steps (or a list of custom time points as 4th arg)
timeOptions = [0, 25*3600, 100]
nucleusMaxRadius = 4.65
diffusionModel = 'free'
sim = Simulator(damage, timeOptions, diffusionModel, nucleusMaxRadius)