#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 4:56 PM

@author: alejandrobertolet
"""
from RadDamDNA.bioStage.running import Simulator

# Time options is a list with initial, final times and number of steps (or a list of custom time points as 4th arg)
timeOptions = [0, 25*3600, 100]
nucleusMaxRadius = 4.65
diffusionModel = 'free'
nRuns = 10
# Damage to be incorporated
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/'
maxDose = 0.25
version = '1.0'

sim = Simulator(timeOptions, diffusionModel, nucleusMaxRadius)
sim.ReadDamage(basepath, maxDose, version)
sim.Run(nRuns, rereadDamageForNewRuns=True, basepath=basepath, maxDose=maxDose, version=version)