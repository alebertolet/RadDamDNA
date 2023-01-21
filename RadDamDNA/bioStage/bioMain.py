#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 4:56 PM

@author: alejandrobertolet
"""
from RadDamDNA.bioStage.running import Simulator

# Time options is a list with initial, final times and number of steps (or a list of custom time points as 4th arg)
timeOptions = [0, 25*3600, 24]
nucleusMaxRadius = 4.65
diffusionModel = 'free'
dsbModel = 'standard'
ssbModel = 'standard'
bdModel = 'standard'
nRuns = 2
irradiationTime = 25*3600
doseratefunction = 'exponential'
doserate = 0.5/3600
halflife = 1*3600


# Damage to be incorporated
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/'
maxDose = 10000
version = '1.0'

sim = Simulator(timeOptions, diffusionModel, dsbModel, ssbModel, bdModel, nucleusMaxRadius,
                irradiationTime, doseratefunction, [doserate, halflife])
sim.ReadDamage(basepath, maxDose, version)
sim.Run(nRuns, rereadDamageForNewRuns=False, basepath=basepath, maxDose=maxDose, version=version)