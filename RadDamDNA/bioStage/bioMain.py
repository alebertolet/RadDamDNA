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
diffusionparams = {'D': 2.0e-6, 'Dunits': 'um^2/s'}
dsbModel = 'standard'
dsbparams = {'NEHJ': True, 'rNCNC': 5.833e-4, 'rNCNCunits': 'rep/s', 'rComplex': 7.222e-5, 'rComplexunits': 'rep/s',
             'rMMEJ': 2.361e-6, 'rMMEJunits': 'rep/s', 'sigma': 0.25, 'sigmaUnits': 'um'}
ssbModel = 'standard'
ssbparams = {'rNC': 5.833e-4, 'rNCunits': 'rep/s', 'rC': 7.222e-5, 'rCunits': 'rep/s'}
bdModel = 'standard'
bdparams = {'r': 5.833e-4, 'runits': 'rep/s'}
nRuns = 2
irradiationTime = 25*3600
#doseratefunction = 'exponential'
doseratefunction = None
doserate = 0.5/3600
halflife = 1*3600


# Damage to be incorporated
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/'
maxDose = 0.5
version = '1.0'

sim = Simulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, irradiationTime=irradiationTime, doseratefunction=doseratefunction, doseratefunctionargs=[doserate, halflife],
                diffusionparams=diffusionparams, dsbparams=dsbparams, ssbparams=ssbparams, bdparams=bdparams)
sim.ReadDamage(basepath, maxDose, version)
sim.Run(nRuns, rereadDamageForNewRuns=False, basepath=basepath, maxDose=maxDose, version=version)