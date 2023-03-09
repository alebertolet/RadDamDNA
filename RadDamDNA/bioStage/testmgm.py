#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 3/5/23 5:00 PM

@author: alejandrobertolet
"""
import numpy as np
from mgm import mgm
from running import Simulator

# Data file, format is 'microdose' (3 columns with energy, specific energy and lineal energy)
data_file = '/Users/ai925/source/MGM/scripts/xray250keV.phsp'

# Other supported formats are:
# - A list of lineal energy values; each value in a new row
# - A list of bins with number of counts in each bin; each bin in a new row, separated by a space

# Initialize calculator.

calc = mgm.MicrodosimetryGammaCalculator(data_file, format='microdose', subsample=1000)
calc.CalculateDamage()
#calc.PlotComplexityDistribution(density=True)

# Gets the number of sites with DSBs
print(calc.getNumberOfSitesWithDSB(perTrack=True))

# Distribute damages over a nucleus with radius 3.5 um. Dose scales the number of damage sites.
# This returns a list of damage sites with position (x, y, z), and complexity
dose = 4 #Gy
damages = calc.DistributeDamageOverNucleus(dose=dose, radius=3.5, inTracks=True)
print('Number of sites for ' + str(dose) + ' Gy: ' + str(len(damages)))
#calc.PlotDistributedDamageOverNucleus()

# Set time options for simulation
initTime = 0
finalTime = 24 * 3600
nSteps = 24
timeOptions = [initTime, finalTime, nSteps]

# Parameters for dsb repair as a function of complexity. These are expected times for repair/misrepair
mu = np.linspace(1, 36, 19) * 3600 # Means 1 h for complexity 2, 36 hours for complexity 20
sigma = np.linspace(0.05, 0.35, 19) * 3600 # Sigma of the lognormal dist as a function of complexity
dsbpars = {'mu': mu, 'sigma': sigma}

# Another option is the 'standard' repair model with complex and non-complex dsb repair
# dsbpars = {'rComplex': 1.0e-5, 'rNCNC': 2.0e-4, '}

# Probability for a misrepaired dsb to become a lethal chromosome aberration
p_lethal_aberration = 0.08

# Probability for a cell to undergo apoptosis after non-repair damage happening at G1/S or G2/M checkpoints
p_apoptosis = 0.05
time_to_G1S_checkpoint = 10 * 3600 # hours to seconds
time_to_G2M_checkpoint = 24 * 3600 # hours to seconds

# Initialize Simulator
sim = Simulator(timeOptions=timeOptions, diffusionmodel='free', dsbmodel='lognormtimes', ssbmodel='none', bdmodel='none', nucleusMaxRadius = 3.5,
                 irradiationTime=0, doseratefunction=None, doseratefunctionargs=None, diffusionparams=None, dsbparams=dsbpars, ssbparams=None, bdparams=None,
                p_lethal_aberration=p_lethal_aberration, p_apoptosis=p_apoptosis, time_to_G1S_checkpoint=time_to_G1S_checkpoint, time_to_G2M_checkpoint=time_to_G2M_checkpoint)
# Load damage from mgm
sim.LoadDamageFromMGM(damages)
# Run simulation
nCellsSimulated = 50
sim.Run(nCellsSimulated, rereadDamageForNewRuns=False, basepath=None, maxDose=-1, version=None, plot=True, outputnorm=True, verbose=2)
survivalfraction = sim.GetSurvivalFraction()

print(survivalfraction)
