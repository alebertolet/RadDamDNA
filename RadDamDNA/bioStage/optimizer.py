#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 1/22/23 3:32 PM

@author: alejandrobertolet
"""
import csv
import numpy as np
import matplotlib.pyplot as plt

from RadDamDNA.bioStage.running import Simulator
from scipy.optimize import minimize, Bounds, curve_fit, least_squares
from scipy.interpolate import interp1d

# Time options is a list with initial, final times and number of steps (or a list of custom time points as 4th arg)
timeOptions = [0, 5*3600, 11]
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
nRuns = 10

# Configuration of the experiment
maxDose = 1.0 # Gy
doseratefunction = None # If None, none of the rest are important (instantaneous dose)
#doseratefunction = 'exponential'
irradiationTime = 25*3600
doserate = 0.5/3600
halflife = 1*3600

# Path for the damage files
basepath = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/xray/sims/250keV.txt/'
version = '1.0'

# Read experimental data
exptimes = np.array([])
expdsb = np.array([])
with open('/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/Medras/medras_normal.csv', 'r') as csvfile:
    data_reader = csv.reader(csvfile)
    for row in data_reader:
        exptimes = np.append(exptimes, float(row[0]))
        expdsb = np.append(expdsb, float(row[1]))

sortindices = np.argsort(exptimes)
exptimes = exptimes[sortindices]
expdsb = expdsb[sortindices]

listresiduals = []
rComplex_0 = 4.222e-4

# Define the function to be optimized
def optimization_function(params):
    #rNCN, rComplex = params
    rNCN = params
    print('Trying rNCN = ' + str(rNCN))# + ', rComplex = ' + str(rComplex) )
    dsbparams = {'NEHJ': True, 'rNCNC': rNCN, 'rNCNCunits': 'rep/s', 'rComplex': rComplex_0, 'rComplexunits': 'rep/s',
                 'rMMEJ': 2.361e-6, 'rMMEJunits': 'rep/s', 'sigma': 0.25, 'sigmaUnits': 'um'}
    sim = Simulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel, bdmodel=bdModel,
                    nucleusMaxRadius=nucleusMaxRadius, irradiationTime=irradiationTime, doseratefunction=doseratefunction, doseratefunctionargs=[doserate, halflife],
                    diffusionparams=diffusionparams, dsbparams=dsbparams, ssbparams=ssbparams, bdparams=bdparams)
    sim.ReadDamage(basepath, maxDose, version)
    sim.Run(nRuns, rereadDamageForNewRuns=False, basepath=basepath, maxDose=maxDose, version=version, plot=False, verbose=1)
    sim_output = sim.avgRemainingDSBOverTime
    sim_times = sim_output.times
    sim_avgDSBremaining = sim_output.avgyvalues / sim_output.avgyvalues[0]
    sim_interp = interp1d(sim_times, sim_avgDSBremaining, kind='cubic', fill_value='extrapolate')
    simdata = sim_interp(exptimes)
    simdata = simdata[exptimes <= timeOptions[1]/3600]
    #return simdata
    expdata = expdsb[exptimes <= timeOptions[1]/3600]
    residuals = (simdata - expdata)**2
    print('Residuals: ' + str(sum(residuals)))
    return sum(residuals)

# Define the initial guesses for the parameters
rNCN_0 = 1e-3
rComplex_0 = 4e-4
initial_guess = [rNCN_0]#, rComplex_0]

# Run the optimization
lower_bounds = [0]#, 0]
upper_bounds = [np.inf]#, np. inf]
result = minimize(optimization_function, initial_guess, bounds=Bounds(lower_bounds, upper_bounds))
#popt, pcov = curve_fit(optimization_function, exptimes, expdsb, p0=initial_guess)#, bounds=(lower_bounds, upper_bounds))

# Print the optimal values of the parameters
#print(popt)
print("Optimal values of rNCN and rComplex: ", result.x)
dsbparams = {'NEHJ': True, 'rNCNC': result.x[0], 'rNCNCunits': 'rep/s', 'rComplex': rComplex_0, 'rComplexunits': 'rep/s',
             'rMMEJ': 2.361e-6, 'rMMEJunits': 'rep/s', 'sigma': 0.25, 'sigmaUnits': 'um'}
sim = Simulator(timeOptions=timeOptions, diffusionmodel=diffusionModel, dsbmodel=dsbModel, ssbmodel=ssbModel,
                bdmodel=bdModel,
                nucleusMaxRadius=nucleusMaxRadius, irradiationTime=irradiationTime, doseratefunction=doseratefunction,
                doseratefunctionargs=[doserate, halflife],
                diffusionparams=diffusionparams, dsbparams=dsbparams, ssbparams=ssbparams, bdparams=bdparams)
sim.ReadDamage(basepath, maxDose, version)
sim.Run(nRuns, rereadDamageForNewRuns=False, basepath=basepath, maxDose=maxDose, version=version, verbose=1)
#plt.scatter(exptimes, expdsb, label='Experimental data')
