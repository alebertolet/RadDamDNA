#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 1/22/23 2:47 PM

@author: alejandrobertolet
"""

class DiffusionModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'free':
            if pars is None:
                self.diffusionCoefficient = 5.0e-8
                self.diffusionCoefficientUnits = 'um^2/s'
            else:
                self.diffusionCoefficient = pars['D']
                try:
                    self.diffusionCoefficientUnits = pars['Dunits']
                except:
                    self.diffusionCoefficientUnits = 'um^2/s'


class DSBRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None:
                self.competentInNEHJ = True
                self.repairRateNCNC = 2.0e-4
                self.repairRateNCNCUnits = 'rep/s'
                self.repairRateComplex = 1.0e-5
                self.repairRateComplexUnits = 'rep/s'
                self.repairMMEJ = 2.361e-7
                self.repairMMEJUnits = 'rep/s'
                self.sigmaDistance = 0.25
                self.sigmaDistanceUnits = 'um'
            else:
                self.competentInNEHJ = pars['NEHJ']
                self.repairRateNCNC = pars['rNCNC']
                try:
                    self.repairRateNCNCUnits = pars['rNCNCunits']
                except:
                    self.repairRateNCNCUnits = 'rep/s'
                self.repairRateComplex = pars['rComplex']
                try:
                    self.repairRateComplexUnits = pars['rComplexunits']
                except:
                    self.repairRateComplexUnits = 'rep/s'
                self.repairMMEJ = pars['rMMEJ']
                try:
                    self.repairMMEJUnits = pars['rMMEJunits']
                except:
                    self.repairMMEJUnits = 'rep/s'
                self.sigmaDistance = pars['sigma']
                try:
                    self.sigmaDistanceUnits = pars['sigmaUnits']
                except:
                    self.sigmaDistanceUnits = 'um'

class SSBRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None:
                self.repairRateNoComplex = 1.774e-3
                self.repairRateNoComplexUnits = 'rep/s'
                self.repairRateComplex = 2.247e-16
                self.repairRateComplexUnits = 'rep/s'
            else:
                self.repairRateNoComplex = pars['rNC']
                try:
                    self.repairRateNoComplexUnits = pars['rNCunits']
                except:
                    self.repairRateNoComplexUnits = 'rep/s'
                self.repairRateComplex = pars['rC']
                try:
                    self.repairRateComplexUnits = pars['rCunits']
                except:
                    self.repairRateComplexUnits = 'rep/s'

class BDRepairModel:
    def __init__(self, name, pars=None):
        self.Model = name
        self.SetParameters(pars)

    def SetParameters(self, pars=None):
        # pars is a dictionary with the parameters
        if self.Model == 'standard':
            if pars is None:
                self.repairRate = 1.774e-3
                self.repairRateUnits = 'rep/s'
            else:
                self.repairRate = pars['r']
                try:
                    self.repairRateUnits = pars['runits']
                except:
                    self.repairRateUnits = 'rep/s'