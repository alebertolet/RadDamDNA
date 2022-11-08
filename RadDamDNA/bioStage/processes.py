#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:00 PM

@author: alejandrobertolet
"""
import numpy as np

class TrackPairProcess:
    def __init__(self):
        pass

    @property
    def InteractionRadius(self):
        return self._intradius
    @InteractionRadius.setter
    def InteractionRadius(self, r):
        self._intradius = r

class Diffusion:
    def __init__(self, model='free'):
        if model == 'free':
            self.ActivateFreeDiffusion()

    def ActivateFreeDiffusion(self):
        self.model = 0
        self.Dcoeff = 1.7e-4 # um^2/s. Diffusion coefficient used by PARTRAC

    def Diffuse(self, pos, timestep):
        if self.model == 0:
            pos = np.array(pos)
            Ns = np.random.normal(0, 1, 3)
            return pos + Ns * np.sqrt(self.Dcoeff * timestep)
        else:
            return pos

class DSBRepair(TrackPairProcess):
    def __init__(self, model='standard'):
        if model == 'standard':
            self.ActivateDSBRepairStandard()

    def ActivateDSBRepairStandard(self):
        self.model = 0
        self.InteractionRadius = 1e8

    def Repair(self, betrack1, betrack2):
        pass

class SSBRepair:
    def __init__(self, model='standard'):
        if model == 'standard':
            self.ActivateSSBRepairStandard()

    def ActivateSSBRepairStandard(self):
        self.model = 0

    def Repair(self, damtrack, timestep):
        if self.model == 0:
            # Two repair rates
            self.repairRateNoComplex = 1.774e-3 # rep/s # Fitted to data from Schiplers and Iliakis 2013
            self.repairRateComplex = 2.247e-16 # Complex SSB are way more difficult to repair
            if damtrack.Complexity >= 10.0:
                prob = 1 - np.exp(-self.repairRateComplex * timestep)
            else:
                prob = 1 - np.exp(-self.repairRateNoComplex * timestep)
            r = np.random.random()
            return 2 * (r <= prob)

class BDRepair:
    def __init__(self, model='standard'):
        if model == 'standard':
            self.ActivateBDRepairStandard()

    def ActivateBDRepairStandard(self):
        self.model = 0

    def Repair(self, damtrack, timestep):
        if self.model == 0:
            # Single exponential repair
            self.repairRate = 4.297e-4  # Fitted data from Rahmanian, Taleei and Nikjoo 2014
            prob = 1 - np.exp(-self.repairRate * timestep)
            r = np.random.random()
            return 2 * (r <= prob)