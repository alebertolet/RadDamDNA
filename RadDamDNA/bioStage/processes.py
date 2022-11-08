#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:00 PM

@author: alejandrobertolet
"""
import numpy as np

DAMAGED = 1
REPAIRED = 2
MISREPAIRED = 3

class TrackPairProcess:
    def __init__(self):
        pass

    @property
    def InteractionRadius(self):
        return self._intradius
    @InteractionRadius.setter
    def InteractionRadius(self, r):
        self._intradius = r

    def _getDistance(self, betrack1, betrack2):
        pos1 = betrack1.GetLastStep().Position
        pos2 = betrack2.GetLastStep().Position
        return np.sqrt(np.power(pos1[0] - pos2[0], 2) + np.power(pos1[1] - pos2[1], 2) + np.power(pos1[2] - pos2[2], 2))

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

    def Repair(self, tracklist, timestep):
        if self.model == 0:
            self.CompetentInNEHJ = True
            self.distanceRegularization = 1.0  # um
            # Rates for different complexity combinations
            # Fast repair for not complex breaks
            self.repairRateNCNC = 5.833e-5 # rep/s/um^2 # MEDRAS parameter adapted by distance
            # Slow repair for complex breaks
            self.repairRateComplex = 7.222e-7 # rep/s # MEDRAS parameter adapted by distance
            # Repair for when NHEJ is non active
            self.repairMMEJ = 2.361e-6 # rep/s # MEDRAS parameter adapted by distance
            ntracks = len(tracklist)
            probmatrix = np.zeros([ntracks, ntracks])
            for i in range(ntracks):
                for j in range(i+1, ntracks):
                    if tracklist[i].GetLastStep().Status == DAMAGED and tracklist[j].GetLastStep().Status == DAMAGED:
                        probmatrix[i, j] = self.GetPairwiseProbability(tracklist[i], tracklist[j], timestep)
                    else:
                        probmatrix[i, j] = 0.0
            rs = np.random.random([ntracks, ntracks])
            rs[probmatrix == 0.0] = 1.0
            repaired = rs <= probmatrix
            # Avoid multiple uses of the same and: pick from maximum to minimum probability
            for i in range(ntracks):
                maxFori = 0.0
                repj = 0
                for j in range(i+1, ntracks):
                    if repaired[i, j] and probmatrix[i, j] > maxFori:
                        maxFori = probmatrix[i, j]
                        repj = j
                if repj > 0:
                    for j in range(i+1, ntracks):
                        if j != repj:
                            repaired[i, j] = False
            return repaired

    def GetPairwiseProbability(self, betrack1, betrack2, timestep):
        if self.model == 0:
            distance = self._getDistance(betrack1, betrack2)
            if distance > self.InteractionRadius:
                return 0.0
            else:
                if distance < self.distanceRegularization:
                    distance = self.distanceRegularization
                fd = np.power(self.distanceRegularization/distance, 2)
                if self.CompetentInNEHJ:
                    if betrack1.GetLastStep().Complexity == 10.0 and betrack2.GetLastStep().Complexity == 10:
                        return (1.0 - np.exp(-self.repairRateNCNC * timestep)) * fd
                    else:
                        return (1.0 - np.exp(-self.repairRateComplex * timestep)) * fd
                else:
                    return (1.0 - np.exp(-self.repairMMEJ * timestep)) * fd

    @property
    def CompetentInNEHJ(self):
        return self._nhej
    @CompetentInNEHJ.setter
    def CompetentInNEHJ(self, v):
        self._nhej = v

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
            return r <= prob

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
            return r <= prob