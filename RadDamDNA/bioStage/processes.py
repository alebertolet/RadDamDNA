#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:00 PM

@author: alejandrobertolet
"""
import numpy as np
import scipy.stats as stats
from RadDamDNA.bioStage.models import DiffusionModel, DSBRepairModel, SSBRepairModel, BDRepairModel

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
    def __init__(self, model='free', pars=None):
        if model == 'free':
            self.ActivateFreeDiffusion(pars)

    def ActivateFreeDiffusion(self, pars):
        self.model = DiffusionModel('free', pars)

    def Diffuse(self, pos, timestep):
        if self.model.Model == 'free':
            pos = np.array(pos)
            Ns = np.random.normal(0, 1, 3)
            return pos + Ns * np.sqrt(self.model.diffusionCoefficient * timestep)
        else:
            return pos

class DSBRepair(TrackPairProcess):
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateDSBRepairStandard(pars)
        elif model == 'lognormtimes':
            self.ActivateDSBRepairLogNormTimes(pars)

    def ActivateDSBRepairStandard(self, pars=None):
        self.model = DSBRepairModel('standard', pars)
        self.InteractionRadius = 1e8

    def ActivateDSBRepairLogNormTimes(self, pars=None):
        self.model = DSBRepairModel('lognormtimes', pars)
        self.InteractionRadius = 1e8

    def Repair(self, tracklist, timestep):
        ntracks = len(tracklist)
        probmatrix = np.zeros([ntracks, ntracks])
        for i in range(ntracks):
            for j in range(i+1, ntracks):
                if tracklist[i].GetLastStep().Status == DAMAGED and tracklist[j].GetLastStep().Status == DAMAGED:
                    probmatrix[i, j] = self.GetStandardPairwiseProbability(tracklist[i], tracklist[j], timestep)
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
                for k in range(ntracks):
                    repaired[repj, k] = False
                for j in range(i+1, ntracks):
                    if j != repj:
                        repaired[i, j] = False
        return repaired

    def GetStandardPairwiseProbability(self, betrack1, betrack2, timestep):
        if self.model.Model == 'standard':
            distance = self._getDistance(betrack1, betrack2)
            if distance > self.InteractionRadius:
                return 0.0
            else:
                fd = np.exp(-distance**2 / (2*self.model.sigmaDistance**2))
                if self.model.competentInNEHJ:
                    if betrack1.GetLastStep().Complexity < 10.03 and betrack2.GetLastStep().Complexity < 10.03:
                        return (1.0 - np.exp(-self.model.repairRateNCNC * timestep)) * fd
                    else:
                        return (1.0 - np.exp(-self.model.repairRateComplex * timestep)) * fd
                else:
                    return (1.0 - np.exp(-self.model.repairMMEJ * timestep)) * fd
        elif self.model.Model == 'lognormtimes':
            distance = self._getDistance(betrack1, betrack2)
            if distance > self.InteractionRadius:
                return 0.0
            else:
                fd = np.exp(-distance**2 / (2*self.model.sigmaDistance**2))
                complexity = int(np.round((betrack1.GetLastStep().Complexity + betrack2.GetLastStep().Complexity) / 2))
                dist = stats.lognorm(self.model.sigma[complexity-2], scale=self.model.mu[complexity-2])
                prob = dist.cdf(timestep) * fd
                return prob

    @property
    def CompetentInNEHJ(self):
        return self._nhej
    @CompetentInNEHJ.setter
    def CompetentInNEHJ(self, v):
        self._nhej = v

class SSBRepair:
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateSSBRepairStandard(pars)

    def ActivateSSBRepairStandard(self, pars):
        self.model = SSBRepairModel('standard', pars)

    def Repair(self, damtrack, timestep):
        if self.model.Model == 'standard':
            if damtrack.Complexity >= 10.0:
                prob = 1 - np.exp(-self.model.repairRateComplex * timestep)
            else:
                prob = 1 - np.exp(-self.model.repairRateNoComplex * timestep)
            r = np.random.random()
            return r <= prob

class BDRepair:
    def __init__(self, model='standard', pars=None):
        if model == 'standard':
            self.ActivateBDRepairStandard(pars)

    def ActivateBDRepairStandard(self, pars):
        self.model = BDRepairModel('standard', pars)

    def Repair(self, timestep):
        if self.model.Model == 'standard':
            prob = 1 - np.exp(-self.model.repairRate * timestep)
            r = np.random.random()
            return r <= prob