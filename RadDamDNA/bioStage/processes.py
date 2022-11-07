#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:00 PM

@author: alejandrobertolet
"""
import numpy as np

class TrackPairProcess:
    def __init__(self, trackPair):
        self.TrackPair = trackPair

    @property
    def InteractionRadius(self):
        return self._intradius
    @InteractionRadius.setter
    def InteractionRadius(self, r):
        self._intradius = r

    @property
    def InteractionRate(self):
        return self._intrate
    @InteractionRate.setter
    def InteractionRate(self, r):
        self._intrate = r

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

class Repair(TrackPairProcess):
    def __init__(self, trackPair):
        super().__init__(trackPair)

class Misrepair(TrackPairProcess):
    def __init__(self, trackPair):
        super().__init__(trackPair)