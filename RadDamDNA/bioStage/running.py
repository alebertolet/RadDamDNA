#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:15 PM

@author: alejandrobertolet
"""

import numpy as np

from RadDamDNA.bioStage import processes
from RadDamDNA.bioStage import tracking
from RadDamDNA import damage

class Simulator:
    def __init__(self, originalDamage, timeOptions=[], diffusionmodel='free', nucleusMaxRadius = None):
        self.runManager = RunManager(originalDamage, timeOptions, diffusionmodel, nucleusMaxRadius)
        self.runManager.Run()

class RunManager:
    def __init__(self, dam, timeOptions = [], diffusionmodel='free', nucleusMaxRadius = None):
        if diffusionmodel == 'free':
            self.DiffusionActivated = True
            self.diffusionModel = processes.Diffusion()
        if len(timeOptions) > 3:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2], timeOptions[3])
        else:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2])
        self.nucleusMaxRadius = nucleusMaxRadius
        self.InitializeBeTracks(dam)

    def InitializeBeTracks(self, dam):
        trackid = 0
        self.betracks = []
        for iCh in dam.DSBMap:
            for iBp in dam.DSBMap[iCh]:
                for iCo in dam.DSBMap[iCh][iBp]:
                    if dam.DSBMap[iCh][iBp][iCo].type > 0:
                        pos = dam.DSBMap[iCh][iBp][iCo].position
                        time = dam.DSBMap[iCh][iBp][iCo].particletime
                        newBeStep = tracking.BeStep(pos, time)
                        newBeTrack = tracking.BeTrack(trackid)
                        newBeTrack.AddNewStep(newBeStep)
                        self.betracks.append(newBeTrack)
                        trackid += 1

    def Run(self):
        while self.clock.CurrentTime != self.clock.FinalTime:
            self.clock.AdvanceTimeStep()
            if self.DiffusionActivated:
                self.DoDiffusion()

    def DoOneStep(self):
        if self.DiffusionActivated:
            self.DoDiffusion()

    def DoDiffusion(self):
        for i, t in enumerate(self.betracks):
            newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            while self._checkPosWithinNucleus(newpos) is False:
                newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            newstep = tracking.BeStep(newpos, self.clock.CurrentTime)
            self.betracks[i].AddNewStep(newstep)

    def _checkPosWithinNucleus(self, pos):
        if self.nucleusMaxRadius is None:
            return True
        else:
            pos = np.array(pos)
            if np.sqrt(np.sum(np.power(pos, 2))) > self.nucleusMaxRadius:
                return False

    @property
    def DiffusionActivated(self):
        if self._diffusionactivated is False:
            return self._diffusionactivated
        else:
            return True
    @DiffusionActivated.setter
    def DiffusionActivated(self, b):
        self._diffusionactivated = b

    def _getDistance(self, betrack1, betrack2):
        pos1 = betrack1.GetLastStep().Position
        pos2 = betrack2.GetLastStep().Position
        return np.sqrt(np.power(pos1[0] - pos2[0], 2) + np.power(pos1[1] - pos2[1], 2) + np.power(pos1[2] - pos2[2], 2))

class Clock:
    def __init__(self, initialTime, finalTime, nSteps, listOfTimePoints = None):
        self.CurrentIndex = 0
        if listOfTimePoints is None:
            self.timepoints = np.linspace(initialTime, finalTime, nSteps)
        else:
            self.timepoints = listOfTimePoints

    def AdvanceTimeStep(self):
        self.CurrentIndex = self._currentindex + 1

    @property
    def CurrentIndex(self):
        return self._currentindex
    @CurrentIndex.setter
    def CurrentIndex(self, i):
        self._currentindex = i

    @property
    def CurrentTime(self):
        return self.timepoints[self.CurrentIndex]

    @property
    def CurrentTimeStep(self):
        if self.CurrentIndex < len(self.timepoints) - 1:
            return self.timepoints[self.CurrentIndex + 1] - self.timepoints[self.CurrentIndex]
        else:
            return self.timepoints[-1] - self.timepoints[-2]

    @property
    def InitialTime(self):
        return self.timepoints[0]

    @property
    def FinalTime(self):
        return self.timepoints[-1]

