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
    def __init__(self, dam, timeOptions = [], diffusionmodel='free', nucleusMaxRadius = None, dsbrepairmodel='standard', ssbrepairmodel='standard', bdrepairmodel='standard'):
        if diffusionmodel == 'free':
            self.DiffusionActivated = True
            self.diffusionModel = processes.Diffusion()
        if dsbrepairmodel == 'standard':
            self.DSBRepairActivated = True
            self.dsbRepairModel = processes.DSBRepair()
        if ssbrepairmodel == 'standard':
            self.SSBRepairActivated = True
            self.ssbRepairModel = processes.SSBRepair()
        if bdrepairmodel == 'standard':
            self.BDRepairActivated = True
            self.bdRepairModel = processes.BDRepair()
        if len(timeOptions) > 3:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2], timeOptions[3])
        else:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2])
        self.nucleusMaxRadius = nucleusMaxRadius
        self.damage = dam
        self.InitializeTracks(dam)

    def InitializeTracks(self, dam):
        trackid = 0
        self.betracks = []
        self.ssbdamages = []
        self.bdamages = []
        for iCh in dam.DSBMap:
            for iBp in dam.DSBMap[iCh]:
                for iCo in dam.DSBMap[iCh][iBp]:
                    if dam.DSBMap[iCh][iBp][iCo].type > 0:
                        pos = dam.DSBMap[iCh][iBp][iCo].position
                        time = dam.DSBMap[iCh][iBp][iCo].particletime
                        newBeStep = tracking.BeStep(pos, time, complexity=dam.DSBMap[iCh][iBp][iCo].complexity)
                        newBeTrack = tracking.BeTrack(trackid)
                        newBeTrack.ChromosomeID = iCh
                        newBeTrack.BasePairID = iBp
                        newBeTrack.StrandID = iCo
                        newBeTrack.AddNewStep(newBeStep)
                        self.betracks.append(newBeTrack)
                        trackid += 1
        for iCh in dam.SSBMap:
            for iBp in dam.SSBMap[iCh]:
                for iCo in dam.SSBMap[iCh][iBp]:
                    if dam.SSBMap[iCh][iBp][iCo].type > 0:
                        pos = dam.SSBMap[iCh][iBp][iCo].position
                        time = dam.SSBMap[iCh][iBp][iCo].particletime
                        newSSBDamage = tracking.DamageTrack(trackid)
                        newSSBDamage.Time = time
                        newSSBDamage.Position = pos
                        newSSBDamage.ChromosomeID = iCh
                        newSSBDamage.BasePairID = iBp
                        newSSBDamage.StrandID = iCo
                        newSSBDamage.Complexity = dam.SSBMap[iCh][iBp][iCo].complexity
                        self.ssbdamages.append(newSSBDamage)
                        trackid += 1
        for iCh in dam.BDMap:
            for iBp in dam.BDMap[iCh]:
                for iCo in dam.BDMap[iCh][iBp]:
                    if dam.BDMap[iCh][iBp][iCo].type > 0:
                        pos = dam.BDMap[iCh][iBp][iCo].position
                        time = dam.BDMap[iCh][iBp][iCo].particletime
                        newBDDamage = tracking.DamageTrack(trackid)
                        newBDDamage.Time = time
                        newBDDamage.Position = pos
                        newBDDamage.ChromosomeID = iCh
                        newBDDamage.BasePairID = iBp
                        newBDDamage.StrandID = iCo
                        newBDDamage.Complexity = dam.BDMap[iCh][iBp][iCo].complexity
                        self.bdamages.append(newBDDamage)
                        trackid += 1

    def Run(self):
        while self.clock.CurrentTime != self.clock.FinalTime:
            self.clock.AdvanceTimeStep()
            self.DoOneStep()

    def DoOneStep(self):
        if self.DiffusionActivated:
            self.DoDiffusion()
        if self.DSBRepairActivated:
            self.DoDSBRepair()
        if self.SSBRepairActivated:
            self.DoSSBRepair()
        if self.BDRepairActivated:
            self.DoBDRepair()
        self.UpdateDamageMaps()
        self.damage.printDamageCount()

    def DoDiffusion(self):
        for i, t in enumerate(self.betracks):
            newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            while self._checkPosWithinNucleus(newpos) is False:
                newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            newstep = tracking.BeStep(newpos, self.clock.CurrentTime)
            self.betracks[i].AddNewStep(newstep)

    def DoDSBRepair(self):
        for i in range(len(self.betracks)):
            if self.betracks[i].GetLastStep().Status == 1:
                for j in range(i + 1, len(self.betracks)):
                    distance = self._getDistance(self.betracks[i], self.betracks[j])
                    if distance <= self.dsbRepairModel.InteractionRadius:
                        newstatus = self.dsbRepairModel.Repair(self.betracks[i], self.betracks[j])
                        self.betracks[i].GetLastStep().Status = newstatus
                        self.betracks[j].GetLastStep().Status = newstatus

    def DoSSBRepair(self):
        for i in range(len(self.ssbdamages)):
            if self.ssbdamages[i].Status == 1:
                newstatus = self.ssbRepairModel.Repair(self.ssbdamages[i], self.clock.CurrentTimeStep)
                self.ssbdamages[i].Status = newstatus

    def DoBDRepair(self):
        for i in range(len(self.bdamages)):
            if self.bdamages[i].Status == 1:
                newstatus = self.bdRepairModel.Repair(self.bdamages[i], self.clock.CurrentTimeStep)
                self.bdamages[i].Status = newstatus

    def UpdateDamageMaps(self):
        for be in self.betracks:
            if be.GetLastStep().Status == 2:
                iCh = be.ChromosomeID
                iBp = be.BasePairID
                iCo = be.StrandID + 1
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        for ssb in self.ssbdamages:
            if ssb.Status == 2:
                iCh = ssb.ChromosomeID
                iBp = ssb.BasePairID
                iCo = ssb.StrandID + 1
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        for bd in self.bdamages:
            if bd.Status == 2:
                iCh = bd.ChromosomeID
                iBp = bd.BasePairID
                iCo = bd.StrandID
                if iCo == 2:
                    iCo = 4
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        self.damage.recomputeDamagesFromReadSites(stopAtTime=self.clock.CurrentTime)

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

    @property
    def DSBRepairActivated(self):
        if self._dsbrepactivated is False:
            return self._dsbrepactivated
        else:
            return True
    @DSBRepairActivated.setter
    def DSBRepairActivated(self, b):
        self._dsbrepactivated = b

    @property
    def SSBRepairActivated(self):
        if self._ssbrepactivated is False:
            return self._ssbrepactivated
        else:
            return True
    @SSBRepairActivated.setter
    def SSBRepairActivated(self, b):
        self._ssbrepactivated = b

    @property
    def BDRepairActivated(self):
        if self._bdrepactivated is False:
            return self._bdrepactivated
        else:
            return True
    @BDRepairActivated.setter
    def BDRepairActivated(self, b):
        self._bdrepactivated = b

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

