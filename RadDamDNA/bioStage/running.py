#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 3:15 PM

@author: alejandrobertolet
"""

import numpy as np
from copy import deepcopy, copy

from RadDamDNA.bioStage import processes
from RadDamDNA.bioStage import tracking
from RadDamDNA.bioStage import output
from RadDamDNA.damage import DamageToDNA

import random
import os

DAMAGED = 1
REPAIRED = 2
MISREPAIRED = 3

# Outputs
DSB = 0
MISREPDSB = 1
SSB = 2
BD = 3

class Simulator:
    def __init__(self, timeOptions=[], diffusionmodel='free', dsbmodel='standard', ssbmodel='standard', bdmodel='standard', nucleusMaxRadius = None,
                 irradiationTime=0, doseratefunction=None, doseratefunctionargs=None, diffusionparams=None, dsbparams=None, ssbparams=None, bdparams=None):
        self.messages = []
        self.irradiationTime = irradiationTime
        self.doseratefunction = doseratefunction
        self.doseratefunctionargs = doseratefunctionargs
        self.runManager = RunManager(timeOptions=timeOptions, diffusionmodel=diffusionmodel, dsbrepairmodel=dsbmodel, ssbrepairmodel=ssbmodel,
                                     bdrepairmodel=bdmodel, nucleusMaxRadius=nucleusMaxRadius, messages=self.messages, diffusionParameters=diffusionparams,
                                     dsbrepairParameters=dsbparams, ssbrepairParameters=ssbparams, bdrepairParameters=bdparams)

    def Run(self, nRuns, rereadDamageForNewRuns=True, basepath=None, maxDose=-1, version=None, plot=True, outputnorm=True, verbose=0):
        self.nRuns = nRuns
        self.runManager.nRuns = nRuns
        self.runManager.maxDose = maxDose
        self.runManager.plotflag = False
        if plot:
            self.runManager.plotflag = True
        if not rereadDamageForNewRuns:
            self.runManager.Run(verbose=verbose, outputnorm=outputnorm)
        else:
            self.runManager.TotalRuns = self.nRuns
            self.runManager.plotflag = False
            for i in range(self.nRuns):
                self.runManager.nRuns = 1
                if i == self.nRuns - 1 and plot:
                    self.runManager.plotflag = True
                if i == 0:
                    self.runManager.Run(verbose=verbose, outputnorm=outputnorm)
                else:
                    self.ReadDamage(basepath, maxDose, version)
                    self.runManager.Run(verbose=verbose, outputnorm=outputnorm)
        self.avgRemainingDSBOverTime = self.runManager.runoutputDSB
        self.avgMisrepairedDSBOverTime = self.runManager.runoutputMisrepairedDSB

    def ReadDamage(self, basepath, maxDose=2.0, version='2.0', recalculatePerEachTrack=False):
        damage = DamageToDNA(messages=self.messages)
        nfiles = len(os.listdir(basepath))
        # Section to get what directories actually contains both dose and SDD. Disregard others!
        listOfAvailableDirs = []
        for j in range(nfiles):
            newpath = basepath + str(j) + '/'
            if str(j) in os.listdir(basepath):
                files = os.listdir(newpath)
                if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
                    if os.path.getsize(newpath + 'DNADamage_sdd.txt') > 0:  # only those with actual data
                        listOfAvailableDirs.append(j)
        neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
        accumulatedose = 0
        for i, e in enumerate(neworder):
            accumulatedose = damage.accumulateDose
            time = self._getTimeForDose(accumulatedose)
            if 0 < self.irradiationTime < time:
                time = 1e20
            path = basepath + str(e) + '/'
            damage.readSDDAndDose(path, version=version, particleTime=time, lesionTime=time)
        damage.populateDamages(getVideo=False, stopAtDose=maxDose, stopAtTime=0, recalculatePerEachTrack=recalculatePerEachTrack)
        self.runManager.damage = damage

    def _getTimeForDose(self, d):
        if self.irradiationTime == 0:
            return 0
        if self.doseratefunction is None:
            return 0
        if self.doseratefunction == 'uniform':
            return d / self.doseratefunctionargs[0]
        if self.doseratefunction == 'linear':
            return (-self.doseratefunction[0] + np.sqrt(self.doseratefunction[0] ** 2 + 4 * self.doseratefunction[1] * d)) / (2 * self.doseratefunction[1])
        if self.doseratefunction == 'exponential':
            constant = np.log(2) / self.doseratefunctionargs[1]
            initialdoserate = self.doseratefunctionargs[0]
            return -1/constant * np.log(1-d/(initialdoserate/constant))

class RunManager:
    def __init__(self, timeOptions = [], diffusionmodel='free', dsbrepairmodel='standard', ssbrepairmodel='standard',
                 bdrepairmodel='standard', nucleusMaxRadius = None, outputs=[DSB, MISREPDSB, SSB, BD], messages=[],
                 diffusionParameters=None, dsbrepairParameters=None, ssbrepairParameters=None, bdrepairParameters=None):
        self.messages = messages
        self.maxDose = -1
        self._diffusionactivated = False
        self._dsbrepactivated = False
        self._ssbrepactivated = False
        self._bdrepactivated = False
        self.trackid = 0
        if diffusionmodel is not None or diffusionmodel.lower() != 'none':
            self.DiffusionActivated = True
            self.diffusionModel = processes.Diffusion(diffusionmodel, diffusionParameters)
        if dsbrepairmodel is not None or dsbrepairmodel.lower() != 'none':
            self.DSBRepairActivated = True
            self.dsbRepairModel = processes.DSBRepair(dsbrepairmodel, dsbrepairParameters)
        if ssbrepairmodel is not None or ssbrepairmodel.lower() != 'none':
            self.SSBRepairActivated = True
            self.ssbRepairModel = processes.SSBRepair(ssbrepairmodel, ssbrepairParameters)
        if bdrepairmodel is not None or bdrepairmodel.lower() != 'none':
            self.BDRepairActivated = True
            self.bdRepairModel = processes.BDRepair(bdrepairmodel, bdrepairParameters)
        if len(timeOptions) > 3:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2], timeOptions[3])
        else:
            self.clock = Clock(timeOptions[0], timeOptions[1], timeOptions[2])
        self.nucleusMaxRadius = nucleusMaxRadius
        self.outputs = outputs
        self.runoutputDSB = output.AverageTimeCurveOverRuns()
        self.runoutputMisrepairedDSB = output.AverageTimeCurveOverRuns()
        self.plotflag = True
        self.currentrun = 0

    def InitializeNewRun(self):
        self.betracks = []
        self.ssbdamages = []
        self.bdamages = []

    def InitializeNewTracks(self, dam):
        for iCh in dam.DSBMap:
            for iBp in dam.DSBMap[iCh]:
                for iCo in dam.DSBMap[iCh][iBp]:
                    if dam.DSBMap[iCh][iBp][iCo].type > 0:
                        pos = dam.DSBMap[iCh][iBp][iCo].position
                        time = dam.DSBMap[iCh][iBp][iCo].particletime
                        if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                            newBeStep = tracking.BeStep(pos, time, complexity=dam.DSBMap[iCh][iBp][iCo].complexity)
                            newBeTrack = tracking.BeTrack(self.trackid, dam.DSBMap[iCh][iBp][iCo].dsbID)
                            newBeTrack.ChromosomeID = iCh
                            newBeTrack.BasePairID = iBp
                            newBeTrack.StrandID = iCo
                            newBeTrack.AddNewStep(newBeStep)
                            self.betracks.append(newBeTrack)
                            self.trackid += 1
        for iCh in dam.SSBMap:
            for iBp in dam.SSBMap[iCh]:
                for iCo in dam.SSBMap[iCh][iBp]:
                    if dam.SSBMap[iCh][iBp][iCo].type > 0:
                        pos = dam.SSBMap[iCh][iBp][iCo].position
                        time = dam.SSBMap[iCh][iBp][iCo].particletime
                        if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                            newSSBDamage = tracking.DamageTrack(self.trackid)
                            newSSBDamage.Time = time
                            newSSBDamage.Position = pos
                            newSSBDamage.ChromosomeID = iCh
                            newSSBDamage.BasePairID = iBp
                            newSSBDamage.StrandID = iCo
                            newSSBDamage.Complexity = dam.SSBMap[iCh][iBp][iCo].complexity
                            self.ssbdamages.append(newSSBDamage)
                            self.trackid += 1
        for iCh in dam.BDMap:
            for iBp in dam.BDMap[iCh]:
                for iCo in dam.BDMap[iCh][iBp]:
                    if dam.BDMap[iCh][iBp][iCo].type > 0:
                        pos = dam.BDMap[iCh][iBp][iCo].position
                        time = dam.BDMap[iCh][iBp][iCo].particletime
                        if time <= self.clock.CurrentTime and time > self.clock.CurrentTime - self.clock.CurrentTimeStep:
                            newBDDamage = tracking.DamageTrack(self.trackid)
                            newBDDamage.Time = time
                            newBDDamage.Position = pos
                            newBDDamage.ChromosomeID = iCh
                            newBDDamage.BasePairID = iBp
                            newBDDamage.StrandID = iCo
                            newBDDamage.Complexity = dam.BDMap[iCh][iBp][iCo].complexity
                            self.bdamages.append(newBDDamage)
                            self.trackid += 1

    def Run(self, verbose=0, outputnorm=True):
        self.originaldamage = deepcopy(self.damage)
        for i in range(self.nRuns):
            self.InitializeNewRun()
            self.repairedList = []
            self.misrepairedlist = []
            self.currentrun += 1
            self.messages.append('Repair run ' + str(self.currentrun) + ' of ' + str(self.TotalRuns) + '...')
            if verbose > 0:
                print(self.messages[-1])
            if DSB in self.outputs:
                self.outDSB = output.TimeCurveForSingleRun('Remaining DSB')
                self.misrepDSB = output.TimeCurveForSingleRun('Misrepaired DSB')
                self.DSBEvolution()
            while self.clock.CurrentTime != self.clock.FinalTime:
                if verbose > 1:
                    print("Time " + str(round(self.clock.CurrentTime/3600,2)) + " h - Dose: " + str(round(self.damage.cumulativeDose, 2)) + " Gy. Number of DSB: " + str(self.damage.numDSB))
                self.InitializeNewTracks(self.damage)
                self.clock.AdvanceTimeStep()
                self.DoOneStep()
            if DSB in self.outputs:
                self.runoutputDSB.AddTimeCurveForSingleRun(self.outDSB)
                self.runoutputMisrepairedDSB.AddTimeCurveForSingleRun(self.misrepDSB)
                #self.outDSB.Plot()
                #self.outDSB.WriteCSV()
            self.messages.append('Repaired: ' + str(len(self.repairedList)) + ' - Misrepaired: ' + str(len(self.misrepairedlist)))
            if verbose > 0:
                print(self.messages[-1])
            self.clock.Reset()
            self.resetDamage()
        if DSB in self.outputs:
            self.runoutputDSB.DoStatistics(outputnorm)
            self.runoutputMisrepairedDSB.DoStatistics(False)
            if self.plotflag:
                self.runoutputDSB.Plot()
                self.runoutputDSB.WriteCSV()
                self.runoutputMisrepairedDSB.Plot()
                self.runoutputMisrepairedDSB.WriteCSV()


    def resetDamage(self):
        self.damage = deepcopy(self.originaldamage)

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
        if DSB in self.outputs:
            self.DSBEvolution()

    def DoDiffusion(self):
        for i, t in enumerate(self.betracks):
            newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            while self._checkPosWithinNucleus(newpos) is False:
                newpos = self.diffusionModel.Diffuse(t.GetLastStep().Position, self.clock.CurrentTimeStep)
            newstep = tracking.BeStep(newpos, self.clock.CurrentTime, complexity=self.betracks[i].GetLastStep().Complexity)
            newstep.Status = self.betracks[i].GetLastStep().Status
            self.betracks[i].AddNewStep(newstep)

    def DoDSBRepair(self):
        repair = self.dsbRepairModel.Repair(self.betracks, self.clock.CurrentTimeStep)
        for i in range(repair.shape[0]):
            for j in range(i+1, repair.shape[1]):
                if repair[i, j]:
                    if self.betracks[i].DSBid == self.betracks[j].DSBid:
                        self.betracks[i].GetLastStep().Status = REPAIRED
                        self.betracks[j].GetLastStep().Status = REPAIRED
                        self.repairedList.append([self.betracks[i], self.betracks[j]])
                    else:
                        self.betracks[i].GetLastStep().Status = MISREPAIRED
                        self.betracks[j].GetLastStep().Status = MISREPAIRED
                        self.misrepairedlist.append([self.betracks[i], self.betracks[j]])

    def DoSSBRepair(self):
        for i in range(len(self.ssbdamages)):
            if self.ssbdamages[i].Status == DAMAGED:
                if self.ssbRepairModel.Repair(self.ssbdamages[i], self.clock.CurrentTimeStep):
                    self.ssbdamages[i].Status = REPAIRED

    def DoBDRepair(self):
        for i in range(len(self.bdamages)):
            if self.bdamages[i].Status == DAMAGED:
                if self.bdRepairModel.Repair(self.clock.CurrentTimeStep):
                    self.bdamages[i].Status = REPAIRED

    def UpdateDamageMaps(self):
        for be in self.betracks:
            if be.GetLastStep().Status == REPAIRED or be.GetLastStep().Status == MISREPAIRED:
                iCh = be.ChromosomeID
                iBp = be.BasePairID
                iCo = be.StrandID + 1
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        for ssb in self.ssbdamages:
            if ssb.Status == REPAIRED:
                iCh = ssb.ChromosomeID
                iBp = ssb.BasePairID
                iCo = ssb.StrandID + 1
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        for bd in self.bdamages:
            if bd.Status == REPAIRED:
                iCh = bd.ChromosomeID
                iBp = bd.BasePairID
                iCo = bd.StrandID
                if iCo == 2:
                    iCo = 4
                self.damage.damageMap[iCh][iBp][iCo].type = 0
        self.damage.recomputeDamagesFromReadSites(stopAtTime=self.clock.CurrentTime, stopAtDose=self.maxDose)

    def DSBEvolution(self):
        self.outDSB.AddTimePoint(self.clock.CurrentTime, self.damage.numDSB)
        self.misrepDSB.AddTimePoint(self.clock.CurrentTime, len(self.misrepairedlist))

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

    @property
    def TotalRuns(self):
        try:
            return self._totalruns
        except:
            return self.nRuns
    @TotalRuns.setter
    def TotalRuns(self, t):
        self._totalruns = t


class Clock:
    def __init__(self, initialTime, finalTime, nSteps, listOfTimePoints = None):
        self.CurrentIndex = 0
        if listOfTimePoints is None:
            self.timepoints = np.linspace(initialTime, finalTime, nSteps)
        else:
            self.timepoints = listOfTimePoints

    def AdvanceTimeStep(self):
        self.CurrentIndex = self._currentindex + 1

    def Reset(self):
        self.CurrentIndex = 0

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

