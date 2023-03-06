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
SURVIVAL = 4

class Simulator:
    def __init__(self, timeOptions=[], diffusionmodel='free', dsbmodel='standard', ssbmodel='standard', bdmodel='standard', nucleusMaxRadius = None,
                 irradiationTime=0, doseratefunction=None, doseratefunctionargs=None, diffusionparams=None, dsbparams=None, ssbparams=None, bdparams=None,
                 p_lethal_aberration=0.1, p_apoptosis=0.1, time_to_G1S_checkpoint=10 * 3600, time_to_G2M_checkpoint=24 * 3600):
        self.messages = []
        self.irradiationTime = irradiationTime
        self.doseratefunction = doseratefunction
        self.doseratefunctionargs = doseratefunctionargs
        self.runManager = RunManager(timeOptions=timeOptions, diffusionmodel=diffusionmodel, dsbrepairmodel=dsbmodel, ssbrepairmodel=ssbmodel,
                                     bdrepairmodel=bdmodel, nucleusMaxRadius=nucleusMaxRadius, messages=self.messages, diffusionParameters=diffusionparams,
                                     dsbrepairParameters=dsbparams, ssbrepairParameters=ssbparams, bdrepairParameters=bdparams,
                                     p_lethal_aberration=p_lethal_aberration, p_apoptosis=p_apoptosis, time_to_G1S_checkpoint=time_to_G1S_checkpoint,
                                     time_to_G2M_checkpoint=time_to_G2M_checkpoint)

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
        npaths = 1
        if type(basepath) == list:
            npaths = len(basepath)
            if type(maxDose) is not list or len(maxDose) != npaths:
                print('Doses have to match the number of radiations. The same dose will be used for all')
                if type(maxDose) is not list:
                    maxDose = [maxDose] * npaths
                else:
                    maxDose = [maxDose[0]] * npaths
        totalaccumulateddose = 0
        for ib, basepath in enumerate(basepath):
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
            for i, e in enumerate(neworder):
                accumulatedose = damage.accumulateDose - totalaccumulateddose
                if accumulatedose > maxDose[ib]:
                    break
                time = self._getTimeForDose(accumulatedose)
                if 0 < self.irradiationTime < time:
                    time = 1e20
                path = basepath + str(e) + '/'
                damage.readSDDAndDose(path, version=version, particleTime=time, lesionTime=time)
            totalaccumulateddose += accumulatedose
        damage.populateDamages(getVideo=False, stopAtDose=-1, stopAtTime=0, recalculatePerEachTrack=recalculatePerEachTrack)
        self.runManager.damage = damage

    def LoadDamageFromMGM(self, listOfDamageSites):
        self.runManager.mgmFlag = True
        self.runManager.mgmDamage = listOfDamageSites
        self.runManager.damage = None

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

    def GetSurvivalFraction(self):
        return self.runManager.GetSurvivalFraction()

class RunManager:
    def __init__(self, timeOptions = [], diffusionmodel='free', dsbrepairmodel='standard', ssbrepairmodel='standard',
                 bdrepairmodel='standard', nucleusMaxRadius = None, outputs=[DSB, MISREPDSB, SSB, BD, SURVIVAL], messages=[],
                 diffusionParameters=None, dsbrepairParameters=None, ssbrepairParameters=None, bdrepairParameters=None,
                 p_lethal_aberration=0.1, p_apoptosis=0.1, time_to_G1S_checkpoint=10 * 3600, time_to_G2M_checkpoint=24 * 3600):
        self.messages = messages
        self.maxDose = -1
        self._diffusionactivated = False
        self._dsbrepactivated = False
        self._ssbrepactivated = False
        self._bdrepactivated = False
        self.trackid = 0
        self.mgmFlag = False
        self.mgmDamage = None
        self.p_lethal_aberration = p_lethal_aberration
        self.p_apoptosis = p_apoptosis
        self.time_to_G1S_checkpoint = time_to_G1S_checkpoint
        self.time_to_G2M_checkpoint = time_to_G2M_checkpoint
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
        self.outputsurvival = []
        self.plotflag = True
        self.currentrun = 0

    def InitializeNewRun(self):
        self.betracks = []
        self.ssbdamages = []
        self.bdamages = []

    def InitializeNewTracks(self, dam):
        if not self.mgmFlag:
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
        else:
            if self.clock.CurrentTime == self.clock.InitialTime:
                for i in range(len(self.mgmDamage)):
                    # First break end. Move it randomly in a 3 nm radius
                    position = self.mgmDamage[i][0] + np.random.normal(0, 3, 3) * 1e-3
                    complexity = self.mgmDamage[i][1]
                    if self.dsbRepairModel.model == 'standard':
                        if complexity > 3:
                            complexity = 11
                        else:
                            complexity = 10
                    time = 0.0
                    newBeStep = tracking.BeStep(position, time, complexity=complexity)
                    newBeTrack = tracking.BeTrack(self.trackid, i)
                    newBeTrack.ChromosomeID = np.random.randint(0, 46)
                    newBeTrack.BasePairID = i
                    newBeTrack.StrandID = 1
                    newBeTrack.AddNewStep(newBeStep)
                    self.betracks.append(newBeTrack)
                    self.trackid += 1
                    # Second break end. Move it randomly in a 3 nm radius
                    position = self.mgmDamage[i][0] + np.random.normal(0, 3, 3) * 1e-3
                    complexity = self.mgmDamage[i][1]
                    time = 0.0
                    newBeStep = tracking.BeStep(position, time, complexity=complexity)
                    newBeTrack = tracking.BeTrack(self.trackid, i)
                    newBeTrack.ChromosomeID = np.random.randint(0, 46)
                    newBeTrack.BasePairID = i
                    newBeTrack.StrandID = 0
                    newBeTrack.AddNewStep(newBeStep)
                    self.betracks.append(newBeTrack)
                    self.trackid += 1

    def Run(self, verbose=0, outputnorm=True):
        if not self.mgmFlag:
            self.originaldamage = deepcopy(self.damage)
        for i in range(self.nRuns):
            self.InitializeNewRun()
            if self.mgmFlag:
                self.InitializeNewTracks(self.damage)
                self.CountCurrentDSB()
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
                    if not self.mgmFlag:
                        print("Time " + str(round(self.clock.CurrentTime/3600,2)) + " h - Dose: " + str(round(self.damage.cumulativeDose, 2)) + " Gy. Number of DSB: " + str(self.damage.numDSB))
                    else:
                        print("Time " + str(round(self.clock.CurrentTime/3600,2)) + " h. Number of DSB: " + str(self.numDSB))
                if not self.mgmFlag:
                    self.InitializeNewTracks(self.damage)
                self.clock.AdvanceTimeStep()
                self.DoOneStep()
            if DSB in self.outputs:
                self.runoutputDSB.AddTimeCurveForSingleRun(self.outDSB)
                self.runoutputMisrepairedDSB.AddTimeCurveForSingleRun(self.misrepDSB)
                #self.outDSB.Plot()
                #self.outDSB.WriteCSV()
            if SURVIVAL in self.outputs:
                self.outputsurvival.append(self.DetermineCellFateForRun())
                if self.outputsurvival[-1] == 0:
                    print('Cell died')
                else:
                    print('Cell survived')
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

    def DetermineCellFateForRun(self):
        # Get number of misrepaired DSB
        numMisrepDSB = int(self.misrepDSB.yvalues[-1])
        # Check if cell is dead
        for i in range(numMisrepDSB):
            if np.random.rand() < self.p_lethal_aberration:
                return 0
        # Get number of unresolved DSB at the checkpoints
        numDSBcheckpoint1 = int(self.outDSB.GetValueForTimePoint(self.time_to_G1S_checkpoint))
        for i in range(numDSBcheckpoint1):
            if np.random.rand() < self.p_apoptosis:
                return 0
        numDSBcheckpoint2 = int(self.outDSB.GetValueForTimePoint(self.time_to_G2M_checkpoint))
        for i in range(numDSBcheckpoint2):
            if np.random.rand() < self.p_apoptosis:
                return 0
        # The cell survived!
        return 1

    def GetSurvivalFraction(self):
        return np.sum(self.outputsurvival)/len(self.outputsurvival)

    def resetDamage(self):
        if not self.mgmFlag:
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

    def CountCurrentDSB(self):
        self.numDSB = 0
        for be in self.betracks:
            if be.GetLastStep().Status == DAMAGED:
                self.numDSB += 0.5
        self.numDSB = int(self.numDSB)

    def UpdateDamageMaps(self):
        if not self.mgmFlag:
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
        self.CountCurrentDSB()
        self.outDSB.AddTimePoint(self.clock.CurrentTime, self.numDSB)
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

