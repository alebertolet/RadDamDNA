#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 9/1/22 10:04 AM

@author: alejandrobertolet
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
from RadDamDNA.SDD import *
from collections import defaultdict

class DamageToDNA:
    def __init__(self, messages=[]):
        self.messages = messages
        self.initializeStructures()
        self.initializeCounters()
        self.SetUpColorMap()
        self.stopAtDose = -1.0
        self.accumulateDose = 0
        self.cumulativeDose = 0
        self.numtracks = 0
        self.damagespertrack = {}
        self.damagespertrack['DSB'] = []
        self.damagespertrack['SB'] = []
        self.damagespertrack['BD'] = []
        self.damagespertrack['Dose'] = []
        self.damagespertrack['NumSites'] = []

    def initializeStructures(self):
        self.damageSites = []
        self.damageSitesAtCurrentTime = []
        self.doses = np.array([])
        self.damageMap = {}
        self.Darray = []
        self.DSBarray = []; self.DSBdirectarray = []; self.DSBindirectarray = []; self.DSBhybridarray = []
        self.SSBarray = []; self.SSBdirectarray = []; self.SSBindirectarray = []
        self.SBarray = []; self.SBdirectarray = []; self.SBindirectarray = []
        self.BDarray = []; self.BDdirectarray = []; self.BDindirectarray = []
        self.numSitesarray = []
        self.DSBPositions = []

    def initializeCounters(self):
        self.numSB = 0; self.numSBDirect = 0; self.numSBIndirect = 0
        self.numSSB = 0; self.numSSBDirect = 0; self.numSSBIndirect = 0
        self.numDSB = 0; self.numDSBDirect = 0; self.numDSBIndirect = 0; self.numDSBHybrid = 0
        self.numBD = 0; self.numBDDirect = 0; self.numBDIndirect = 0
        self.numSSBPlus = 0; self.numDSBPlus = 0; self.numDSBComplex = 0
        self.provnsites = 0

    def readSDDAndDose(self, path, namessd = 'DNADamage_sdd.txt', namephsp = 'DNADamage.phsp', version='2.0', particleTime=0, lesionTime=0):
        if namephsp is not None:
            dosepath = path + namephsp
            f = open(dosepath, 'r')
            lines = f.readlines()
            for l in lines:
                split = l.split()
                self.doses = np.append(self.doses, float(split[1]))
                self.accumulateDose += self.doses[-1]
                self.damagespertrack['Dose'].append(self.doses[-1])
        sddpath = path + namessd
        self.readFromSDD(sddpath, version, particleTime, lesionTime)
        f.close()
        self.namessd = namessd
        self.namephsp = namephsp

    def readFromSDD(self, path, version='2.0', particleTime=0, lesionTime=0):
        reader = SDDReader(path)
        for key in reader.headerProperties:
            setattr(self, key.replace("/", "_"), reader.headerProperties[key])
        self.nbpForDSB = int(self.Damagedefinition[2])
        for d in reader.damages:
            dsite = SDDDamageSite(version)
            for key in d:
                setattr(dsite, key, d[key])
            dsite.initialBp = np.round(dsite.ChromosomePosition * self.Chromosomesizes[dsite.chromosomeNumber]*1e6)
            dsite.ParticleTime = particleTime
            dsite.LesionTime = lesionTime
            self.damageSites.append(dsite)
            if int(dsite.newExposure) > 0:
                self.numtracks += 1
                self.damagespertrack['DSB'].append(0)
                self.damagespertrack['SB'].append(0)
                self.damagespertrack['BD'].append(0)
                self.damagespertrack['NumSites'].append(0)
            self.damagespertrack['DSB'][self.numtracks-1] += dsite.numberofDSBs
            self.damagespertrack['SB'][self.numtracks-1] += dsite.numberofstrandbreaks
            self.damagespertrack['BD'][self.numtracks-1] += dsite.numberofbasedamages
            self.damagespertrack['NumSites'][self.numtracks-1] += 1

    def recomputeDamagesFromReadSites(self, stopAtDose = -1, stopAtTime = -1, recalculatePerEachTrack=False):
        self.recomputeSitesAtCurrentTime()
        self.populateDamages(stopAtDose, stopAtTime, recalculatePerEachTrack)
        self.computeStrandBreaks()

    def populateDamages(self, stopAtDose = -1, stopAtTime = -1, recalculatePerEachTrack=False, recalculateEveryQuarterOfGray=False,
                        classifySites=True, getVideo=False):
        # This will allow to read the SDD file and populate the damage sites
        # Initial dose to stop
        if stopAtDose >= 0:
            self.stopAtDose = stopAtDose
        iExposure = -1
        iEvent = 0
        # Handling a parallel list of damage sites considering time
        damageSites = self.damageSites.copy()
        if stopAtTime == -1:
            self.damageSitesAtCurrentTime = damageSites.copy()
        else:
            for d in damageSites:
                if d.LesionTime <= stopAtTime:
                    self.damageSitesAtCurrentTime.insert(0, d)
        del(damageSites)
        # Initialize the damage map
        dosefromlastrecalculation = 0
        for id, damage in enumerate(self.damageSitesAtCurrentTime):
            if (0 <= self.stopAtDose < self.cumulativeDose) or damage.isRepaired:
                self.computeStrandBreaks()
                break
            iCh = damage.chromosomeNumber
            if iCh not in self.damageMap.keys():
                self.damageMap[iCh] = {}
            for bpdamage in damage.individualdamages:
                iBp = int(damage.initialBp+bpdamage['basepairID'])
                if iBp not in self.damageMap[iCh].keys():
                    self.damageMap[iCh][iBp] = {}
                extension = np.abs(np.array([damage.maxX - damage.minX, damage.maxY - damage.minY, damage.maxZ - damage.minZ]))
                position = np.array([damage.centerX, damage.centerY, damage.centerZ]) + (-5+int(bpdamage['basepairID']))/10 * extension + (-2.5+bpdamage['subcomponent']) * np.array([1e-3, 0, 0])
                self.damageMap[iCh][damage.initialBp+bpdamage['basepairID']][bpdamage['subcomponent']] = \
                    SubcomponentLesion(bpdamage['type'], position, damage.LesionTime, damage.ParticleTime, iEvent)
            if (id < len(self.damageSitesAtCurrentTime) - 1 and self.damageSitesAtCurrentTime[id + 1].newExposure > 1) or id == len(self.damageSitesAtCurrentTime) - 1:
                if id < len(self.damageSitesAtCurrentTime) - 1 and self.damageSitesAtCurrentTime[id + 1].newExposure > 1:
                    iExposure += 1
                    iEvent += 1
                if recalculatePerEachTrack or id == len(self.damageSitesAtCurrentTime) - 1:
                    self.computeStrandBreaks(classifySites=classifySites)
                    dosefromlastrecalculation = 0
                if len(self.doses) > 0:
                    self.cumulativeDose += self.doses[iExposure]
                    dosefromlastrecalculation += self.doses[iExposure]
                    if recalculateEveryQuarterOfGray and dosefromlastrecalculation >= 0.25:
                        self.computeStrandBreaks(classifySites=classifySites)
                        dosefromlastrecalculation = 0
                    self.Darray.append(self.cumulativeDose)
                    self.DSBarray.append(self.numDSB); self.DSBdirectarray.append(self.numDSBDirect); self.DSBindirectarray.append(self.numDSBIndirect); self.DSBhybridarray.append(self.numDSBHybrid)
                    self.SSBarray.append(self.numSSB); self.SSBdirectarray.append(self.numSSBDirect); self.SSBindirectarray.append(self.numSSBIndirect)
                    self.SBarray.append(self.numSB); self.SBdirectarray.append(self.numSBDirect); self.SBindirectarray.append(self.numSBIndirect)
                    self.BDarray.append(self.numBD); self.BDdirectarray.append(self.numBDDirect); self.BDindirectarray.append(self.numBDIndirect)
                    self.numSitesarray.append(self.provnsites)
                    if getVideo:
                        self.produce3DImage(show=False)
            if damage.newExposure == 1:
                iEvent += 1
            self.damageSites.remove(damage)  # Removing to avoid adding the same damage again

    def computeStrandBreaks(self, classifySites=True):
        self.DSBMap = {iCh: {iBp: {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)} for iBp in
                  self.damageMap[iCh]} for iCh in self.damageMap}
        DSBPairs = {iCh: {} for iCh in self.damageMap}
        self.SSBMap = {}
        self.BDMap = {}
        self.DSBPositions = []
        numberofdsbs = 0
        for iCh in self.damageMap:
            bp_keys = {iBp for iBp in self.damageMap[iCh]}
            for iBp in bp_keys:
                # Checks if there is damage in backbone 1 (iComp = 2) and there is not already a DSB identified in strand 1
                if 2 in self.damageMap[iCh][iBp] and self.damageMap[iCh][iBp][2].type > 0 and (1 not in self.DSBMap[iCh][iBp] or self.DSBMap[iCh][iBp][1].type == 0):
                    dsbFound = False
                    for i2 in range(0, self.nbpForDSB):
                        if iBp+i2 in self.damageMap[iCh] and 3 in self.damageMap[iCh][iBp+i2]:
                            typeDamageInStrand2 = self.damageMap[iCh][iBp+i2][3].type
                            if typeDamageInStrand2 > 0:
                                if iBp+i2 not in self.DSBMap[iCh]:
                                    self.DSBMap[iCh][iBp+i2] = {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                if 2 not in self.DSBMap[iCh][iBp + i2] or self.DSBMap[iCh][iBp + i2][2].type == 0:
                                    adjustedTypeDamageInStrand2 = 2
                                    # Checks if damage in backbone 2 is direct and keeps it in that case
                                    # If it is multiple damage, considers direct (since it would not happen without chemistry)
                                    if typeDamageInStrand2 == 1 or typeDamageInStrand2 == 3:
                                        adjustedTypeDamageInStrand2 = 1

                                    # Once damage is identified in strand 2 (and damage ensured in strand 1), looks for the closest damage in strand 1
                                    closestPosInStrand1 = iBp
                                    adjustedTypeDamageInStrand1 = 2
                                    for i1 in range(0, self.nbpForDSB):
                                        if iBp + i2 + i1 in bp_keys and 2 in self.damageMap[iCh][iBp + i2 + i1]:
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2+i1][2].type
                                            if typeDamageInStrand1 > 0:
                                                if iBp+i2+i1 not in self.DSBMap[iCh]:
                                                    self.DSBMap[iCh][iBp + i2 + i1] = {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                                if 1 not in self.DSBMap[iCh][iBp+i2+i1] or self.DSBMap[iCh][iBp+i2+i1][1].type == 0:
                                                    closestPosInStrand1 = iBp + i2 + i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp+i2-i1 in bp_keys and 2 in self.damageMap[iCh][iBp+i2-i1]:
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2-i1][2].type
                                            if iBp+i2-i1 >= 0 and typeDamageInStrand1 > 0:
                                                if iBp+i2-i1 not in self.DSBMap[iCh]:
                                                    self.DSBMap[iCh][iBp + i2 - i1] =  {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                                if 1 not in self.DSBMap[iCh][iBp+i2-i1] or self.DSBMap[iCh][iBp+i2-i1][1].type == 0:
                                                    closestPosInStrand1 = iBp + i2 - i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                if adjustedTypeDamageInStrand1 > 0:
                                    posSB1 = np.array(self.damageMap[iCh][closestPosInStrand1][2].position)
                                    posSB2 = np.array(self.damageMap[iCh][iBp + i2][3].position)
                                    self.DSBMap[iCh][closestPosInStrand1][1] = SubcomponentLesion(adjustedTypeDamageInStrand1, posSB1, 0,
                                                                                                  particletime=self.damageMap[iCh][closestPosInStrand1][2].particletime, dsbid=numberofdsbs)
                                    self.DSBMap[iCh][iBp + i2][2] = SubcomponentLesion(adjustedTypeDamageInStrand2, posSB2, 0,
                                                                                                  particletime=self.damageMap[iCh][iBp + i2][3].particletime, dsbid=numberofdsbs)
                                    numberofdsbs += 1
                                    pos = (closestPosInStrand1, iBp + i2)
                                    DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                    self.DSBPositions.append((posSB1 + posSB2)/2)
                                    dsbFound = True
                        if dsbFound:
                            break
                        if iBp-i2 >= 0 and iBp-i2 in self.damageMap[iCh].keys() and 3 in self.damageMap[iCh][iBp-i2].keys():
                            typeDamageInStrand2 = self.damageMap[iCh][iBp - i2][3].type
                            if i2 > 0 and typeDamageInStrand2 > 0:
                                if iBp-i2 not in self.DSBMap[iCh]:
                                    self.DSBMap[iCh][iBp-i2] = {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                if 2 not in self.DSBMap[iCh][iBp-i2] or self.DSBMap[iCh][iBp-i2][2].type == 0:
                                    adjustedTypeDamageInStrand2 = 2
                                    if typeDamageInStrand2 == 1 or typeDamageInStrand2 == 3:
                                        adjustedTypeDamageInStrand2 = 1
                                    closestPosInStrand1 = iBp
                                    adjustedTypeDamageInStrand1 = 2
                                    for i1 in range(0, self.nbpForDSB):
                                        if iBp-i2+i1 in self.damageMap[iCh] and 2 in self.damageMap[iCh][iBp-i2+i1]:
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2+i1][2].type
                                            if typeDamageInStrand1 > 0:
                                                if iBp-i2+i1 not in self.DSBMap[iCh]:
                                                    self.DSBMap[iCh][iBp - i2 + i1] = {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                                if 1 not in self.DSBMap[iCh][iBp-i2+i1] or self.DSBMap[iCh][iBp-i2+i1][1].type == 0:
                                                    closestPosInStrand1 = iBp-i2+i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp-i2-i1 in bp_keys and 2 in self.damageMap[iCh][iBp - i2 - i1]:
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2-i1][2].type
                                            if iBp-i2-i1 >= 0 and i1 > 0 and typeDamageInStrand1 > 0:
                                                if iBp-i2-i1 not in self.DSBMap[iCh]:
                                                    self.DSBMap[iCh][iBp - i2 - i1] = {1: SubcomponentLesion(0, None, 0, 0), 2: SubcomponentLesion(0, None, 0, 0)}
                                                if 1 not in self.DSBMap[iCh][iBp-i2-i1].keys() or self.DSBMap[iCh][iBp-i2-i1][1].type == 0:
                                                    closestPosInStrand1 = iBp-i2-i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                    if adjustedTypeDamageInStrand1 > 0:
                                        posSB1 = np.array(self.damageMap[iCh][closestPosInStrand1][2].position)
                                        posSB2 = np.array(self.damageMap[iCh][iBp - i2][3].position)
                                        self.DSBMap[iCh][closestPosInStrand1][1] = SubcomponentLesion(adjustedTypeDamageInStrand1, posSB1, 0,
                                                                                                      particletime=self.damageMap[iCh][closestPosInStrand1][2].particletime, dsbid=numberofdsbs)
                                        self.DSBMap[iCh][iBp - i2][2] = SubcomponentLesion(adjustedTypeDamageInStrand2, posSB2, 0,
                                                                                                      particletime=self.damageMap[iCh][iBp - i2][3].particletime, dsbid=numberofdsbs)
                                        numberofdsbs += 1
                                        pos = (closestPosInStrand1, iBp - i2)
                                        DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                        self.DSBPositions.append((posSB1 + posSB2)/2)
                                        dsbFound = True
                        if dsbFound:
                            break

        #Loops through damage map to exclude DSBs
        for iCh in self.damageMap.keys():
            if iCh not in self.SSBMap.keys():
                self.SSBMap[iCh] = {}
            if iCh not in self.BDMap.keys():
                self.BDMap[iCh] = {}
            for iBp in self.damageMap[iCh]:
                if iBp not in self.SSBMap[iCh].keys():
                    self.SSBMap[iCh][iBp] = {}
                if iBp not in self.BDMap[iCh].keys():
                    self.BDMap[iCh][iBp] = {}
                # Single Strand breaks
                if (2 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][2].type > 0 and (1 not in self.DSBMap[iCh][iBp].keys() or
                                                                                                  not self.DSBMap[iCh][iBp][1].type > 0)) or (3 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][3].type > 0 and
                                                                                                (2 not in self.DSBMap[iCh][iBp].keys() or not self.DSBMap[iCh][iBp][2].type > 0)):
                    if 2 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][2].type > 0:
                        self.SSBMap[iCh][iBp][1] = SubcomponentLesion(self.damageMap[iCh][iBp][2].type, self.damageMap[iCh][iBp][2].position, self.damageMap[iCh][iBp][2].lesiontime, self.damageMap[iCh][iBp][2].particletime)
                    if 3 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][3].type > 0:
                        self.SSBMap[iCh][iBp][2] = SubcomponentLesion(self.damageMap[iCh][iBp][3].type, self.damageMap[iCh][iBp][3].position, self.damageMap[iCh][iBp][3].lesiontime, self.damageMap[iCh][iBp][3].particletime)
                # Base damages
                if 1 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][1].type > 0:
                    self.BDMap[iCh][iBp][1] = SubcomponentLesion(self.damageMap[iCh][iBp][1].type, self.damageMap[iCh][iBp][1].position, self.damageMap[iCh][iBp][1].lesiontime, self.damageMap[iCh][iBp][1].particletime)
                if 4 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][4].type > 0:
                    self.BDMap[iCh][iBp][2] = SubcomponentLesion(self.damageMap[iCh][iBp][4].type, self.damageMap[iCh][iBp][4].position, self.damageMap[iCh][iBp][4].lesiontime, self.damageMap[iCh][iBp][4].particletime)
        # Remove all entries with no damage
        # DSB Map
        keys_to_remove = []
        for iCh in self.DSBMap.keys():
            for iBp in self.DSBMap[iCh].keys():
                if self.DSBMap[iCh][iBp][1].type == 0 and self.DSBMap[iCh][iBp][2].type == 0:
                    keys_to_remove.append((iCh, iBp))
        for key in keys_to_remove:
            del self.DSBMap[key[0]][key[1]]
        keys_to_remove = []
        for iCh in self.DSBMap.keys():
            if len(self.DSBMap[iCh]) == 0:
                keys_to_remove.append(iCh)
        for key in keys_to_remove:
            del self.DSBMap[key]

        # SSB Map
        keys_to_remove = []
        for iCh in self.SSBMap.keys():
            for iBp in self.SSBMap[iCh].keys():
                if 1 in self.SSBMap[iCh][iBp].keys() and self.SSBMap[iCh][iBp][1].type == 0:
                    keys_to_remove.append((iCh, iBp, 1))
                if 2 in self.SSBMap[iCh][iBp].keys() and self.SSBMap[iCh][iBp][2].type == 0:
                    keys_to_remove.append((iCh, iBp, 2))
        for key in keys_to_remove:
            del self.SSBMap[key[0]][key[1]][key[2]]
        keys_to_remove = []
        for iCh in self.SSBMap.keys():
            for iBp in self.SSBMap[iCh].keys():
                if len(self.SSBMap[iCh][iBp]) == 0:
                    keys_to_remove.append((iCh, iBp))
        for key in keys_to_remove:
            del self.SSBMap[key[0]][key[1]]
        keys_to_remove = []
        for iCh in self.SSBMap.keys():
            if len(self.SSBMap[iCh]) == 0:
                keys_to_remove.append(iCh)
        for key in keys_to_remove:
            del self.SSBMap[key]

        # BD Map
        keys_to_remove = []
        for iCh in self.BDMap.keys():
            for iBp in self.BDMap[iCh].keys():
                if 1 in self.BDMap[iCh][iBp].keys() and self.BDMap[iCh][iBp][1].type == 0:
                    keys_to_remove.append((iCh, iBp, 1))
                if 2 in self.BDMap[iCh][iBp].keys() and self.BDMap[iCh][iBp][2].type == 0:
                    keys_to_remove.append((iCh, iBp, 2))
        for key in keys_to_remove:
            del self.BDMap[key[0]][key[1]][key[2]]
        keys_to_remove = []
        for iCh in self.BDMap.keys():
            for iBp in self.BDMap[iCh].keys():
                if len(self.BDMap[iCh][iBp]) == 0:
                    keys_to_remove.append((iCh, iBp))
        for key in keys_to_remove:
            del self.BDMap[key[0]][key[1]]
        keys_to_remove = []
        for iCh in self.BDMap.keys():
            if len(self.BDMap[iCh]) == 0:
                keys_to_remove.append(iCh)
        for key in keys_to_remove:
            del self.BDMap[key]

        # DSB pairs
        keys_to_remove = []
        for iCh in DSBPairs.keys():
            for iBp in DSBPairs[iCh].keys():
                if DSBPairs[iCh][iBp] == 0:
                    keys_to_remove.append((iCh, iBp))
        for key in keys_to_remove:
            del DSBPairs[key[0]][key[1]]
        keys_to_remove = []
        for iCh in DSBPairs.keys():
            if DSBPairs[iCh] == {}:
                keys_to_remove.append(iCh)
        for key in keys_to_remove:
            del DSBPairs[key]
        self.quantifyDamage(DSBPairs, self.SSBMap, self.BDMap)
        if classifySites:
            self.classifiedDamage = self.classifyDamageSites()

    def quantifyDamage(self, DSBPairs, SSBMap, BDMap):
        self.initializeCounters()
        for iCh in DSBPairs:
            for pos in DSBPairs[iCh]:
                typeDamage = DSBPairs[iCh][pos]
                if typeDamage > 0:
                    self.numDSB += 1
                    self.numSB += 2
                    if typeDamage == 2:
                        self.numDSBDirect += 1; self.numSBDirect += 2
                    if typeDamage == 3:
                        self.numDSBHybrid += 1; self.numSBDirect += 1; self.numSBIndirect += 1
                    if typeDamage == 4:
                        self.numDSBIndirect += 1; self.numSBIndirect += 2
        for iCh in SSBMap:
            for iBp in SSBMap[iCh]:
                for iCo in SSBMap[iCh][iBp]:
                    typeDamage = SSBMap[iCh][iBp][iCo].type
                    if typeDamage > 0:
                        self.numSSB += 1
                        self.numSB += 1
                        if typeDamage == 1 or typeDamage == 3:
                            self.numSSBDirect += 1; self.numSBDirect += 1
                        if typeDamage == 2:
                            self.numSSBIndirect += 1; self.numSBIndirect += 1
        for iCh in BDMap:
            for iBp in BDMap[iCh]:
                for iCo in BDMap[iCh][iBp]:
                    typeDamage = BDMap[iCh][iBp][iCo].type
                    if typeDamage > 0:
                        self.numBD += 1
                        if typeDamage == 1 or typeDamage == 3:
                            self.numBDDirect += 1
                        if typeDamage == 2:
                            self.numBDIndirect += 1

    def classifyDamageSites(self):
        damageScore = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        for iCh, bdMap in self.BDMap.items():
            for iBp, bdCoMap in bdMap.items():
                for iCo, bd in bdCoMap.items():
                    if bd.type > 0:
                        damageScore[iCh][iBp][iCo] += 0.001
        for iCh, ssbMap in self.SSBMap.items():
            for iBp, ssbCoMap in ssbMap.items():
                for iCo, ssb in ssbCoMap.items():
                    if ssb.type > 0:
                        damageScore[iCh][iBp][iCo] += 0.1
        for iCh, dsbMap in self.DSBMap.items():
            for iBp, dsbCoMap in dsbMap.items():
                for iCo, dsb in dsbCoMap.items():
                    if dsb.type > 0:
                        damageScore[iCh][iBp][iCo] += 5.0

        chromosomeAndInitialBpIdOfDamageSites = {}
        maxDamage = 10 * self.nbpForDSB
        minDamage = 1e-8
        # Optimize
        for iCh, damageScoreCh in damageScore.items():
            startingBpIdForDamageSites = []
            basePairs = list(damageScoreCh.keys())
            while maxDamage > minDamage and len(basePairs) != 0:
                linearAccumulatedDamage = (sum(damageScoreCh[bp + i][1] + damageScoreCh[bp + i][2] for i in range(self.nbpForDSB))
                    for bp in basePairs)
                initialBpIds = (bp for bp in basePairs)
                maxDamage, index = max((d, i) for d, i in zip(linearAccumulatedDamage, initialBpIds))
                startingBpIdForDamageSites.append(index)
                self.provnsites += 1
                for i in range(self.nbpForDSB):
                    self.assignComplexities(iCh, index + i, maxDamage)
                basePairs = [bp for bp in basePairs if bp < index or bp > index + self.nbpForDSB]
            chromosomeAndInitialBpIdOfDamageSites[iCh] = startingBpIdForDamageSites
        return chromosomeAndInitialBpIdOfDamageSites

    def assignComplexities(self, iChr, iBp, damage):
        for i in range(self.nbpForDSB):
            if iChr in self.BDMap.keys():
                if iBp + i in self.BDMap[iChr].keys():
                    if 1 in self.BDMap[iChr][iBp + i].keys():
                        self.BDMap[iChr][iBp+i][1].complexity = damage
                    if 2 in self.BDMap[iChr][iBp + i].keys():
                        self.BDMap[iChr][iBp+i][2].complexity = damage
            if iChr in self.SSBMap.keys():
                if iBp + i in self.SSBMap[iChr].keys():
                    if 1 in self.SSBMap[iChr][iBp + i].keys():
                        self.SSBMap[iChr][iBp+i][1].complexity = damage
                    if 2 in self.SSBMap[iChr][iBp + i].keys():
                        self.SSBMap[iChr][iBp+i][2].complexity = damage
            if iChr in self.DSBMap.keys():
                if iBp + i in self.DSBMap[iChr].keys():
                    if 1 in self.DSBMap[iChr][iBp + i].keys():
                        self.DSBMap[iChr][iBp+i][1].complexity = damage
                    if 2 in self.DSBMap[iChr][iBp + i].keys():
                        self.DSBMap[iChr][iBp+i][2].complexity = damage

    def recomputeSitesAtCurrentTime(self):
        sitesatcurrentime = self.classifiedDamage
        numSites = 0
        firstExposure = True
        newEvent = True
        self.damageSitesAtCurrentTime.clear()
        for iCh in sitesatcurrentime.keys():
            ibpsTakenForThisChromosome = []
            for i in range(len(sitesatcurrentime[iCh])):
                newsite = SDDDamageSite()
                newsite.isRepaired = True
                initialBpId = sitesatcurrentime[iCh][i]
                dir = 0; indir = 0
                bd = 0; sb = 0; dsb = 0
                for j in range(self.nbpForDSB):
                    if initialBpId + j not in ibpsTakenForThisChromosome:
                        if iCh in self.SSBMap.keys():
                            if initialBpId + j in self.SSBMap[iCh].keys():
                                if 1 in self.SSBMap[iCh][initialBpId + j].keys():
                                    if self.SSBMap[iCh][initialBpId + j][1] == 1 or self.SSBMap[iCh][initialBpId + j][1].type == 3:
                                        dir += 1
                                        sb += 1
                                    if self.SSBMap[iCh][initialBpId + j][1] == 2:
                                        indir += 1
                                        sb += 1
                                if 2 in self.SSBMap[iCh][initialBpId + j].keys():
                                    if self.SSBMap[iCh][initialBpId + j][2] == 1 or self.SSBMap[iCh][initialBpId + j][2].type == 3:
                                        dir += 1
                                        sb += 1
                                    if self.SSBMap[iCh][initialBpId + j][2] == 2:
                                        indir += 1
                                        sb += 1
                        if iCh in self.DSBMap.keys():
                            if initialBpId + j in self.DSBMap[iCh].keys():
                                if 1 in self.DSBMap[iCh][initialBpId + j].keys():
                                    if self.DSBMap[iCh][initialBpId + j][1].type == 1 or self.DSBMap[iCh][initialBpId + j][1].type == 3:
                                        dir += 1
                                        sb += 1
                                        dsb += 1
                                    if self.DSBMap[iCh][initialBpId + j][1].type == 2:
                                        indir += 1
                                        sb += 1
                                        dsb += 1
                                if 2 in self.DSBMap[iCh][initialBpId + j].keys():
                                    if self.DSBMap[iCh][initialBpId + j][2].type == 1 or self.DSBMap[iCh][initialBpId + j][2].type == 3:
                                        dir += 1
                                        sb += 1
                                    if self.DSBMap[iCh][initialBpId + j][2].type == 2:
                                        indir += 1
                                        sb += 1
                        if iCh in self.BDMap.keys():
                            if initialBpId + j in self.BDMap[iCh].keys():
                                if 1 in self.BDMap[iCh][initialBpId + j].keys():
                                    if self.BDMap[iCh][initialBpId + j][1].type == 1 or self.BDMap[iCh][initialBpId + j][1].type == 3:
                                        dir += 1
                                        bd += 1
                                    if self.BDMap[iCh][initialBpId + j][1].type == 2:
                                        indir += 1
                                        bd += 1
                                if 2 in self.BDMap[iCh][initialBpId + j].keys():
                                    if self.BDMap[iCh][initialBpId + j][2].type == 1 or self.BDMap[iCh][initialBpId + j][2].type == 3:
                                        dir += 1
                                        bd += 1
                                    if self.BDMap[iCh][initialBpId + j][2].type == 2:
                                        indir += 1
                                        bd += 1
                    ibpsTakenForThisChromosome.append(initialBpId + j)
                numSites += 1
                #Field 1
                newExposureFlag = str(0)
                if newEvent:
                    newExposureFlag = str(1)
                    newEvent = False
                if firstExposure:
                    newExposureFlag = str(2)
                    firstExposure = False
                newsite.Classification = [newExposureFlag, 0]
                #Field 2
                damagePositions = []
                for j in range(self.nbpForDSB):
                    if iCh in self.damageMap.keys():
                        if initialBpId + j in self.damageMap[iCh].keys():
                            if 1 in self.damageMap[iCh][initialBpId + j].keys():
                                damagePositions.append(self.damageMap[iCh][initialBpId + j][1].position)
                            if 2 in self.damageMap[iCh][initialBpId + j].keys():
                                damagePositions.append(self.damageMap[iCh][initialBpId + j][2].position)
                            if 3 in self.damageMap[iCh][initialBpId + j].keys():
                                damagePositions.append(self.damageMap[iCh][initialBpId + j][3].position)
                            if 4 in self.damageMap[iCh][initialBpId + j].keys():
                                damagePositions.append(self.damageMap[iCh][initialBpId + j][4].position)
                if len(damagePositions) > 0:
                    damageCenterMaxMin = self.getDamageCenterAndBoundaries(np.array(damagePositions))
                    center = damageCenterMaxMin[0]
                    max = damageCenterMaxMin[1]
                    min = damageCenterMaxMin[2]
                    newsite.SpatialCoordinatesAndExtent = [center[0], center[1], center[2], max[0], max[1], max[2], min[0], min[1], min[2]]
                # Field 3
                newsite.ChromosomeID = [1, iCh, 1, 0]
                # Field 4
                chromosomeLength = self.Chromosomesizes[iCh]
                damageChromPos = initialBpId / chromosomeLength / 1e6
                newsite.ChromosomePosition = float("{:.12f}".format(damageChromPos))
                # Field 5
                typeDamage = 0
                if dir == 0 and indir > 0:
                    typeDamage = 1
                if dir > 0 and indir > 0:
                    typeDamage = 2
                newsite.Cause = [typeDamage, dir, indir]
                # Field 6
                newsite.DamageTypes = [bd, sb, dsb]
                # Field 7
                damageSpec = []
                for j in range(0, self.nbpForDSB):
                    if iCh in self.damageMap.keys():
                        if initialBpId + j in self.damageMap[iCh].keys():
                            if 2 in self.damageMap[iCh][initialBpId + j].keys():
                                typeDamage = self.damageMap[iCh][initialBpId + j][2].type
                                if typeDamage <= 0:
                                    typeDamage = 0
                                damageSpec.append(1)
                                damageSpec.append(j+1)
                                damageSpec.append(typeDamage)
                for j in range(0, self.nbpForDSB):
                    if iCh in self.damageMap.keys():
                        if initialBpId + j in self.damageMap[iCh].keys():
                            if 1 in self.damageMap[iCh][initialBpId + j].keys():
                                typeDamage = self.damageMap[iCh][initialBpId + j][1].type
                                if typeDamage <= 0:
                                    typeDamage = 0
                                damageSpec.append(2)
                                damageSpec.append(j+1)
                                damageSpec.append(typeDamage)
                for j in range(0, self.nbpForDSB):
                    if iCh in self.damageMap.keys():
                        if initialBpId + j in self.damageMap[iCh].keys():
                            if 4 in self.damageMap[iCh][initialBpId + j].keys():
                                typeDamage = self.damageMap[iCh][initialBpId + j][4].type
                                if typeDamage <= 0:
                                    typeDamage = 0
                                damageSpec.append(3)
                                damageSpec.append(j+1)
                                damageSpec.append(typeDamage)
                for j in range(0, self.nbpForDSB):
                    if iCh in self.damageMap.keys():
                        if initialBpId + j in self.damageMap[iCh].keys():
                            if 3 in self.damageMap[iCh][initialBpId + j].keys():
                                typeDamage = self.damageMap[iCh][initialBpId + j][3].type
                                if typeDamage <= 0:
                                    typeDamage = 0
                                damageSpec.append(4)
                                damageSpec.append(j+1)
                                damageSpec.append(typeDamage)
                newsite.FullBreakSpec = damageSpec
                self.damageSitesAtCurrentTime.append(newsite)

    def writeSDD(self, filename, version='2.0'):
        self.outputSDDHeader(filename)
        self.outputSDDFile(filename, version)

    def outputSDDHeader(self, filename):
        dataEntries = ""
        for i, e in enumerate(self.Dataentries):
            dataEntries = dataEntries + str(int(e))
            if i < len(self.Dataentries) - 1:
                dataEntries += ", "
        with open(filename, 'w') as f:
            f.write("SDD Version, " + str(self.SDDVersion) + ';\n')
            f.write("Software, " + str(self.Software) + ';\n')
            f.write("Author, " + str(self.Author) + ';\n')
            f.write("Simulation Details, " + str(self.SimulationDetails) + ';\n')
            f.write("Source, " + str(self.Source) + ';\n')
            f.write("Source type, " + str(int(self.Sourcetype)) + ';\n')
            if self.Incidentparticles != '':
                f.write("Incident particles, " + str(int(self.Incidentparticles))  + ';\n')
            else:
                f.write("Incident particles, " + str(self.Incidentparticles) + ';\n')
            f.write("Mean particle energy, " + str(self.Meanparticleenergy)  + ';\n')
            f.write("Energy distribution, " + self.getStringOutOfList(self.Energydistribution) + ';\n')
            f.write("Particle fraction, " + str(self.Particlefraction) + ';\n')
            f.write("Dose or fluence, " + self.getStringOutOfList(self.Doseorfluence) + ';\n')
            f.write("Dose rate, " + str(0.0) + ';\n')
            f.write("Irradiation target, " + str(self.Irradiationtarget) + ';\n')
            f.write("Volumes, " + self.getStringOutOfList(self.Volumes) + ';\n')
            f.write("Chromosome sizes, " + self.getStringOutOfList(self.Chromosomesizes) + ';\n')
            f.write("DNA Density, " + str(self.DNADensity) + ';\n')
            f.write("Cell Cycle Phase, " + str(self.CellCyclePhase) + ';\n')
            f.write("DNA Structure, " + self.getStringOutOfList(self.DNAStructure) + ';\n')
            f.write("In vitro / in vivo, " + str(self.Invitro_invivo) + ';\n')
            f.write("Proliferation status, " + str(int(self.Proliferationstatus)) + ';\n')
            f.write("Microenvironment, " + self.getStringOutOfList(self.Microenvironment) + ';\n')
            f.write("Damage definition, " + self.getStringOutOfList(self.Damagedefinition) + ';\n')
            f.write("Time, " + str(self.Time) + ';\n')
            f.write("Damage and primary count, " + str(self.Damageandprimarycount) + ';\n')
            f.write("Data entries, " + self.getStringOutOfList(self.Dataentries, 'int') + ';\n')
            f.write("Additional information, " + str(self.Additionalinformation) + ';\n')
            f.write("***EndOfHeader***;\n")

    def outputSDDFile(self, filename, version='2.0'):
        damageSitesPerTime = self.classifiedDamageSites
        numSites = 0
        firstExposure = True
        for damageSites in damageSitesPerTime:
            newEvent = True
            with open(filename, 'a') as f:
                for iCh in damageSites.keys():
                    ibpsTakenForThisChromosome = []
                    for i in range(len(damageSites[iCh])):
                        initialBpId = damageSites[iCh][i]
                        dir = 0; indir = 0
                        bd = 0; sb = 0; dsb = 0
                        for j in range(self.nbpForDSB):
                            if initialBpId + j not in ibpsTakenForThisChromosome:
                                if iCh in self.SSBMap.keys():
                                    if initialBpId + j in self.SSBMap[iCh].keys():
                                        if 1 in self.SSBMap[iCh][initialBpId + j].keys():
                                            if self.SSBMap[iCh][initialBpId + j][1] == 1 or self.SSBMap[iCh][initialBpId + j][1].type == 3:
                                                dir += 1
                                                sb += 1
                                            if self.SSBMap[iCh][initialBpId + j][1] == 2:
                                                indir += 1
                                                sb += 1
                                        if 2 in self.SSBMap[iCh][initialBpId + j].keys():
                                            if self.SSBMap[iCh][initialBpId + j][2] == 1 or self.SSBMap[iCh][initialBpId + j][2].type == 3:
                                                dir += 1
                                                sb += 1
                                            if self.SSBMap[iCh][initialBpId + j][2] == 2:
                                                indir += 1
                                                sb += 1
                                if iCh in self.DSBMap.keys():
                                    if initialBpId + j in self.DSBMap[iCh].keys():
                                        if 1 in self.DSBMap[iCh][initialBpId + j].keys():
                                            if self.DSBMap[iCh][initialBpId + j][1].type == 1 or self.DSBMap[iCh][initialBpId + j][1].type == 3:
                                                dir += 1
                                                sb += 1
                                                dsb += 1
                                            if self.DSBMap[iCh][initialBpId + j][1].type == 2:
                                                indir += 1
                                                sb += 1
                                                dsb += 1
                                        if 2 in self.DSBMap[iCh][initialBpId + j].keys():
                                            if self.DSBMap[iCh][initialBpId + j][2].type == 1 or self.DSBMap[iCh][initialBpId + j][2].type == 3:
                                                dir += 1
                                                sb += 1
                                            if self.DSBMap[iCh][initialBpId + j][2].type == 2:
                                                indir += 1
                                                sb += 1
                                if iCh in self.BDMap.keys():
                                    if initialBpId + j in self.BDMap[iCh].keys():
                                        if 1 in self.BDMap[iCh][initialBpId + j].keys():
                                            if self.BDMap[iCh][initialBpId + j][1].type == 1 or self.BDMap[iCh][initialBpId + j][1].type == 3:
                                                dir += 1
                                                bd += 1
                                            if self.BDMap[iCh][initialBpId + j][1].type == 2:
                                                indir += 1
                                                bd += 1
                                        if 2 in self.BDMap[iCh][initialBpId + j].keys():
                                            if self.BDMap[iCh][initialBpId + j][2].type == 1 or self.BDMap[iCh][initialBpId + j][2].type == 3:
                                                dir += 1
                                                bd += 1
                                            if self.BDMap[iCh][initialBpId + j][2].type == 2:
                                                indir += 1
                                                bd += 1
                            ibpsTakenForThisChromosome.append(initialBpId + j)
                        numSites += 1
                        #Field 1
                        newExposureFlag = str(0)
                        if newEvent:
                            newExposureFlag = str(1)
                            newEvent = False
                        if firstExposure:
                            newExposureFlag = str(2)
                            firstExposure = False
                        f.write(newExposureFlag + ', 0; ')
                        #Field 2
                        damagePositions = []
                        for j in range(self.nbpForDSB):
                            if iCh in self.damageMap.keys():
                                if initialBpId + j in self.damageMap[iCh].keys():
                                    if 1 in self.damageMap[iCh][initialBpId + j].keys():
                                        damagePositions.append(self.damageMap[iCh][initialBpId + j][1].position)
                                    if 2 in self.damageMap[iCh][initialBpId + j].keys():
                                        damagePositions.append(self.damageMap[iCh][initialBpId + j][2].position)
                                    if 3 in self.damageMap[iCh][initialBpId + j].keys():
                                        damagePositions.append(self.damageMap[iCh][initialBpId + j][3].position)
                                    if 4 in self.damageMap[iCh][initialBpId + j].keys():
                                        damagePositions.append(self.damageMap[iCh][initialBpId + j][4].position)
                        if len(damagePositions) > 0:
                            damageCenterMaxMin = self.getDamageCenterAndBoundaries(np.array(damagePositions))
                            center = damageCenterMaxMin[0]
                            f.write(str(center[0]) + ", " + str(center[1]) + ", " + str(center[2]))
                            max = damageCenterMaxMin[1]
                            f.write(" / " + str(max[0]) + ", " + str(max[1]) + ", " + str(max[2]))
                            min = damageCenterMaxMin[2]
                            f.write(" / " + str(min[0]) + ", " + str(min[1]) + ", " + str(min[2]) + "; ")
                        # Field 3
                        f.write("1, " + str(iCh) + ", 1, 0; ")
                        # Field 4
                        chromosomeLength = self.Chromosomesizes[iCh]
                        damageChromPos = initialBpId / chromosomeLength / 1e6
                        f.write("{:.12f}".format(damageChromPos) + "; ")
                        # Field 5
                        typeDamage = 0
                        if dir == 0 and indir > 0:
                            typeDamage = 1
                        if dir > 0 and indir > 0:
                            typeDamage = 2
                        f.write(str(typeDamage) + ", " + str(dir) + ", " + str(indir) + "; ")
                        # Field 6
                        f.write(str(bd) + ", " + str(sb) + ", " + str(dsb) + "; ")
                        # Field 7
                        damageSpec = ''
                        if float(version) < 2.0:
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 1 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][1].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "1, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 2 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][2].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "2, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 3 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][3].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "3, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 4 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][4].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "4, " + str(j+1) + ", " + str(typeDamage) + " / "
                        else:
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 2 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][2].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "1, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 1 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][1].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "2, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 4 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][4].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "3, " + str(j+1) + ", " + str(typeDamage) + " / "
                            for j in range(0, self.nbpForDSB):
                                if iCh in self.damageMap.keys():
                                    if initialBpId + j in self.damageMap[iCh].keys():
                                        if 3 in self.damageMap[iCh][initialBpId + j].keys():
                                            typeDamage = self.damageMap[iCh][initialBpId + j][3].type
                                            if typeDamage <= 0:
                                                typeDamage = 0
                                            damageSpec += "4, " + str(j+1) + ", " + str(typeDamage) + " / "
                        damageSpec = damageSpec[:-2]
                        f.write(damageSpec + ";")
                        f.write('\n')
        return numSites

    def getParticleTimes(self):
        times = []
        for d in self.damageSitesAtCurrentTime:
            time = d.ParticleTime
            if time not in times:
                times.append(time)
        return times

    def getDamageCenterAndBoundaries(self, pos):
        xmax = pos[0][0]
        xmin = pos[0][0]
        ymax = pos[0][1]
        ymin = pos[0][1]
        zmax = pos[0][2]
        zmin = pos[0][2]
        for p in pos:
            if p[0] > xmax:
                xmax = p[0]
            if p[0] < xmin:
                xmin = p[0]
            if p[1] > ymax:
                ymax = p[1]
            if p[1] < ymin:
                ymin = p[1]
            if p[2] > zmax:
                zmax = p[2]
            if p[2] < zmin:
                zmin = p[2]
        xcenter = (xmax + xmin) / 2
        ycenter = (ymax + ymin) / 2
        zcenter = (zmax + zmin) / 2
        centerMaxMin = []
        centerMaxMin.append([xcenter, ycenter, zcenter])
        centerMaxMin.append([xmax, ymax, zmax])
        centerMaxMin.append([xmin, ymin, zmin])
        return centerMaxMin

    def getStringOutOfList(self, list, type='float'):
        string = ''
        for i, e in enumerate(list):
            if type == 'int':
                string += str(int(e))
            if type == 'float':
                string += str(e)
            if i < len(list) - 1:
                string += ", "
        return string

    def printDamageCount(self):
        self.messages.append("Summary of damage"); print(self.messages[-1])
        self.messages.append("-----------------"); print(self.messages[-1])
        self.messages.append("Dose " + str(self.cumulativeDose) + " Gy"); print(self.messages[-1])
        self.messages.append("DSB " + str(self.numDSB)); print(self.messages[-1])
        self.messages.append("DSB_Direct " + str(self.numDSBDirect)); print(self.messages[-1])
        self.messages.append("DSB_Indirect " + str(self.numDSBIndirect)); print(self.messages[-1])
        self.messages.append("DSB_Hybrid " + str(self.numDSBHybrid)); print(self.messages[-1])
        self.messages.append("SSB " + str(self.numSSB)); print(self.messages[-1])
        self.messages.append("SSB_Direct " + str(self.numSSBDirect)); print(self.messages[-1])
        self.messages.append("SSB_Indirect " + str(self.numSSBIndirect)); print(self.messages[-1])
        self.messages.append("SB " + str(self.numSB)); print(self.messages[-1])
        self.messages.append("SB_Direct " + str(self.numSBDirect)); print(self.messages[-1])
        self.messages.append("SB_Indirect " + str(self.numSBIndirect)); print(self.messages[-1])
        self.messages.append("BD " + str(self.numBD)); print(self.messages[-1])
        self.messages.append("BD_Direct " + str(self.numBDDirect)); print(self.messages[-1])
        self.messages.append("BD_Indirect " + str(self.numBDIndirect)); print(self.messages[-1])
        self.messages.append("DSB positions " + str(len(self.DSBPositions))); print(self.messages[-1])
        self.messages.append("Number of foci " + str(self.getNumberOfFoci(0.5))); print(self.messages[-1])

    def getDoseResponseCurve(self, plot=True, q='dsb'):
        if q.lower() == 'dsb':
            y = self.DSBarray
            ylabel = 'DSB'
        elif q.lower() == 'ssb':
            y = self.SSBarray
            ylabel = 'SSB'
        elif q.lower() == 'sb':
            y = self.SBarray
            ylabel = 'SB'
        elif q.lower() == 'bd':
            y = self.BDarray
            ylabel = 'BD'
        elif q.lower() == 'dsbd':
            y = self.DSBdirectarray
            ylabel = 'Direct DSB'
        elif q.lower() == 'dsbi':
            y = self.DSBindirectarray
            ylabel = 'Indirect DSB'
        elif q.lower() == 'dsbh':
            y = self.DSBhybridarray
            ylabel = 'Hybrid DSB'
        elif q.lower() == 'ssbd':
            y = self.SSBdirectarray
            ylabel = 'Direct SSB'
        elif q.lower() == 'ssbi':
            y = self.SSBindirectarray
            ylabel = 'Indirect SSB'
        elif q.lower() == 'sbd':
            y = self.SBdirectarray
            ylabel = 'Direct SB'
        elif q.lower() == 'sbi':
            y = self.SBindirectarray
            ylabel = 'Indirect SB'
        elif q.lower() == 'bdd':
            y = self.BDdirectarray
            ylabel = 'Direct BD'
        elif q.lower() == 'bdi':
            y = self.BDindirectarray
            ylabel = 'Indirect BD'
        elif q.lower() == 'nsites':
            y = self.numSitesarray
        if plot:
            fig = plt.figure()
            fig.set_size_inches((4, 4))
            ax = fig.add_subplot(111)
            ax.plot(self.Darray, y)
            ax.set_xlim(0, None)
            ax.set_ylim(0, None)
            ax.set_xlabel('Dose (Gy)')
            ax.set_ylabel(ylabel)
            ax.grid()
            fig.tight_layout()
            plt.show()
        return self.Darray, y

    def getDamageClusterDistribution(self, plot=True, savefile='', fromRecomputedSites=True):
        n = np.array([], dtype=int)
        if fromRecomputedSites:
            damageSitesPerTime = self.classifiedDamageSites
            for damageSites in damageSitesPerTime:
                for iCh in damageSites.keys():
                    ibpsTakenForThisChromosome = []
                    for i in range(len(damageSites[iCh])):
                        initialBpId = damageSites[iCh][i]
                        bd = 0; sb = 0; dsb = 0;
                        for j in range(self.nbpForDSB):
                            if initialBpId + j not in ibpsTakenForThisChromosome:
                                if iCh in self.SSBMap.keys():
                                    if initialBpId + j in self.SSBMap[iCh].keys():
                                        if 1 in self.SSBMap[iCh][initialBpId + j].keys():
                                            if self.SSBMap[iCh][initialBpId + j][1].type >= 1:
                                                sb += 1
                                        if 2 in self.SSBMap[iCh][initialBpId + j].keys():
                                            if self.SSBMap[iCh][initialBpId + j][2].type >= 1:
                                                sb += 1
                                if iCh in self.DSBMap.keys():
                                    if initialBpId + j in self.DSBMap[iCh].keys():
                                        if 1 in self.DSBMap[iCh][initialBpId + j].keys():
                                            if self.DSBMap[iCh][initialBpId + j][1].type >= 1:
                                                sb += 1
                                                dsb += 1
                                        if 2 in self.DSBMap[iCh][initialBpId + j].keys():
                                            if self.DSBMap[iCh][initialBpId + j][2].type >= 1:
                                                sb += 1
                                if iCh in self.BDMap.keys():
                                    if initialBpId + j in self.BDMap[iCh].keys():
                                        if 1 in self.BDMap[iCh][initialBpId + j].keys():
                                            if self.BDMap[iCh][initialBpId + j][1].type >= 1:
                                                bd += 1
                                        if 2 in self.BDMap[iCh][initialBpId + j].keys():
                                            if self.BDMap[iCh][initialBpId + j][2].type >= 1:
                                                bd += 1
                            ibpsTakenForThisChromosome.append(initialBpId + j)
                        if dsb > 0:
                            numdamages = sb + bd
                            n = np.append(n, numdamages-2)
        else:
            for damage in self.damageSitesAtCurrentTime:
                if damage.numberofDSBs > 0:
                    numdamages =  damage.numberofbasedamages + damage.numberofstrandbreaks
                    n = np.append(n, numdamages-2)

        if plot:
            # Plot the histogram
            plt.hist(n, 'auto', density=True, alpha=0.5, color='b')
            plt.title('Distribution of clusters')
            plt.xlabel('Number of damages in a cluster')
            plt.ylabel('Frequency')
            if not savefile == '':
                plt.savefig(savefile)
            else:# Show the plot
                plt.show()
        return n

    def getNumberOfFoci(self, fociSize = 0.4):
        indexIsAvailable = []
        for i in range(len(self.DSBPositions)):
            indexIsAvailable.append(True)
        vectorOfDSBsinEachFocus = []
        dsbIdsInThisFocus = []
        for i in range(len(self.DSBPositions)):
            if indexIsAvailable[i]:
                indexIsAvailable[i] = False
                dsbIdsInThisFocus.append(i)
                for j in range(len(self.DSBPositions)):
                    if indexIsAvailable[j] and self.getDistance(self.DSBPositions[i], self.DSBPositions[j]) < fociSize / 2:
                        indexIsAvailable[j] = False
                        dsbIdsInThisFocus.append(j)
                vectorOfDSBsinEachFocus.append(dsbIdsInThisFocus)
            dsbIdsInThisFocus = []
        return len(vectorOfDSBsinEachFocus)

    def getDistance(self, a, b):
        return np.sqrt(np.sum(np.power(a-b, 2)))

    def produce3DImage(self, show=True, microscopePSFWidth = 0.4, resolution = 0.4, xmin = -5, xmax = 5, ymin = -5, ymax = 5, zmin = -5, zmax = 5,
                       indQuantity='Dose', indValueString=''):
        dx = resolution
        dy = resolution
        dz = resolution
        nx = int(np.floor((xmax - xmin) / dx + 1e-14) + 1)
        ny = int(np.floor((ymax - ymin) / dy + 1e-14) + 1)
        nz = int(np.floor((zmax - zmin) / dz + 1e-14) + 1)
        img3d = np.zeros([nx, ny, nz])
        # Gets the PSF (Gaussian)
        halfsize = np.floor(3 * microscopePSFWidth / resolution)
        if halfsize < 3:
            halfsize = 3
        nkernel = int(2 * halfsize + 1)
        psf = np.zeros([nkernel, nkernel, nkernel])
        minkernel = -halfsize * resolution
        sum = 0
        for ix in range(nkernel):
            xpos = minkernel + ix * dx
            for iy in range(nkernel):
                ypos = minkernel + iy * dy
                for iz in range(nkernel):
                    zpos = minkernel + iz * dz
                    v = self.Gaussian3D(xpos, ypos, zpos, microscopePSFWidth)
                    psf[ix, iy, iz] = v
                    sum = sum + v
        # Normalize
        psf = psf / sum
        # Convolves PSF for each DSB position
        for iDsb in range(len(self.DSBPositions)):
            idx = int(np.floor((self.DSBPositions[iDsb][0] - xmin) / dx + 1e-14))
            idy = int(np.floor((self.DSBPositions[iDsb][1] - xmin) / dy + 1e-14))
            idz = int(np.floor((self.DSBPositions[iDsb][2] - xmin) / dz + 1e-14))
            startx = int(idx - halfsize)
            starty = int(idy - halfsize)
            startz = int(idz - halfsize)
            if startx < 0:
                startx = 0
            if starty < 0:
                starty = 0
            if startz < 0:
                startz = 0
            endx = int(idx + halfsize)
            endy = int(idy + halfsize)
            endz = int(idz + halfsize)
            if endx > nx:
                endx = nx
            if endy > ny:
                endy = ny
            if endz > nz:
                endz = nz
            for i in range(startx, endx):
                for j in range(starty, endy):
                    for k in range(startz, endz):
                        img3d[i, j, k] = img3d[i, j, k] + psf[i - startx, j - starty, k - startz]
        x = []
        y = []
        z = []
        v = []
        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    x.append(i)
                    y.append(j)
                    z.append(k)
                    v.append(img3d[i,j,k])
        x = np.array(x)
        y = np.array(y)
        z = np.array(z)
        v = np.array(v)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        if indQuantity == 'Dose':
            ax.set_title('Dose = ' + str(np.round(self.accumulateDose, 2)) + ' Gy')
        elif indValueString == '':
            ax.set_title()
        else:
            ax.set_title(indQuantity + ' = ' + indValueString)
        img = ax.scatter(x, y, z, c=v, cmap='MyColorMapAlpha', marker='s', s=20)
        fig.colorbar(img)
        if show:
            plt.show()
        else:
            d = os.getcwd() + '/RadDamDNA/video/'
            nfiles = len(os.listdir(d))
            plt.savefig(d + str(nfiles) + '.png')

    def produce2DImages(self, microscopePSFWidth = 0.8, resolution = 0.1, xmin = -5, xmax = 5, ymin = -5, ymax = 5):
        halfSize = int(np.floor(3*microscopePSFWidth / resolution))
        if halfSize < 3:
            halfSize = 3
        nkernel = 2 * halfSize + 1
        psf = np.zeros([nkernel, nkernel])
        minkernel = -halfSize * resolution
        sum = 0
        for ix in range(nkernel):
            xpos = minkernel + ix * resolution
            for iy in range(nkernel):
                ypos = minkernel + iy * resolution
                v = self.Gaussian2D(xpos, ypos, microscopePSFWidth)
                psf[ix, iy] = v
                sum = sum + v
        psf = psf / sum
        fig = plt.figure()
        # Different planes (0,1), (0,2) and (1,2)
        ids1 = [0, 0, 1]
        ids2 = [1, 2, 2]
        pos = [131, 132, 133]
        titles = ['Z-plane', 'Y-plane', 'X-plane']
        xaxislbl = [r'x ($\mu$m)', r'x ($\mu$m)', r'y ($\mu$m)']
        yaxislbl = [r'y ($\mu$m)', r'z ($\mu$m)', r'z ($\mu$m)']
        for id in range(len(ids1)):
            dx = resolution
            dy = resolution
            nx = int(np.floor((xmax - xmin) / dx + 1e-14) + 1)
            ny = int(np.floor((ymax - ymin) / dy + 1e-14) + 1)
            img2d = np.zeros([nx, ny])
            for iDsb in range(len(self.DSBPositions)):
                idx = int(np.floor((self.DSBPositions[iDsb][ids1[id]] - xmin)) / dx + 1e-14)
                idy = int(np.floor((self.DSBPositions[iDsb][ids2[id]] - ymin)) / dy + 1e-14)
                startx = int(idx - halfSize)
                starty = int(idy - halfSize)
                if startx < 0:
                    startx = 0
                if starty < 0:
                    starty = 0
                endx = int(idx + halfSize)
                endy = int(idy + halfSize)
                if endx > nx:
                    endx = nx
                if endy > ny:
                    endy = ny
                for i in range(startx, endx):
                    for j in range(starty, endy):
                        img2d[i, j] = img2d[i, j] + psf[i - startx, j - starty]
            ax = fig.add_subplot(pos[id])
            ax.imshow(img2d, cmap='MyColorMapAlpha', extent=[xmin, xmax, ymin, ymax], origin='lower')
            ax.set_title(titles[id])
            #ax.set_xticks([j for j in np.linspace(xmin, xmax, 5)])
            #ax.set_yticks([j for j in np.linspace(ymin, ymax, 5)])
            ax.set_xlabel(xaxislbl[id])
            ax.set_ylabel(yaxislbl[id])
            ax.set_aspect('equal')
        plt.show()

    def Gaussian3D(self, x, y, z, sigma):
        return np.exp(-(x ** 2 + y ** 2 + z ** 2) / (2 * sigma ** 2))

    def Gaussian2D(self, x, y, sigma):
        return np.exp(-(x ** 2 + y ** 2) / (2 * sigma ** 2))

    def SetUpColorMap(self):
        ncolors = 256
        color_array = plt.get_cmap('Reds')(range(ncolors))
        color_array[:, -1] = np.linspace(0.01, 0.99, ncolors)
        map_object = LinearSegmentedColormap.from_list(name='MyColorMapAlpha', colors=color_array)
        try:
            plt.register_cmap(cmap=map_object)
        except:
            pass

class SubcomponentLesion:
    def __init__(self, type, position, lesiontime, particletime, eventid=0, dsbid=None, complexity=None):
        self.type = type
        self.position = position
        self.lesiontime = lesiontime
        self.particletime = particletime
        self.eventID = eventid
        self.dsbID = dsbid
        self.complexity = complexity