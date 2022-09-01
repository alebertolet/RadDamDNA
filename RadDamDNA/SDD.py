#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 10:32 AM

@author: alejandrobertolet
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import os
import random

class DamageToDNA:
    def __init__(self):
        self.initializeStructures()
        self.initializeCounters()
        self.SetUpColorMap()

    def initializeStructures(self):
        self.damageSites = []
        self.doses = np.array([])
        self.accumulateDose = 0
        self.damageMap = {}
        self.Darray = []
        self.DSBarray = []
        self.SSBarray = []
        self.SBarray = []
        self.BDarray = []
        self.DSBPositions = []

    def initializeCounters(self):
        self.numSB = 0; self.numSBDirect = 0; self.numSBIndirect = 0
        self.numSSB = 0; self.numSSBDirect = 0; self.numSSBIndirect = 0
        self.numDSB = 0; self.numDSBDirect = 0; self.numDSBIndirect = 0; self.numDSBHybrid = 0
        self.numBD = 0; self.numBDDirect = 0; self.numBDIndirect = 0
        self.numSSBPlus = 0; self.numDSBPlus = 0; self.numDSBComplex = 0

    def readSDDAndDose(self, path, namessd = 'DNADamage_sdd.txt', namephsp = 'DNADamage.phsp', defectiveChromosomeNumber=False, particleTime=0, lesionTime=0):
        if namephsp is not None:
            dosepath = path + namephsp
            f = open(dosepath, 'r')
            lines = f.readlines()
            for l in lines:
                split = l.split()
                self.doses = np.append(self.doses, float(split[1]))
        sddpath = path + namessd
        self.readFromSDD(sddpath, defectiveChromosomeNumber, particleTime, lesionTime)
        f.close()
        self.namessd = namessd
        self.namephsp = namephsp

    def readFromSDD(self, path, dcn=False, particleTime=0, lesionTime=0):
        reader = SDDReader(path)
        for key in reader.headerProperties:
            setattr(self, key.replace("/", "_"), reader.headerProperties[key])
        self.nbpForDSB = int(self.Damagedefinition[2])
        for d in reader.damages:
            dsite = SDDDamageSite(dcn)
            for key in d:
                setattr(dsite, key, d[key])
            dsite.initialBp = np.round(dsite.ChromosomePosition * self.Chromosomesizes[dsite.chromosomeNumber]*1e6)
            dsite.ParticleTime = particleTime
            dsite.LesionTime = lesionTime
            self.damageSites.append(dsite)

    def populateDamages(self, getVideo=False, stopAtDose=-1):
        # Key of damage map is the chromosome number
        iExposure = -1
        if getVideo:
            d = os.getcwd() + '/RadDamDNA/video/'
            for f in os.listdir(d):
                os.remove(os.path.join(d, f))
        iEvent = 0
        for damage in self.damageSites:
            if stopAtDose > 0 and self.accumulateDose > stopAtDose:
                break
            iCh = damage.chromosomeNumber
            if iCh not in self.damageMap.keys():
                self.damageMap[iCh] = {}
            for bpdamage in damage.individualdamages:
                iBp = int(damage.initialBp+bpdamage['basepairID'])
                if iBp not in self.damageMap[iCh].keys():
                    self.damageMap[iCh][iBp] = {}
                self.damageMap[iCh][damage.initialBp+bpdamage['basepairID']][bpdamage['subcomponent']] = \
                    SubcomponentLesion(bpdamage['type'], [damage.centerX, damage.centerY, damage.centerZ], damage.LesionTime, damage.ParticleTime, iEvent)
            if damage.newExposure > 1:
                iExposure += 1
                iEvent += 1
                self.computeStrandBreaks()
                if len(self.doses) > 0:
                    self.accumulateDose += self.doses[iExposure]
                    self.Darray.append(self.accumulateDose)
                    self.DSBarray.append(self.numDSB)
                    self.SSBarray.append(self.numSSB)
                    self.SBarray.append(self.numSB)
                    self.BDarray.append(self.numBD)
                    if getVideo:
                        self.produce3DImage(show=False)
            if damage.newExposure == 1:
                iEvent += 1

    def computeStrandBreaks(self):
        self.DSBMap = {}
        DSBPairs = {}
        self.SSBMap = {}
        self.BDMap = {}
        self.DSBPositions = []
        for iCh in self.damageMap.keys():
            if iCh not in self.DSBMap.keys():
                self.DSBMap[iCh] = {}
                DSBPairs[iCh] = {}
            for iBp in self.damageMap[iCh].keys():
                if iBp not in self.DSBMap[iCh].keys():
                    self.DSBMap[iCh][iBp] = {}
                    self.DSBMap[iCh][iBp][1] = SubcomponentLesion(0, None, 0, 0)
                    self.DSBMap[iCh][iBp][2] = SubcomponentLesion(0, None, 0, 0)
                # Checks if there is damage in backbone 1 (iComp = 2) and there is not already a DSB identified in strand 1
                if 2 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][2].type > 0 and (1 not in self.DSBMap[iCh][iBp].keys() or self.DSBMap[iCh][iBp][1].type == 0):
                    dsbFound = False
                    for i2 in range(0, self.nbpForDSB):
                        if iBp+i2 in self.damageMap[iCh].keys() and 3 in self.damageMap[iCh][iBp+i2].keys():
                            typeDamageInStrand2 = self.damageMap[iCh][iBp+i2][3].type
                            if typeDamageInStrand2 > 0:
                                if iBp+i2 not in self.DSBMap[iCh].keys():
                                    self.DSBMap[iCh][iBp+i2] = {}
                                if 2 not in self.DSBMap[iCh][iBp+i2].keys() or self.DSBMap[iCh][iBp+i2][2].type == 0:
                                    adjustedTypeDamageInStrand2 = 2
                                    # Checks if damage in backbone 2 is direct and keeps it in that case
                                    # If it is multiple damage, considers direct (since it would not happen without chemistry)
                                    if typeDamageInStrand2 == 1 or typeDamageInStrand2 == 3:
                                        adjustedTypeDamageInStrand2 = 1

                                    # Once damage is identified in strand 2 (and damage ensured in strand 1), looks for the closest damage in strand 1
                                    closestPosInStrand1 = iBp
                                    adjustedTypeDamageInStrand1 = 2
                                    for i1 in range(0, self.nbpForDSB):
                                        if iBp+i2+i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp+i2+i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2+i1][2].type
                                            if typeDamageInStrand1 > 0:
                                                if iBp+i2+i1 not in self.DSBMap[iCh].keys():
                                                    self.DSBMap[iCh][iBp + i2 + i1] = {}
                                                if 1 not in self.DSBMap[iCh][iBp+i2+i1].keys() or self.DSBMap[iCh][iBp+i2+i1][1].type == 0:
                                                    closestPosInStrand1 = iBp + i2 + i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp+i2-i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp+i2-i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2-i1][2].type
                                            if iBp+i2-i1 >= 0 and typeDamageInStrand1 > 0:
                                                if iBp+i2-i1 not in self.DSBMap[iCh].keys():
                                                    self.DSBMap[iCh][iBp + i2 - i1] = {}
                                                if 1 not in self.DSBMap[iCh][iBp+i2-i1].keys() or self.DSBMap[iCh][iBp+i2-i1][1].type == 0:
                                                    closestPosInStrand1 = iBp + i2 - i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                self.DSBMap[iCh][closestPosInStrand1][1] = SubcomponentLesion(0, None, 0, 0)
                                self.DSBMap[iCh][closestPosInStrand1][1].type = adjustedTypeDamageInStrand1
                                self.DSBMap[iCh][closestPosInStrand1][1].particletime = self.damageMap[iCh][closestPosInStrand1][2].particletime
                                self.DSBMap[iCh][iBp+i2][2] = SubcomponentLesion(0, None, 0, 0)
                                self.DSBMap[iCh][iBp+i2][2].type = adjustedTypeDamageInStrand2
                                self.DSBMap[iCh][iBp+i2][2].particletime = self.damageMap[iCh][iBp+i2][3].particletime
                                pos = (closestPosInStrand1, iBp + i2)
                                DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                posSB1 = np.array(self.damageMap[iCh][closestPosInStrand1][2].position)
                                posSB2 = np.array(self.damageMap[iCh][iBp+i2][3].position)
                                self.DSBPositions.append((posSB1 + posSB2)/2)
                                dsbFound = True
                        if iBp-i2 >= 0 and iBp-i2 in self.damageMap[iCh].keys() and 3 in self.damageMap[iCh][iBp-i2].keys():
                            typeDamageInStrand2 = self.damageMap[iCh][iBp - i2][3].type
                            if i2 > 0 and typeDamageInStrand2 > 0:
                                if iBp-i2 not in self.DSBMap[iCh].keys():
                                    self.DSBMap[iCh][iBp-i2] = {}
                                if 2 not in self.DSBMap[iCh][iBp-i2].keys() or self.DSBMap[iCh][iBp-i2][2].type == 0:
                                    adjustedTypeDamageInStrand2 = 2
                                    if typeDamageInStrand2 == 1 or typeDamageInStrand2 == 3:
                                        adjustedTypeDamageInStrand2 = 1
                                    closestPosInStrand1 = iBp
                                    adjustedTypeDamageInStrand1 = 2
                                    for i1 in range(0, self.nbpForDSB):
                                        if iBp-i2+i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp-i2+i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2+i1][2].type
                                            if typeDamageInStrand1 > 0:
                                                if iBp-i2+i1 not in self.DSBMap[iCh].keys():
                                                    self.DSBMap[iCh][iBp - i2 + i1] = {}
                                                if 1 not in self.DSBMap[iCh][iBp-i2+i1].keys() or self.DSBMap[iCh][iBp-i2+i1][1].type == 0:
                                                    closestPosInStrand1 = iBp-i2+i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp-i2-i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp - i2 - i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2-i1][2].type
                                            if iBp-i2-i1 >= 0 and i1 > 0 and typeDamageInStrand1 > 0:
                                                if iBp-i2-i1 not in self.DSBMap[iCh].keys():
                                                    self.DSBMap[iCh][iBp - i2 - i1] = {}
                                                if 1 not in self.DSBMap[iCh][iBp-i2-i1].keys() or self.DSBMap[iCh][iBp-i2-i1][1].type == 0:
                                                    closestPosInStrand1 = iBp-i2-i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                    self.DSBMap[iCh][closestPosInStrand1][1] = SubcomponentLesion(0, None, 0, 0)
                                    self.DSBMap[iCh][closestPosInStrand1][1].type = adjustedTypeDamageInStrand1
                                    self.DSBMap[iCh][closestPosInStrand1][1].particletime = self.damageMap[iCh][closestPosInStrand1][2].particletime
                                    self.DSBMap[iCh][iBp-i2][2] = SubcomponentLesion(0, None, 0, 0)
                                    self.DSBMap[iCh][iBp-i2][2].type = adjustedTypeDamageInStrand2
                                    self.DSBMap[iCh][iBp-i2][2].particletime = self.damageMap[iCh][iBp-i2][3].particletime
                                    pos = (closestPosInStrand1, iBp - i2)
                                    DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                    posSB1 = np.array(self.damageMap[iCh][closestPosInStrand1][2].position)
                                    posSB2 = np.array(self.damageMap[iCh][iBp - i2][3].position)
                                    self.DSBPositions.append((posSB1 + posSB2) / 2)
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
                        self.SSBMap[iCh][iBp][1] = SubcomponentLesion(self.damageMap[iCh][iBp][2].type, None, self.damageMap[iCh][iBp][2].lesiontime, self.damageMap[iCh][iBp][2].particletime)
                    if 3 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][3].type > 0:
                        self.SSBMap[iCh][iBp][2] = SubcomponentLesion(self.damageMap[iCh][iBp][3].type, None, self.damageMap[iCh][iBp][3].lesiontime, self.damageMap[iCh][iBp][3].particletime)
                # Base damages
                if 1 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][1].type > 0:
                    self.BDMap[iCh][iBp][1] = SubcomponentLesion(self.damageMap[iCh][iBp][1].type, None, self.damageMap[iCh][iBp][1].lesiontime, self.damageMap[iCh][iBp][1].particletime)
                if 4 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][4].type > 0:
                    self.BDMap[iCh][iBp][2] = SubcomponentLesion(self.damageMap[iCh][iBp][4].type, None, self.damageMap[iCh][iBp][4].lesiontime, self.damageMap[iCh][iBp][4].particletime)
        self.quantifyDamage(DSBPairs, self.SSBMap, self.BDMap)

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
        times = self.getParticleTimes()
        listOfChrAndIniBpIdPertTime = []
        for time in times:
            damageScore = {}
            for iCh in self.BDMap.keys():
                if iCh not in damageScore.keys():
                    damageScore[iCh] = {}
                for iBp in self.BDMap[iCh].keys():
                    if iBp not in damageScore[iCh].keys():
                        damageScore[iCh][iBp] = {}
                    for iCo in self.BDMap[iCh][iBp].keys():
                        if iCo not in damageScore[iCh][iBp].keys():
                            damageScore[iCh][iBp][iCo] = 0
                        if self.BDMap[iCh][iBp][iCo].type > 0 and self.BDMap[iCh][iBp][iCo].particletime == time:
                            damageScore[iCh][iBp][iCo] += 0.001
            for iCh in self.SSBMap.keys():
                if iCh not in damageScore.keys():
                    damageScore[iCh] = {}
                for iBp in self.SSBMap[iCh].keys():
                    if iBp not in damageScore[iCh].keys():
                        damageScore[iCh][iBp] = {}
                    for iCo in self.SSBMap[iCh][iBp].keys():
                        if iCo not in damageScore[iCh][iBp].keys():
                            damageScore[iCh][iBp][iCo] = 0
                        if self.SSBMap[iCh][iBp][iCo].type > 0 and self.SSBMap[iCh][iBp][iCo].particletime == time:
                            damageScore[iCh][iBp][iCo] += 0.1
            for iCh in self.DSBMap.keys():
                if iCh not in damageScore.keys():
                    damageScore[iCh] = {}
                for iBp in self.DSBMap[iCh].keys():
                    if iBp not in damageScore[iCh].keys():
                        damageScore[iCh][iBp] = {}
                    for iCo in self.DSBMap[iCh][iBp].keys():
                        if iCo not in damageScore[iCh][iBp].keys():
                            damageScore[iCh][iBp][iCo] = 0
                        if self.DSBMap[iCh][iBp][iCo].type > 0 and self.DSBMap[iCh][iBp][iCo].particletime == time:
                            damageScore[iCh][iBp][iCo] += 5.0
            chromosomeAndInitialBpIdOfDamageSites = {}
            for iCh in damageScore.keys():
                linearAccumulatedDamage = []
                initialBpIds = []
                startingBpIdForDamageSites = []
                maxDamage = 10 * self.nbpForDSB
                minDamage = 1e-8
                while maxDamage > minDamage:
                    for iBp in damageScore[iCh].keys():
                        sumDamage = 0
                        for i in range(self.nbpForDSB):
                            if iBp + i in damageScore[iCh].keys():
                                if 1 in damageScore[iCh][iBp + i].keys():
                                    sumDamage += damageScore[iCh][iBp + i][1]
                                if 2 in damageScore[iCh][iBp + i].keys():
                                    sumDamage += damageScore[iCh][iBp + i][2]
                        linearAccumulatedDamage.append(sumDamage)
                        initialBpIds.append(iBp)
                    m = 0
                    index = 0
                    for i in range(len(linearAccumulatedDamage)):
                        if m < linearAccumulatedDamage[i]:
                            m = linearAccumulatedDamage[i]
                            index = initialBpIds[i]
                    maxDamage = m
                    if maxDamage < minDamage:
                        break
                    startingBpIdForDamageSites.append(index)
                    for i in range(self.nbpForDSB):
                        if index + i in damageScore[iCh].keys():
                            damageScore[iCh][index + i][1] = 0
                            damageScore[iCh][index + i][2] = 0
                    linearAccumulatedDamage = []
                    initialBpIds = []
                chromosomeAndInitialBpIdOfDamageSites[iCh] = startingBpIdForDamageSites
            listOfChrAndIniBpIdPertTime.append(chromosomeAndInitialBpIdOfDamageSites)
        return listOfChrAndIniBpIdPertTime

    def writeSDD(self, filename):
        self.outputSDDHeader(filename)
        self.outputSDDFile(filename)

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
            f.write("Incident particles, " + str(int(self.Incidentparticles))  + ';\n')
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

    def outputSDDFile(self, filename):
        damageSitesPerTime = self.classifyDamageSites()
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
                        damageSpec = damageSpec[:-2]
                        f.write(damageSpec + ";")
                        f.write('\n')
        return numSites

    def computeMedrasBreaks(self):
        self.medrasBreaks = []
        self.populateDamages()
        self.computeStrandBreaks()
        self.emptySets = 0
        self.complexities = []
        damageSitesPerEvent = self.classifyDamageSites()
        numSites = 0
        firstExposure = True
        for damageSites in damageSitesPerEvent:
            complexity = 0
            setDSB = 0
            newEvent = True
            self.medrasBreaks.append([])
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
                                            setDSB += 1
                                            pos1 = self.damageMap[iCh][initialBpId + j][2].position
                                            lesionTime = self.damageMap[iCh][initialBpId + j][2].lesiontime
                                        if self.DSBMap[iCh][initialBpId + j][1].type == 2:
                                            indir += 1
                                            sb += 1
                                            dsb += 1
                                            setDSB += 1
                                            pos1 = self.damageMap[iCh][initialBpId + j][2].position
                                            lesionTime = self.damageMap[iCh][initialBpId + j][2].lesiontime
                                    if 2 in self.DSBMap[iCh][initialBpId + j].keys():
                                        if self.DSBMap[iCh][initialBpId + j][2].type == 1 or self.DSBMap[iCh][initialBpId + j][2].type == 3:
                                            dir += 1
                                            sb += 1
                                            pos2 = self.damageMap[iCh][initialBpId + j][3].position
                                        if self.DSBMap[iCh][initialBpId + j][2].type == 2:
                                            indir += 1
                                            sb += 1
                                            pos2 = self.damageMap[iCh][initialBpId + j][3].position
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
                    if dsb == 0:
                        break
                    if dsb + sb + bd > 2:
                        complexBreak = True
                        complexity += 1
                    else:
                        complexBreak = False
                    newExposureFlag = 0
                    if newEvent:
                        newExposureFlag = 1
                        newEvent = False
                    if firstExposure:
                        newExposureFlag = 2
                        firstExposure = False
                    chromosomeLength = self.Chromosomesizes[iCh]
                    damageChromPos = initialBpId / chromosomeLength / 1e6
                    chromID = [iCh, 1, 0]
                    typedamage = 0
                    if dir == 0 and indir > 0:
                        typedamage = 1
                    if dir > 0 and indir > 0:
                        typedamage = 2
                    cause = [typedamage, dir, indir]
                    # Break is: index, position, complexity, chromosome ID, upstream/downstream, new event status, time and cause
                    medrasbreak1 = [numSites, pos1, complexBreak, chromID[:], damageChromPos, -1, newExposureFlag, lesionTime, cause]
                    medrasbreak2 = [numSites, pos2, complexBreak, chromID[:], damageChromPos, 1, 0, lesionTime, cause]
                    self.medrasBreaks[-1] += [medrasbreak1, medrasbreak2]
                    numSites += 1
            self.complexities.append(complexity)
            if setDSB < 1:
                self.medrasBreaks.pop()
                self.complexities.pop()
                self.emptySets += 1
        complexFrac = [1.0 * c / (len(b) / 2.0) for c,b in zip(self.complexities, self.medrasBreaks)]
        self.meanComplexity = np.mean(complexFrac)
        return numSites

    def getParticleTimes(self):
        times = []
        for iCh in self.damageMap.keys():
            for iBp in self.damageMap[iCh].keys():
                for iCo in self.damageMap[iCh][iBp].keys():
                    time = self.damageMap[iCh][iBp][iCo].particletime
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
        print("Summary of damage")
        print("-----------------")
        print("Dose", self.accumulateDose, "Gy")
        print("DSB", self.numDSB)
        print("DSB_Direct", self.numDSBDirect)
        print("DSB_Indirect", self.numDSBIndirect)
        print("DSB_Hybrid", self.numDSBHybrid)
        print("SSB", self.numSSB)
        print("SSB_Direct", self.numSSBDirect)
        print("SSB_Indirect", self.numSSBIndirect)
        print("SB", self.numSB)
        print("SB_Direct", self.numSBDirect)
        print("SB_Indirect", self.numSBIndirect)
        print("BD", self.numBD)
        print("BD_Direct", self.numBDDirect)
        print("BD_Indirect", self.numBDIndirect)
        print("DSB positions", len(self.DSBPositions))
        print("Number of foci", self.getNumberOfFoci(0.5))

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

    def produce3DImage(self, show=True, microscopePSFWidth = 0.4, resolution = 0.4, xmin = -5, xmax = 5, ymin = -5, ymax = 5, zmin = -5, zmax = 5):
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
        ax.set_title('Dose = ' + str(np.round(self.accumulateDose, 2)) + ' Gy')
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
        plt.register_cmap(cmap=map_object)

class SubcomponentLesion:
    def __init__(self, type, position, lesiontime, particletime, eventid=0):
        self.type = type
        self.position = position
        self.lesiontime = lesiontime
        self.particletime = particletime
        self.eventID = eventid

class SDDDamageSite:
    def __init__(self, defectiveChromosomeNumber=False):
        self.initialBp = 0
        self.defectiveChromosomeNumber = defectiveChromosomeNumber
        self.particles = []

    def PrintProperties(self):
        temp = vars(self)
        for v in temp:
            print(v, " : ", temp[v])

    @property
    def Classification(self):
        return self._classification
    @Classification.setter
    def Classification(self, c):
        self._classification = c
        self.newExposure = int(c[0])
        self.eventID = int(c[1])

    @property
    def SpatialCoordinatesAndExtent(self):
        return self._coords
    @SpatialCoordinatesAndExtent.setter
    def SpatialCoordinatesAndExtent(self, v):
        self._coords = v
        if len(v) >= 3:
            self.centerX = v[0]
            self.centerY = v[1]
            self.centerZ = v[2]
        if len(v) >= 6:
            self.maxX = v[3]
            self.maxY = v[4]
            self.maxZ = v[5]
        if len(v) == 9:
            self.minX = v[6]
            self.minY = v[7]
            self.minZ = v[8]

    @property
    def ChromosomeID(self):
        return self._chromosomeID
    @ChromosomeID.setter
    def ChromosomeID(self, v):
        self._chromosomeID = v
        self.typeOfChromatine = int(v[0])
        if self.defectiveChromosomeNumber:
            self.chromosomeNumber = int(v[1])-1  ### -1 IS FOR THE OLD VERSION OF THE SCORER; THIS HAS BEEN CORRECTED
        else:
            self.chromosomeNumber = int(v[1])
        self.chromatideNumber = int(v[2])
        self.chromosomeArm = int(v[3])

    @property
    def ChromosomePosition(self):
        return self._chromosomePos
    @ChromosomePosition.setter
    def ChromosomePosition(self, v):
        self._chromosomePos = v

    @property
    def Cause(self):
        return self._cause
    @Cause.setter
    def Cause(self, v):
        self._cause = v
        self.damagetype = int(v[0])
        self.numberofdirectdamages = int(v[1])
        self.numberofindirectdamages = int(v[2])

    @property
    def DamageTypes(self):
        return self._damageTypes
    @DamageTypes.setter
    def DamageTypes(self, v):
        self._damageTypes = v
        self.numberofbasedamages = int(v[0])
        self.numberofstrandbreaks = int(v[1])
        self.numberofDSBs = int(v[2])

    @property
    def FullBreakSpec(self):
        return self._fullBreakSpec
    @FullBreakSpec.setter
    def FullBreakSpec(self, v):
        self._fullBreakSpec = v
        self.ndamages = int(len(v)/3)
        self.individualdamages = []
        for i in range(0, self.ndamages):
            damage = {}
            damage['subcomponent'] = int(v[i*3])
            damage['basepairID'] = int(v[i*3+1])
            damage['type'] = int(v[i*3+2])
            if damage['type'] == 4:
                damage['type'] = 1
            if damage['type'] == 5:
                damage['type'] = 3
            self.individualdamages.append(damage)

    @property
    def DNASequence(self):
        return self._dnaseq
    @DNASequence.setter
    def DNASequence(self, v):
        self._dnaseq = v

    @property
    def LesionTime(self):
        return self._lesiontime
    @LesionTime.setter
    def LesionTime(self, v):
        self._lesiontime = v
        if type(v) is list:
            self.lesiontimes = []
            for lt in v:
                self.lesiontimes.append(lt)

    @property
    def ParticleTypes(self):
        return self._particles
    @ParticleTypes.setter
    def ParticleTypes(self, v):
        self._particles = v
        if type(v) is list:
            self.particles = []
            for pt in v:
                particle = {}
                particle['Type'] = pt
                self.particles.append(particle)

    @property
    def Energies(self):
        return self._energies
    @Energies.setter
    def Energies(self, v):
        self._energies = v
        if type(v) is list and len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Energy'] = v[i]

    @property
    def Translation(self):
        return self._translation
    @Translation.setter
    def Translation(self, v):
        self._translation = v
        if type(v) is list and len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Translation'] = v[i]

    @property
    def Direction(self):
        return self._direction
    @Direction.setter
    def Direction(self, v):
        self._direction = v
        if type(v) is list and len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Direction'] = v[i]
    @property
    def ParticleTime(self):
        return self._particleTime
    @ParticleTime.setter
    def ParticleTime(self, v):
        self._particleTime = v
        if type(v) is list and len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Time'] = v[i]

class SDDReader:
    def __init__(self, path):
        self.file = open(path, 'r')
        self.__splitLines()
        self.file.close()
        self.__readHeader()
        self.__readDamage()

    def __splitLines(self):
        lines = self.file.readlines()
        self.headerlines = []
        self.damagelines = []
        iLine = 1
        for line in lines:
            if '***EndOfHeader***' in line:
                break
            else:
                self.headerlines.append(line)
                iLine += 1
        for i in range(iLine, len(lines)):
            if lines[i] != '\n':
                self.damagelines.append(lines[i])

    def __readHeader(self):
        self.headerProperties = {}
        for l in self.headerlines:
            split = l.split(',')
            if len(split) == 2:
                try:
                    self.headerProperties[split[0].strip().replace(' ', '')] = float(split[1].strip().replace(';',''))
                except:
                    self.headerProperties[split[0].strip().replace(' ', '')] = split[1].strip().replace(';', '')
            else:
                fields = []
                for i in range(1, len(split)):
                    try:
                        fields.append(float(split[i].strip().replace(';','')))
                    except:
                        fields.append(split[i].strip().replace(';', ''))
                self.headerProperties[split[0].strip().replace(' ', '')] = fields
        self.headerProperties['ScoringVolume'] = [self.headerProperties['Volumes'][7]]+list(map(float,self.headerProperties['Volumes'][8:]))
        self.dataEntries = []
        if self.headerProperties['Dataentries'][0] == 1:
            self.dataEntries.append('Classification')
        if self.headerProperties['Dataentries'][1] == 1:
            self.dataEntries.append('SpatialCoordinatesAndExtent')
        if self.headerProperties['Dataentries'][2] == 1:
            self.dataEntries.append('ChromosomeID')
        if self.headerProperties['Dataentries'][3] == 1:
            self.dataEntries.append('ChromosomePosition')
        if self.headerProperties['Dataentries'][4] == 1:
            self.dataEntries.append('Cause')
        if self.headerProperties['Dataentries'][5] == 1:
            self.dataEntries.append('DamageTypes')
        if self.headerProperties['Dataentries'][6] == 1:
            self.dataEntries.append('FullBreakSpec')
        if self.headerProperties['Dataentries'][7] == 1:
            self.dataEntries.append('DNASequence')
        if self.headerProperties['Dataentries'][8] == 1:
            self.dataEntries.append('LesionTime')
        if self.headerProperties['Dataentries'][9] == 1:
            self.dataEntries.append('ParticleTypes')
        if self.headerProperties['Dataentries'][10] == 1:
            self.dataEntries.append('Energies')
        if self.headerProperties['Dataentries'][11] == 1:
            self.dataEntries.append('Translation')
        if self.headerProperties['Dataentries'][12] == 1:
            self.dataEntries.append('Direction')
        if self.headerProperties['Dataentries'][13] == 1:
            self.dataEntries.append('ParticleTime')

    def __readDamage(self):
        self.damages = []
        for l in self.damagelines:
            newDamage = {}
            split = l.replace(';\n','').split(';')
            if len(split) != len(self.dataEntries):
                print("Inconsistency between data entries and damages read line by line.")
            else:
                for i in range(0, len(split)):
                    splittwice = split[i].split(',')
                    if len(splittwice) == 1:
                        try:
                            newDamage[self.dataEntries[i]] = float(splittwice[0].strip())
                        except:
                            newDamage[self.dataEntries[i]] = splittwice
                    else:
                        subfields = []
                        for subf in splittwice:
                            splitthree = subf.split('/')
                            if len(splitthree) > 1:
                                for subsubf in splitthree:
                                    try:
                                        subfields.append(float(subsubf.strip()))
                                    except:
                                        subfields.append(subsubf)
                            else:
                                try:
                                    subfields.append(float(subf.strip()))
                                except:
                                    subfields.append(subf)
                        newDamage[self.dataEntries[i]] = subfields
            self.damages.append(newDamage)

    def printHeader(self):
        print(self.headerProperties)

    def printDamages(self):
        print(self.damages)