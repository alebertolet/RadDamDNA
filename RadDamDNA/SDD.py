#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 10:32 AM

@author: alejandrobertolet
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class DamageToDNA:
    def __init__(self):
        self.damageSites = []
        self.doses = np.array([])
        self.accumulateDose = 0
        self.damageMap = {}
        self.Darray = []
        self.DSBarray = []

    def readSDDAndDose(self, path, namessd = 'DNADamage_sdd.txt', namephsp = 'DNADamage.phsp'):
        dosepath = path + namephsp
        f = open(dosepath, 'r')
        lines = f.readlines()
        for l in lines:
            split = l.split()
            self.doses = np.append(self.doses, float(split[1]))
        sddpath = path + namessd
        self.readFromSDD(sddpath)
        f.close()

    def readFromSDD(self, path):
        reader = SDDReader(path)
        for key in reader.headerProperties:
            setattr(self, key, reader.headerProperties[key])
        self.nbpForDSB = int(self.Damagedefinition[2])
        for d in reader.damages:
            dsite = DamageSite()
            for key in d:
                setattr(dsite, key, d[key])
            dsite.initialBp = int(dsite.ChromosomePosition * self.Chromosomesizes[dsite.chromosomeNumber]*1e6)
            self.damageSites.append(dsite)

    def populateDamages(self):
        # Key of damage map is the chromosome number
        iExposure = -1
        for damage in self.damageSites:
            iCh = damage.chromosomeNumber
            if iCh not in self.damageMap.keys():
                self.damageMap[iCh] = {}
            for bpdamage in damage.individualdamages:
                iBp = damage.initialBp+bpdamage['basepairID']
                if iBp not in self.damageMap[iCh].keys():
                    self.damageMap[iCh][iBp] = {}
                self.damageMap[iCh][damage.initialBp+bpdamage['basepairID']][bpdamage['subcomponent']] = bpdamage['type']
            if damage.newExposure > 1:
                iExposure += 1
                self.accumulateDose += self.doses[iExposure]
                self.Darray.append(self.accumulateDose)
                self.DSBarray.append(self.computeStrandBreaks())

    def computeStrandBreaks(self):
        DSBMap = {}
        DSBPairs = {}
        SSBMap = {}
        BDMap = {}
        for iCh in self.damageMap.keys():
            if iCh not in DSBMap.keys():
                DSBMap[iCh] = {}
                DSBPairs[iCh] = {}
            for iBp in self.damageMap[iCh].keys():
                if iBp not in DSBMap[iCh].keys():
                    DSBMap[iCh][iBp] = {}
                    DSBMap[iCh][iBp][1] = 0
                    DSBMap[iCh][iBp][2] = 0
                # Checks if there is damage in backbone 1 (iComp = 2) and there is not already a DSB identified in strand 1
                if 2 in self.damageMap[iCh][iBp].keys() and self.damageMap[iCh][iBp][2] > 0 and (1 not in DSBMap[iCh][iBp].keys() or DSBMap[iCh][iBp][1] == 0):
                    dsbFound = False
                    for i2 in range(0, self.nbpForDSB):
                        if iBp+i2 in self.damageMap[iCh].keys() and 3 in self.damageMap[iCh][iBp+i2].keys():
                            typeDamageInStrand2 = self.damageMap[iCh][iBp+i2][3]
                            if typeDamageInStrand2 > 0:
                                if iBp+i2 not in DSBMap[iCh].keys():
                                    DSBMap[iCh][iBp+i2] = {}
                                if 2 not in DSBMap[iCh][iBp+i2].keys() or DSBMap[iCh][iBp+i2][2] == 0:
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
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2+i1][2]
                                            if typeDamageInStrand1 > 0:
                                                if iBp+i2+i1 not in DSBMap[iCh].keys():
                                                    DSBMap[iCh][iBp + i2 + i1] = {}
                                                if 1 not in DSBMap[iCh][iBp+i2+i1].keys() or DSBMap[iCh][iBp+i2+i1][1] == 0:
                                                    closestPosInStrand1 = iBp + i2 + i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp+i2-i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp+i2-i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp+i2-i1][2]
                                            if iBp+i2-i1 >= 0 and typeDamageInStrand1 > 0:
                                                if iBp+i2-i1 not in DSBMap[iCh].keys():
                                                    DSBMap[iCh][iBp + i2 - i1] = {}
                                                if 1 not in DSBMap[iCh][iBp+i2-i1].keys() or DSBMap[iCh][iBp+i2-i1][1] == 0:
                                                    closestPosInStrand1 = iBp + i2 - i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                DSBMap[iCh][closestPosInStrand1][1] = adjustedTypeDamageInStrand1
                                DSBMap[iCh][iBp+i2][2] = adjustedTypeDamageInStrand2
                                pos = (closestPosInStrand1, iBp + i2)
                                DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                dsbFound = True
                        if iBp-i2 >= 0 and iBp-i2 in self.damageMap[iCh].keys() and 3 in self.damageMap[iCh][iBp-i2].keys():
                            typeDamageInStrand2 = self.damageMap[iCh][iBp - i2][3]
                            if i2 > 0 and typeDamageInStrand2 > 0:
                                if iBp-i2 not in DSBMap[iCh].keys():
                                    DSBMap[iCh][iBp-i2] = {}
                                if 2 not in DSBMap[iCh][iBp-i2].keys() or DSBMap[iCh][iBp-i2][2] == 0:
                                    adjustedTypeDamageInStrand2 = 2
                                    if typeDamageInStrand2 == 1 or typeDamageInStrand2 == 3:
                                        adjustedTypeDamageInStrand2 = 1
                                    closestPosInStrand1 = iBp
                                    adjustedTypeDamageInStrand1 = 2
                                    for i1 in range(0, self.nbpForDSB):
                                        if iBp-i2+i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp-i2+i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2+i1][2]
                                            if typeDamageInStrand1 > 0:
                                                if iBp-i2+i1 not in DSBMap[iCh].keys():
                                                    DSBMap[iCh][iBp - i2 + i1] = {}
                                                if 1 not in DSBMap[iCh][iBp-i2+i1].keys() or DSBMap[iCh][iBp-i2+i1][1] == 0:
                                                    closestPosInStrand1 = iBp-i2+i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                        if iBp-i2-i1 in self.damageMap[iCh].keys() and 2 in self.damageMap[iCh][iBp - i2 - i1].keys():
                                            typeDamageInStrand1 = self.damageMap[iCh][iBp-i2-i1][2]
                                            if iBp-i2-i1 >= 0 and i1 > 0 and typeDamageInStrand1 > 0:
                                                if iBp-i2-i1 not in DSBMap[iCh].keys():
                                                    DSBMap[iCh][iBp - i2 - i1] = {}
                                                if 1 not in DSBMap[iCh][iBp-i2-i1].keys() or DSBMap[iCh][iBp-i2-i1][1] == 0:
                                                    closestPosInStrand1 = iBp-i2-i1
                                                    if typeDamageInStrand1 == 1 or typeDamageInStrand1 == 3:
                                                        adjustedTypeDamageInStrand1 = 1
                                                    break
                                    DSBMap[iCh][closestPosInStrand1][1] = adjustedTypeDamageInStrand1
                                    DSBMap[iCh][iBp-i2][2] = adjustedTypeDamageInStrand2
                                    pos = (closestPosInStrand1, iBp - i2)
                                    DSBPairs[iCh][pos] = adjustedTypeDamageInStrand1 + adjustedTypeDamageInStrand2
                                    dsbFound = True
                        if dsbFound:
                            break
        ndsb = 0
        for iCh in DSBPairs.keys():
            ndsb = ndsb + len(DSBPairs[iCh])
        return ndsb
        #Loops through damage map to exclude

    def plotDoseResponseCurve(self):
        fig = plt.figure()
        fig.set_size_inches((4, 4))
        ax = fig.add_subplot(111)
        ax.plot(self.Darray, self.DSBarray)
        ax.set_xlim(0, None)
        ax.set_ylim(0, None)
        ax.grid()
        fig.tight_layout()
        plt.show()


class DamageSite:
    def __init__(self):
        self.initialBp = 0

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
        return self._chromosomeId
    @ChromosomeID.setter
    def ChromosomeID(self, v):
        self._chromosomeID = v
        self.typeOfChromatine = int(v[0])
        self.chromosomeNumber = int(v[1])-1  ### -1 IS FOR THE OLD VERSION OF THE SCORER; THIS HAS BEEN CORRECTED
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
        self.lesiontimes = []
        for lt in v:
            self.lesiontimes.append(lt)

    @property
    def ParticleTypes(self):
        return self._particleTypes
    @ParticleTypes.setter
    def ParticleTypes(self, v):
        self._particles = v
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
        if len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Energy'] = v[i]

    @property
    def Translation(self):
        return self._translation
    @Translation.setter
    def Translation(self, v):
        self._translation = v
        if len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Translation'] = v[i]

    @property
    def Direction(self):
        return self._direction
    @Direction.setter
    def Direction(self, v):
        self._direction = v
        if len(self.particles) == len(v):
            for i, p in enumerate(self.particles):
                p['Direction'] = v[i]
    @property
    def ParticleTime(self):
        return self._particleTime
    @ParticleTime.setter
    def ParticleTime(self, v):
        self._particleTime = v
        if len(self.particles) == len(v):
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