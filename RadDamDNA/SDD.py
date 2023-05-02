#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/16/22 10:32 AM

@author: alejandrobertolet
"""

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

class SDDDamageSite:
    def __init__(self, version='2.0'):
        self.initialBp = 0
        self.version = version
        self.particles = []
        self.LesionTime = 0
        self.ParticleTime = 0
        self.isRepaired = False

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
        if float(self.version) < 2.0:
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
        self.isThereDSB = int(v[2])

    @property
    def FullBreakSpec(self):
        return self._fullBreakSpec
    @FullBreakSpec.setter
    def FullBreakSpec(self, v):
        self._fullBreakSpec = v
        self.ndamages = int(len(v)/3)
        self.individualdamages = []
        nbreaksStrand1 = 0
        nbreaksStrand2 = 0
        for i in range(0, self.ndamages):
            damage = {}
            subid = int(v[i*3])
            if float(self.version) < 2.0:
                damage['subcomponent'] = subid
            else:
                if subid == 1:
                    damage['subcomponent'] = 2
                    nbreaksStrand1 += 1
                elif subid == 2:
                    damage['subcomponent'] = 1
                elif subid == 3:
                    damage['subcomponent'] = 4
                elif subid == 4:
                    damage['subcomponent'] = 3
                    nbreaksStrand2 += 1
            damage['basepairID'] = int(v[i*3+1])
            damage['type'] = int(v[i*3+2])
            if damage['type'] == 4:
                damage['type'] = 1
            if damage['type'] == 5:
                damage['type'] = 3
            self.individualdamages.append(damage)
        # Count number of DSB
        self.numberofDSBs = min([nbreaksStrand1, nbreaksStrand2])

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