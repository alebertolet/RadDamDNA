#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/31/22 9:26 AM
Reimplementation of the MEDRAS-MC model developed by
Dr. Stephen McMahon (stephen.mcmahon@qub.ac.uk)
@author: alejandrobertolet (abertoletreina@mgh.harvard.edu)
"""

import numpy as np
from scipy import stats, spatial
import copy
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.patches
import os
haveDisplay = "DISPLAY" in os.environ
if not haveDisplay:
	matplotlib.use('Agg')
import re
from RadDamDNA.damage import DamageToDNA

class MedrasRepair:
    def __init__(self,  damage=None, sigma=0.04187, maxExposures=1000, repeats=50, minMisrepSize=0, writeKinetics=True, writeAllKinetics=False,
                 addFociDelay=True, kineticLimit=25, doPlot=False, allFragments=False, listAcentrics=False, simulationLimit=24):
        '''
        :param sigma: Misrepair range, as function of nuclear radius
        :param maxExposures: Maximum exposures per file to simulate
        :param repeats: Number of repeats for each exposure
        :param minMisrepSize: Neglect misrepair events separated by less than this genetic distance (MBP)
        :param writeKinetics: Print kinetics (simplified version)
        :param writeAllKinetics: Print all kinetics
        :param addFociDelay: Adds foci delay
        :param kineticLimit: Hours, maximum time for repair kinetics
        :param doPlot: Plot
        :param allFragments: Consider all fragments
        :param listAcentrics: List acentric aberrations
        :param simulationLimit: Hours, time at which to simulate misrepair
        '''
        self.damage = damage
        self.sigma = sigma
        self.maxExposures = maxExposures
        self.repeats = repeats
        self.minMisrepSize = minMisrepSize
        self.writeKinetics = writeKinetics
        self.writeAllKinetics = writeAllKinetics
        self.addFociDelay = addFociDelay
        self.kineticLimit = kineticLimit
        self.doPlot = doPlot
        self.allFragments = allFragments
        self.listAcentrics = listAcentrics
        self.simulationLimit = simulationLimit
        self.analysisFunctions = {}
        self.analysisFunctions['Fidelity'] = self.repairFidelity
        self.analysisFunctions['Separation'] = self.misrepairSeparation
        self.analysisFunctions['Spectrum'] = self.misrepairSpectrum
        self.analysisFunctions['DSBSeparation'] = self.dsbSeparation
        self.analysisFunctions['DSBRadial'] = self.radialDSBs
        self.analysisFunctions['TrackBreaks'] = self.trackBreaks
        self.fidelityRun = False
        self.separationRun = False
        self.radialRun = False
        self.trackRun = False

    def repairSimulation(self, path='', type='Fidelity', recalculateDamages=True):
        summaries = []
        if self.damage is not None:
            self.damage.computeMedrasBreaks(recalculateDamages)
        elif path != '':
            if path[-1] != '/':
                path += '/'
            fileNames = os.listdir(path)
            self.__sortNicely(fileNames)
            filePaths = [path + f for f in fileNames]
            for filePath in filePaths:
                if os.path.isdir(filePath) or (filePath[-7:] != 'sdd.txt' and filePath[-3:] != 'sdd'):
                    continue
                self.damage = DamageToDNA()
                self.damage.readFromSDD(filePath)
                self.damage.computeMedrasBreaks()
                # Find matching method and run
        else:
            print('Neither a valid path or damage object were provided!')
            return
        for name in self.analysisFunctions.keys():
            if type == name:
                summaries.append(self.analysisFunctions[name]())
        if len(summaries) == 0:
            print('Did not find matching analysis function. Options are:')
            print(', '.join(f[0] for f in self.analysisFunctions))
            return
        if len(summaries) == 0:
            print('No output returned!')
            return
        if summaries[0] is not None:
            for summary in summaries:
                print(summary)

    def repairFidelity(self):
        calcMR = MisrepairCalculator()
        # Print output header
        print('Break Set\tBreak Count\tMisrepair\tStdev\tInter-Chromosome Rate', end='')
        if self.writeAllKinetics:
            print('\t\t', '\t'.join(map(str, [tau/10.0 for tau in range(10*self.kineticLimit)])), end='')
        print()

        # Extract file data
        scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
        radius = scaledSigma / self.sigma
        chromSizes = self.damage.Chromosomesizes[1:]
        complexity = self.damage.meanComplexity

        outputs = [[0,0,0,0,0,0,0]] * self.damage.emptySets
        outputTimes = []
        for n, breakList in enumerate(self.damage.medrasBreaks):
            if n >= self.maxExposures: break
            print(str(n), '\t', len(breakList)/2, '\t', end='')
            if len(breakList) <= 2 or len(breakList) > 20000:
                repairData = [0, 0, 0, [], 0, 0]
                if len(breakList) > 20000:
                    print("Skipping due to memory concerns", '\t', end='')
            else:
                repairData = calcMR.fullRepair(breakList, scaledSigma, repeats=self.repeats,
                                               addFociClearance=self.addFociDelay, radius=radius,
                                               chromSizes=chromSizes, sizeLimit=self.minMisrepSize)
            # Unpack repair data
            mean, stdevRate, interChromRate, repTimes, analyticMisrep, etaSum = repairData
            outputs.append([len(breakList) / 2.0, mean, stdevRate, interChromRate, analyticMisrep, etaSum])
            outputTimes += repTimes
            # Write misrepair and kinetic data
            print(mean, '\t', stdevRate, '\t', interChromRate, end='')
            if self.writeAllKinetics:
                print('\t\t', self.summarizeKinetics(repTimes), end='')
            print()
        print()
        # Build summary
        summary = ''
        if self.fidelityRun is False:
            self.fidelityRun = True
            summary += ('\tBreak Set\tTotal Breaks\tComplexity\tBreaks per Exposure\tBreaks Stdev'
                        '\tMisrepair\tStdev\tInter-Chromosome Rate')
            if self.writeKinetics:
                summary += ('\t\tTime\t' + '\t'.join(str(tau/10.0) for tau in range(10*self.kineticLimit)))
            summary += '\n'
        summary += self.summarizeFidelity(complexity, outputs)
        if self.writeKinetics:
            summary += '\t\t\t'
            summary += self.summarizeKinetics(outputTimes)
        return summary

    def setVideoForOneEvent(self, finalTime=np.inf, recalculateDamages=True):
        calcMR = MisrepairCalculator()
        if self.damage is not None:
            self.damage.computeMedrasBreaks(recalculateDamages)
        for n, breakList in enumerate(self.damage.medrasBreaks):
            scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
            misrepairedPairs, repairEvents, remBreaks = calcMR.singleRepair(breakList.copy(), scaledSigma, finalTime)
            times = np.concatenate((np.linspace(0, 1, 13), np.linspace(1.25, 4, 12), np.linspace(5, 24, 20)))
            misrepairedBreaks = []
            for m in misrepairedPairs:
                misrepairedBreaks.append(m[0])
                misrepairedBreaks.append(m[1])
            for it in range(1, len(times)):
                eventsInThisTime = []
                for e in repairEvents:
                    timeEv = e[0]
                    if timeEv >= times[it - 1] and timeEv < times[it]:
                        eventsInThisTime.append(e)
                for e in eventsInThisTime:
                    endOne = e[1]
                    endTwo = e[2]
                    breakOne = self.damage.medrasBreaks[n][endOne]
                    breakTwo = self.damage.medrasBreaks[n][endTwo]
                    if breakOne not in misrepairedBreaks or breakTwo not in misrepairedBreaks:
                        pos1 = np.array(breakOne[1])
                        pos2 = np.array(breakTwo[1])
                        pos = np.array((pos1 + pos2) / 2)
                        newdsbpos = []
                        for p in self.damage.DSBPositions:
                            if not np.array_equal(p, pos):
                                print(p, pos)
                                newdsbpos.append(p)
                        print(times[it-1], times[it], e, 'pre', len(self.damage.DSBPositions))
                        self.damageDSBPositions = newdsbpos
                        print(times[it-1], times[it], e, 'post', len(self.damage.DSBPositions))

                    #print(times[it-1], times[it], e, breakOne, breakTwo)
                    # endone and endtwo are indexes in breaklist, which can be accessed as self.damage.medrasbreaks[n]
                    # misrepairedPairs include ends that will not ever be repaired. can be accessed as misrepairedPairs[i][0] and misrepairedPairs[i][1] (these are the same as self.damage.medrasbreaks[n])
                    # remember a break in medras format is given by 'index', 'position', 'comlpexitiy', 'chromID', damagecrhompos... (not relevant from here)

        # idea: for each time in repair events, get the repaired dsb and remove them from the
        # list of dsbs in damagetodna class. Create then a 3D picture using the methods there and
        # eventually a video. To do this, damage.DNAPositions needs to be changed. This is an array
        # with the mean position of the DSB. We can calculate that using something like
        # posSB1 = np.array(self.damageMap[iCh][closestPosInStrand1][2].position)
        # posSB2 = np.array(self.damageMap[iCh][iBp+i2][3].position)
        # self.DSBPositions.append((posSB1 + posSB2)/2)

    def misrepairSpectrum(self):
        calcMR = MisrepairCalculator()
        aberrations = Aberrations()
        scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
        # Test to see if chromosome ID data is available, abort if not
        firstBreak = self.damage.medrasBreaks[0][0]
        if firstBreak[3][0] == -1:
            print('No chromosome IDs found in data file, aborting!')
            return None
        noChroms = self.damage.Chromosomesizes[0]
        baseChromosomes = [[n, 0, self.damage.Chromosomesizes][1:][n] for n in range(noChroms)]
        fullFrags = []
        allChroms = []
        allRings = []
        for m, breakList in enumerate(self.damage.medrasBreaks):
            if m >= self.maxExposures:
                break
            misrepList, repList, remBreaks = calcMR.singleRepair(copy.deepcopy(breakList), None, scaledSigma,
                                                                 finalTime=self.simulationLimit)
            trimMisrep, trimRemBreaks = self.prepareDamage(misrepList, remBreaks, baseChromosomes)
            chroms, rings, frags = aberrations.doRepair(baseChromosomes, trimMisrep, remBreaks=trimRemBreaks, index=m,
                                                        breaks=len(breakList)//2, baseBreaks=breakList, plot=self.doPlot,
                                                        allFragments=self.allFragments, outFile='plot'+str(m)+'.png')
            frags = [f + [m] for f in frags]
            fullFrags += frags
            allChroms += chroms
            allRings += rings
        if self.listAcentrics:
            self.listAcentricSizes(aberrations, baseChromosomes, allChroms + allRings)
        return None

    def misrepairSeparation(self):
        calcMR = MisrepairCalculator()
        # Initialize some shared values
        bins = 500
        scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
        radius = scaledSigma / self.sigma
        maxSeparation = 2 * radius
        rBins = [n*(maxSeparation*1.0/bins) for n in range(bins+1)]
        if self.separationRun is False:
            self.separationRun = True
            print('BreakSets\tEmptySets\tSeparation (um):\t', end=' ')
            print('\t'.join(map(str, rBins)))
        print(len(self.damage.medrasBreaks), '\t', self.damage.emptySets, '\t', end=' ')
        misrepairSeps = []
        totBreaks = 0
        for m, breakList in enumerate(self.damage.medrasBreaks):
            if m >= self.maxExposures:
                break
            totBreaks += len(breakList)
            for n in range(self.repeats):
                misrepList, repList, remBreaks = calcMR.singleRepair(copy.deepcopy(breakList), None, scaledSigma)
                misrepairSeps += [mis[2] for mis in misrepList]
        if len(misrepairSeps) > 0:
            print(totBreaks / m / 2, '\t', len(misrepairSeps) / (m * self.repeats), end='\t')
            print('\t'.join(map(str, np.histogram(misrepairSeps, rBins, density=True)[0])))
        else:
            print()
        return None

    def dsbSeparation(self):
        # Initalize some shared values
        bins = 500
        scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
        radius = scaledSigma / self.sigma
        maxSeparation = 2 * radius
        rBins = [n*(maxSeparation*1.0/bins) for n in range(bins+1)]
        if self.separationRun is False:
            self.separationRun = True
            print('BreakSets\tEmptySets\tSeparation (um):\t', end=' ')
            print('\t'.join(map(str, rBins)))
        print(len(self.damage.medrasBreaks), '\t', self.damage.emptySets, '\t', end=' ')
        seps = []
        for m, breakList in enumerate(self.damage.medrasBreaks):
            for i in range(0, len(breakList), 2):
                for j in range(i+2, len(breakList), 2):
                    seps.append(np.sqrt(MisrepairCalculator.distanceToSq((breakList[i][1], breakList[j][1]))))
        if len(seps) > 0:
            print('\t'.join(map(str, np.histogram(seps, rBins, density=True)[0])))
        else:
            print()
        return None

    def radialDSBs(self):
        # Initialize some shared values
        bins = 1600
        maxSeparation = 4.0
        rBins = [n * (maxSeparation * 1.0/bins) for n in range(bins + 1)]
        # For first run, print header
        if self.radialRun is False:
            self.radialRun = True
            print('BreakSets\tEmptySets\tSeparation(um):\t', end=' ')
            print('\t'.join(map(str,rBins)))
        print(len(self.damage.medrasBreaks), '\t', self.damage.emptySets, '\t', end=' ')
        seps = []
        for m, breakList in enumerate(self.damage.medrasBreaks):
            for i in range(len(breakList)):
                seps.append(np.sqrt(pow(breakList[i][1][0], 2) + pow(breakList[i][1][1], 2)))
        if len(seps) > 0:
            print('\t'.join(map(str, np.histogram(seps, rBins, density=True)[0])))
        else:
            print()
        return None

    def trackBreaks(self):
        # Extract file data
        scaledSigma = self.sigma * stats.gmean(self.damage.ScoringVolume[1:4])
        # For first run, print header
        if self.radialRun is False:
            self.radialRun = True
            print('Break count\tBreaks')
        totalTracks = 0
        totalBreaks = 0
        lastTrack = -1
        for m, breakList in enumerate(self.damage.medrasBreaks):
            breakCount = 0
            for i in range(len(breakList)):
                if breakList[i][0] != lastTrack:
                    if lastTrack != 1:
                        print(lastTrack, breakCount)
                        totalTracks += 1
                        totalBreaks += breakCount
                    lastTrack = breakList[i][0]
                    breakCount = 1
                else:
                    breakCount += 1
            if lastTrack != -1:
                print(lastTrack, breakCount)
        retString = '\t'.join(map(str, [totalTracks, totalBreaks, totalBreaks/totalTracks]))
        return retString

    def prepareDamage(self, misrepairList, remainingBreaks, chromosomes):
        # Trim misrepair list - store chromosome, genetic position, and orientation of each break (a&b)
        trimMisrep = [[[a[3][1], a[4], a[5], a[1]], [b[3][1], b[4], b[5], a[1]]] for a,b,c,d in misrepairList]
        # Extract same data for remaining breaks
        trimBreaks = [[a[3][1], a[4], a[5], a[1]] for a in remainingBreaks]
        # Scale by chromosome sizes, if needed
        for misrep in trimMisrep:
            for damage in misrep:
                chromID = damage[0]
                if damage[1] < 1:
                    damage[1] *= chromosomes[chromID][2]
        for damage in trimBreaks:
            chromID = damage[0]
            if damage[1] < 1:
                damage[1] *= chromosomes[chromID][2]
        # Filter out duplicated DNA ends
        for n in range(len(trimBreaks)-1, 0, -1):
            if trimBreaks[n][0] == trimBreaks[n-1][0] and trimBreaks[n][1] == trimBreaks[n-1][1]:
                trimBreaks.pop(n)
        return trimMisrep, trimBreaks

    def listAcentricSizes(self, aberrations, baseChromosomes, finalChromosomes, pos=0.5):
        fullChroms = aberrations.characterizeChroms(finalChromosomes)
        print('\n\nAcentric sizes:')
        for c in fullChroms:
            centromereCount = aberrations.centricCount(c, baseChromosomes, pos)
            if centromereCount == 0:
                theSize = sum(abs(f[1]-f[2]) for f in c[3])
                print(theSize)

    def summarizeKinetics(self, repairTimes):
        repairTimes.sort()
        timeStep = 0.1
        currTime = 0
        kinetic = [1]
        for n, t in enumerate(repairTimes):
            while t > currTime + timeStep:
                kinetic.append(1.0 - 1.0 * n / len(repairTimes))
                currTime += timeStep
                if currTime + timeStep > self.kineticLimit:
                    break
            if currTime + timeStep > self.kineticLimit:
                break
        kinetic = kinetic + [0]*(int(self.kineticLimit / timeStep) - len(kinetic))
        return '\t'.join(map(str, kinetic))

    def summarizeFidelity(self, complexity, outputs):
        totalBreaks = 1.0 * sum(o[0] for o in outputs)
        averageBreaks = np.mean([o[0] for o in outputs])
        breakStdev = np.std([o[0] for o in outputs])
        averageMisrep = np.mean([o[0]*o[1] for o in outputs])/averageBreaks
        misrepStdev = np.std([o[1] for o in outputs])
        fileAverages = [sum([o[n]*o[0] for o in outputs])/totalBreaks for n in range(3,6)]
        smry = '\tSummary\t' + str(totalBreaks) + '\t' + str(complexity) + '\t' + str(averageBreaks) + \
               '\t' + str(breakStdev) + '\t' + str(averageMisrep) + '\t' + str(misrepStdev) + '\t' + str(fileAverages[0])
        return smry

    def __sortNicely(self, l):
        convert = lambda text: int(text) if text.isdigit() else text
        alphanumkey = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
        l.sort(key=alphanumkey)

class MisrepairCalculator:
    def __init__(self, fastRate=2.07, slowRate=0.259, fastFoci=8.12, slowFoci=0.405, mmejFoci=0.0176):
        # Input parameters, physical DSB interaction rates based on break complexity (per hour)
        self.fastRate = fastRate
        self.slowRate = slowRate
        # Foci clearance rates (per hour)
        self.fastFoci = fastFoci
        self.slowFoci = slowFoci
        self.mmejFoci = mmejFoci
        self.rateTable = None

    def fullRepair(self, baseBreaks, sigma, repeats=1, addFociClearance=True, radius=1,
                   chromSizes=None, sizeLimit=0, finalTime=np.inf):
        # Sort breaks in order of time of creation and build table of interaction rates
        baseBreaks.sort(key=lambda x:x[7])
        self.buildRateTable(baseBreaks, sigma)
        fullMisrepairPairs = []
        fullRepairEvents = []
        for n in range(repeats):
            pairRates = []
            breakList = copy.deepcopy(baseBreaks)
            misrepairedPairs, repairEvents, remBreaks = self.singleRepair(breakList, sigma, finalTime)
            # Filter out any intra-chromosome breaks below size limit
            if sizeLimit >= 0 and chromSizes is not None:
                for i in range(len(misrepairedPairs)-1, -1, -1):
                    misrepair = misrepairedPairs[i]
                    if misrepair[0][3][1] == misrepair[1][3][1]:
                        chromID = misrepair[0][3][1]
                        misrepSize = chromSizes[chromID] * abs(misrepair[0][4]-misrepair[1][4])
                        if misrepSize < sizeLimit:
                            misrepairedPairs.pop(i)
            fullMisrepairPairs.append(misrepairedPairs)
            fullRepairEvents += repairEvents
        if addFociClearance:
            fullRepairTimes = []
            for t, p1, p2, complexity in fullRepairEvents:
                repRate = self.fastFoci
                if baseBreaks[p1][0] != baseBreaks[p2][0]:
                    repRate = self.mmejFoci
                else:
                    if complexity > 0:
                        repRate = self.slowFoci
                fullRepairTimes.append(t - np.log(np.random.random()) / repRate)
            fullRepairTimes = sorted(fullRepairTimes)
        else:
            fullRepairTimes = sorted([x[0] for x in fullRepairEvents])
        # Calculate rate of misrepair and stdev
        misrepairCounts = [len(m) for m in fullMisrepairPairs]
        misrepRate = np.mean(misrepairCounts) / (0.5 * len(baseBreaks))
        stdevRate = np.std(misrepairCounts) / (0.5 * len(baseBreaks))
        # Calculate inter-chromosome rate
        interChromEvents = sum(sum(misrep[3] for misrep in repeat) for repeat in fullMisrepairPairs)
        interChromRate = interChromEvents / max(1.0, 1.0 * sum(misrepairCounts))
        analyticRate, errRate = self.analyticRepair(baseBreaks, sigma, radius)
        return misrepRate, stdevRate, interChromRate, fullRepairTimes, analyticRate, errRate

    def singleRepair(self, breakList, sigma=None, finalTime=np.inf):
        # Sort breaks in order of creation in time and set up interaction rates if needed
        if self.rateTable is None:
            breakList.sort(key=lambda x:x[7])
            self.buildRateTable(breakList, sigma)
        rateTable = self.rateTable.copy()
        # Get fast/slow kinetic data from breaklist
        repairRate = np.array([self.fastRate/2 if b[2]==False else self.slowRate/2 for b in breakList])
        # Sample interaction time for every break
        interactionSamples = -np.log(np.random.rand(len(breakList)))
        # Get base time, and initialize lists of live and pending breaks
        baseTime = breakList[0][7]
        liveBreaks = [n for n,b in enumerate(breakList) if b[7] <= baseTime]
        pendingBreaks = [n for n,b in enumerate(breakList) if b[7] > baseTime]
        lastBreak = liveBreaks[-1]
        if len(pendingBreaks) > 0:
            nextBreakTime = breakList[pendingBreaks[0]][7]
        else:
            nextBreakTime = np.inf
        # Reporting variables
        repairEvents = []
        misrepairedPairs = []
        # Simulate repair until no more breaks reman, or we exceed the simulation limit
        while len(liveBreaks) + len(pendingBreaks) > 0:
            if len(liveBreaks) > 0:
                # For every break, baseline rate is repair rate divided by total interaction with all DSB ends
                rateSums = np.sum(rateTable[:, 0:lastBreak+1], axis=1)
                interactionTimes = np.array(interactionSamples)[liveBreaks]/(repairRate[liveBreaks]*rateSums[liveBreaks])
                nextTime = baseTime+min(interactionTimes)
            else:
                nextTime = nextBreakTime
            # Break loop if next event is after the simulation end
            if min(nextTime, nextBreakTime) > finalTime: break
            # If next repair is before next break, log the repair and remove the break ends
            if nextTime < nextBreakTime:
                # Get repaired break ends based on matching repair time, then select a pair
                endOne = liveBreaks[np.argmin(interactionTimes)]
                endTwo = self.pickRepair(rateTable[endOne, 0:lastBreak+1], endOne)
                # Assign complexity for foci clearance based on most complex end
                complexity = max(breakList[endOne][2], breakList[endTwo][2])
                # Identify type of repair, and log details
                if breakList[endOne][0] != breakList[endTwo][0]:
                    # If ends are not from the same break, then it is a misrepair. Log details.
                    if breakList[endOne][3][1] == breakList[endTwo][3][1] and breakList[endOne][3][2] == breakList[endTwo][3][2]:
                        interChrom = 0
                    else:
                        interChrom = 1
                    separation = np.sqrt(self.distanceToSq(breakList[endOne][1], breakList[endTwo][1]))
                    misrepairedPairs.append([breakList[endOne], breakList[endTwo], separation, interChrom])
                    # Track if either end is the second end from a DSB. Log reduction in breaks if so
                    p1Partner = endOne + 1 - 2 * (endOne % 2)
                    if breakList[p1Partner] == 0:
                        repairEvents.append([nextTime, endOne, endTwo, complexity])
                    p2Partner = endTwo + 1 - 2 * (endTwo % 2)
                    if breakList[p2Partner] == 0:
                        repairEvents.append([nextTime, endOne, endTwo, complexity])
                else:
                    # Matched DSB ends always clear 1 DSB, no misrepair
                    repairEvents.append([nextTime, endOne, endTwo, complexity])
                # Tidy up break data
                liveBreaks.pop(liveBreaks.index(endOne))
                liveBreaks.pop(liveBreaks.index(endTwo))
                rateTable[:, endOne] = 0
                rateTable[:, endTwo] = 0
                breakList[endOne] = 0
                breakList[endTwo] = 0
            else:
                # If next DSB is before next repair, add pending DSBs and update interaction
                baseTime = nextBreakTime
                liveBreaks.append(pendingBreaks.pop(0))
                while len(pendingBreaks) > 0 and baseTime >= breakList[pendingBreaks[0]][7]:
                    liveBreaks.append(pendingBreaks.pop(0))
                interactionSamples = -np.log(np.random.rand(len(breakList)))
                if len(pendingBreaks) > 0:
                    nextBreakTime = breakList[pendingBreaks[0]][7]
                else:
                    nextBreakTime = np.inf
                lastBreak = liveBreaks[-1]
        # Return a list of unrepaired breaks, if we halted before all were repaired
        if finalTime <= np.inf:
            remBreaks = [b for b in breakList if b != 0]
            return misrepairedPairs, repairEvents, remBreaks
        # Otherwise, clean up any breaks which were missed and return 'full' repair
        # This is an effectively random pairing, but should hopefully just be tidying up corner cases
        if len(liveBreaks) + len(pendingBreaks) > 0:
            remBreaks = [p for p in range(len(breakList)) if breakList[p] != 0]
            while len(remBreaks) > 0:
                p1 = remBreaks.pop()
                p2 = remBreaks.pop()
                repairEvents.append([1E6, p1, p2, 1])
                misrepairedPairs.append([breakList[p1], breakList[p2], -1, 0])
        # Otherwise, just return 'full' repair
        return misrepairedPairs, repairEvents, []

    def pickRepair(self, interactionArray, n):
        cumulativeInteraction = np.cumsum(interactionArray)
        interactionSample = np.random.uniform() * cumulativeInteraction[-1]
        offsetInteraction = cumulativeInteraction - interactionSample
        chosenInteraction = np.where(offsetInteraction >= 0)[0][0]
        return chosenInteraction

    def analyticRepair(self, breakList, sigma, radius):
        if self.rateTable is None:
            breakList.sort(key=lambda x:x[7])
            self.buildRateTable(breakList, sigma)
        rateTable = self.rateTable.copy()
        correctRate = 0.0
        errRate = 0.0
        for i in range(len(breakList)):
            for j in range(i + 1, len(breakList)):
                if breakList[i][0] != breakList[j][0]:
                    errRate += rateTable[i][j]
        errRate = errRate / (len(breakList) * 0.5)
        base = 0.81026
        rate = 8.51
        skewCorrection = base + (1 - base) * (1 - np.exp(-rate * sigma/radius))
        errRate = errRate / skewCorrection
        if errRate < 1E-10: return 0, errRate
        analyticMisrep = 1 - 2 * (np.arctan(errRate + 1) - np.arctan(1)) / errRate
        return analyticMisrep, errRate

    def buildRateTable(self, baseBreaks, sigma):
        sigmaSq = sigma * sigma
        positions = np.array([b[1] for b in baseBreaks])
        seps = spatial.distance.pdist(positions)
        seps = spatial.distance.squareform(seps)
        self.rateTable = np.exp(-(seps * seps) / sigmaSq)
        self.rateTable = np.clip(self.rateTable, 1E-200, None)
        np.fill_diagonal(self.rateTable, 0)

    def interactionRate(self, p1, p2, sigma):
        return np.exp(self.__logInteractionRate(p1, p2, sigma))

    def __logInteractionRate(self, p1, p2, sigma):
        return -self.distanceToSq(p1, p2) / (2 * sigma * sigma)

    @staticmethod
    def distanceToSq(p1, p2):
        dx = p1[0] - p2[0]
        dy = p1[0] - p2[0]
        dz = p1[0] - p2[0]
        return dx*dx+dy*dy+dz*dz

class Aberrations:
    def __init__(self, minFragment=0, maxFragment=5.7, largeMisrepThreshold=3000000):
        '''
        :param minFragment: MBP, minimum measurable fragment for DNA loss
        :param maxFragment: MBP, maximum measurable fragment for DNA loss
        :param largeMisrepThreshold: BP, size of 'large' misrepair
        '''
        self.minFragment = minFragment
        self.maxFragment = maxFragment
        self.largeMisrepThreshold = largeMisrepThreshold

    def doRepair(self, chromosomes, repairs, remBreaks=None, index=0, breaks=-1, baseBreaks=None, plot=False,
                 allFragments=False, inFile=None, outFile=None):
        # Build our breaklist, and abort if empty
        breakList = [b for r in repairs for b in r]
        if remBreaks is not None:
            breakList += remBreaks
        if len(breakList) == 0:
            print(index, '\tNo misrepair!')
            return [], [], []
        # Split up chromosomes, with repaired and remaining breaks
        breakList.sort()
        chromFrags = self.splitChromosomes(chromosomes, breakList)
        chromList = copy.deepcopy(chromFrags)
        rings = []
        # Rejoin them
        for b1, b2 in repairs:
            c1, p1, d1, l1 = b1
            c2, p2, d2, l2 = b2
            n1, end1 = self.indexChrom(chromList, c1, p1, d1)
            n2, end2 = self.indexChrom(chromList, c2, p2, d2)
            # If ends are both the same chromosome, it's a ring
            if n1 == n2:
                chromOne = chromList.pop(n1)
                rings.append(chromOne)
                continue
            # If not, stick them together and append new fragment
            if n1 > n2:
                chromOne = chromList.pop(n1)
                chromTwo = chromList.pop(n2)
            else:
                chromOne = chromList.pop(n2)
                chromTwo = chromList.pop(n1)
            newFrag = self.appendFragments(chromOne, chromTwo, end1, end2)
            chromList.append(newFrag)
        linearChromosomes = self.characterizeChroms(chromList)
        ringChromosomes = self.characterizeChroms(rings, True)
        # Plot and print stats if requested
        if plot:
            self.drawChroms(chromosomes, linearChromosomes, ringChromosomes, inFile=inFile, outFile=outFile)
        print(index, '\t', breaks, '\t', len(remBreaks), '\t', self.misrepairStats(repairs, chromosomes),
              '\t', self.calculateComplexities(chromList + rings), '\t', self.centricCheck(chromosomes, linearChromosomes),
              '\t', self.centricCheck(chromosomes, ringChromosomes), end='')
        if baseBreaks is not None:
            lostDNA, lostFragments = self.dnaLoss(chromosomes, linearChromosomes + ringChromosomes)
            print('\t', self.fragmentDistribution(chromosomes, baseBreaks), '\t', lostDNA, end='')
            if allFragments:
                print('\t\t', '\t'.join(map(str, sorted(lostFragments))), end='')
        print()
        return chromList, rings, lostFragments

    def splitChromosomes(self, chromosomes, breaks):
        chromStack = copy.deepcopy(chromosomes)
        chromStack = [[c + [None, None]] for c in chromStack]
        for b in breaks:
            breakChrom = b[0]
            breakPos = b[1]
            spatialPos = b[3]
            for n, chrom in enumerate(chromStack):
                if chrom[0][0] == breakChrom and breakPos > chrom[0][1] and breakPos < chrom[0][2]:
                    chromStack.pop(n)
                    topChrom = [chrom[0][0], chrom[0][1], breakPos, chrom[0][3], spatialPos]
                    chromStack.append([topChrom])
                    botChrom = [chrom[0][0], breakPos, chrom[0][2], spatialPos, chrom[0][4]]
                    chromStack.append([botChrom])
        return chromStack

    def appendFragments(self, chromOne, chromTwo, end1, end2):
        # Need to flip one chromosome if it is RHS-RHS or LHS-LHS joining
        if end1 == end2:
            chromTwo = [[c[0], c[2], c[1], c[4], c[3]] for c in reversed(chromTwo)]
        # Line up appropriately
        if end1 == 0:
            return chromTwo + chromOne
        else:
            return chromOne + chromTwo

    def fragmentDistribution(self, baseChromosomes, breakList):
        breakPoints = []
        for damage in breakList[::2]:
            breakPoints.append([damage[3][1], damage[4] * baseChromosomes[damage[3][1]][2]])
        # Sort by chromosome and genetic position
        breakPoints.sort(key=lambda x:(x[0], x[1]))
        currChrom = breakPoints[0][0]
        currPos = 0
        totalFrags = 0
        for damage in breakPoints:
            if currChrom != damage[0]:
                fullChomSize = baseChromosomes[currChrom][2]
                remLength = fullChomSize - currPos
                if remLength < self.maxFragment and remLength > self.minFragment:
                    totalFrags += remLength
                currPos = 0
                currChrom = damage[0]
            fragLength = damage[1] - currPos
            if fragLength < self.maxFragment and fragLength > self.minFragment:
                totalFrags += fragLength
            currPos = damage[1]
        return totalFrags

    def dnaLoss(self, baseChromosomes, finalChromosomes):
        lostDNA = 0
        allFragments = []
        for c in finalChromosomes:
            length = c[0]
            spatialSep = c[-1]
            if length < self.maxFragment:
                lostDNA += length
            allFragments.append([length, spatialSep, c[1]])
        return lostDNA, allFragments

    def calculateComplexities(self, finalFragments):
        simple = 0
        complexes = 0
        for f in finalFragments:
            junctions = 0
            for n in range(1, len(f)):
                if (f[n - 1][0] // 2 != f[n][0] // 2):
                    junctions += 1
            if junctions == 1:
                simple += 1
            else:
                if junctions > 1:
                    complexes += 1
        return str(simple) + "\t" + str(complexes)

    def misrepairStats(self, misreps, chromosomes):
        interChrom = 0
        largeMisrep = 0
        for misrep in misreps:
            if misrep[0][0] != misrep[1][0]:
                interChrom += 1
            else:
                chromID = misrep[0][0]
                chromLength = chromosomes[chromID][2]
                if abs((misrep[0][1] - misrep[1][1])) >= self.largeMisrepThreshold:
                    largeMisrep += 1
        return str(len(misreps)) + "\t" + str(largeMisrep) + "\t" + str(interChrom)

    def indexChrom(self, chromFragments, c, p, d):
        # Tolerate match within 1 part in 10^12, in case of rounding issues
        for n, chrom in enumerate(chromFragments):
            firstChrom = chrom[0]
            if firstChrom[0] == c:
                if abs(firstChrom[1] - p) < 1E-12 * p:
                    if ((firstChrom[2] - firstChrom[1] > 0 and d > 0) or (firstChrom[2] - firstChrom[1] < 0 and d < 0)):
                        return n, 0
            lastChrom = chrom[-1]
            if lastChrom[0] == c:
                if abs(lastChrom[2] - p) < 1E-12 * p:
                    if ((lastChrom[2] - lastChrom[1] > 0 and d > 0) or (lastChrom[2] - lastChrom[1] < 0 and d > 0)):
                        return n, 1
            print('Failed to find matching chromosome fragment! Requested chromosome details:')
            print(c, p, d)
            print('All Fragments:')
            for c in chromFragments:
                print(c)
            return None

    def characterizeChroms(self, chroms, rings=False, doPrint=False):
        retChroms = []
        for c in chroms:
            totLen = 0
            lenDict = dict()
            for f in c:
                totLen += abs(f[2] - f[1])
                if f[0] not in lenDict:
                    lenDict[f[0]] = abs(f[2] - f[1])
                else:
                    lenDict[f[0]] += abs(f[2] - f[1])
                if doPrint:
                    print(f[0], abs(f[2]- f[1]), '\t', end=' ')
            if c[0][3] is None or c[-1][4] is None:
                spatialLen = -1
            else:
                spatialLen = np.linalg.norm(np.array(c[0][3]) - np.array(c[-1][4]))
            if doPrint:
                print('\t', totLen)
            mainChrom = max(iter(lenDict.keys()), key=(lambda key: lenDict[key]))
            retChroms.append([totLen, mainChrom, rings, c, spatialLen])
        return retChroms

    def centricCheck(self, baseChromosomes, finalChromosomes, pos=0.5):
        normal = 0
        acentric = 0
        multicentric = 0
        largeLoss = 0
        for c in finalChromosomes:
            centromereCount = self.centricCount(c, baseChromosomes, pos)
            if centromereCount == 0:
                acentric += 1
                if sum(abs(f[1] - f[2]) for f in c[3]) > self.largeMisrepThreshold:
                    largeLoss += 1
            if centromereCount == 1:
                normal += 1
            if centromereCount > 1:
                multicentric += 1
        return '\t'.join(map(str, [normal, acentric, largeLoss, multicentric]))

    def centricCount(self, c, baseChromosomes, pos=0.5):
        centromereCount = 0
        for f in c[3]:
            start = min(f[1], f[2])
            end = max(f[1], f[2])
            cent = pos * baseChromosomes[f[0]][2]
            if cent >= start and cent <= end:
                centromereCount += 1
        return centromereCount

    def drawChroms(self, baseChromosomes,chromosomes,rings, inFile=None, outFile=None):
        chromSet = chromosomes + rings
        chromSet.sort(key=lambda x: x[1])
        plt.style.use('dark_background')
        ax = plt.axes()
        ax.set_facecolor('black')
        width = 10
        intraOffset = 5
        interOffset = 30
        currX = interOffset
        maxHeight = max(c[0] for c in chromSet)
        yOffset = maxHeight * 1.1 + 10
        refY = -yOffset / 2.0
        lastChrom = -1
        rowStart = -100
        if inFile is not None:
            ax.text(currX, refY - 0.25 * yOffset, inFile, color=(0.5, 0.5, 0.5))
        for c in chromSet:
            if lastChrom // 2 == c[1] // 2:
                currX += intraOffset
            else:
                lastChrom = c[1]
                currX += interOffset
                if c[1] - rowStart > 10:
                    refY -= yOffset
                    currX = interOffset
                    rowStart = c[1]
            currY = refY
            # Draw a linear chromosome
            if c[2] == False:
                ax.text(currX + width / 2.0, currY - yOffset * 0.25,
                        str((c[1] // 2) + 1) + self.centromereTag(c, baseChromosomes), ha="center",
                        color=(0.5, 0.5, 0.5))
                for f in c[3]:
                    fragHeight = abs(f[2] - f[1]) / 2.0
                    rect = matplotlib.patches.Rectangle((currX, currY), width, fragHeight,
                                                        color=self.fetchColor(f[0]))
                    ax.add_patch(rect)
                    currY += fragHeight
                currX += width
            else:
                # Draw a ring chromosome
                totLength = c[0] / 2.0
                circRad = np.sqrt(totLength * width / 3.14159)
                currX += max(circRad, width / 2.0)
                ax.text(currX, currY - yOffset * .25,
                        str((c[1] // 2) + 1) + self.centromereTag(c, baseChromosomes), ha="center",
                        color=(0.5, 0.5, 0.5), fontweight="bold")
                currTheta = 0
                for f in c[3]:
                    fracLength = abs(f[2] - f[1]) / 2.0
                    deltaTheta = 360 * (fracLength / totLength)
                    wedge = matplotlib.patches.Wedge((currX, currY + circRad), circRad, currTheta,
                                                     currTheta + deltaTheta, color=self.fetchColor(f[0]))
                    ax.add_patch(wedge)
                    currTheta = currTheta + deltaTheta
                currX += max(circRad, width / 2.0)
        ax.autoscale()
        ax.axis('off')
        if outFile is None:
            plt.savefig("ChromAberrs.png")
        else:
            plt.savefig(outFile)
        if haveDisplay:
            plt.show()
        else:
            plt.close()

    def fetchColor(self, c):
        colorSet = [(0.847058824, 0.882352941, 0.305882353), (0.749019608, 0.007843137, 0.211764706),
                    (0.517647059, 0.682352941, 0.839215686),
                    (0.019607843, 0.670588235, 0.28627451), (0.68627451, 0.262745098, 0.57254902),
                    (0.88627451, 0.474509804, 0.145098039),
                    (0.901960784, 0.639215686, 0.184313725), (0.48627451, 0.431372549, 0.596078431),
                    (0.552941176, 0.847058824, 0.850980392),
                    (0.607843137, 0.701960784, 0.278431373), (0.168627451, 0.654901961, 0.545098039),
                    (0.854901961, 0.17254902, 0.223529412),
                    (0.380392157, 0.807843137, 0.858823529), (0.949019608, 0.08627451, 0.184313725),
                    (0.588235294, 0.447058824, 0.666666667),
                    (0.996078431, 0.988235294, 0.498039216), (0.270588235, 0.650980392, 0.831372549),
                    (0.964705882, 0.937254902, 0.945098039),
                    (0.168627451, 0.698039216, 0.298039216), (0.301960784, 0.435294118, 0.674509804),
                    (0.894117647, 0.921568627, 0.498039216),
                    (0.968627451, 0.968627451, 0.576470588), (0.462745098, 0.788235294, 0.301960784)]
        return colorSet[c // 2]

    def centromereTag(self, chromosome, baseChromosomes, pos=0.5):
        centromereCount = self.centricCount(chromosome, baseChromosomes, pos)
        if centromereCount == 1:
            return ""
        if centromereCount == 0:
            return "#"
        return "*"