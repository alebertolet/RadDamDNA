#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 8/28/22 7:59 PM

@author: alejandrobertolet
"""

from RadDamDNA.damage import *
import random
from scipy import interpolate
import pickle

picklefile = '/Users/ai925/source/workspace/repair/yieldData.pickle'
readfromdisk = False

fig = plt.figure()
fig.set_size_inches((12, 12))
ax = fig.add_subplot(111)

if not readfromdisk:
    nboot = 40

    maxdose = 6.0

    base = '/Users/ai925/Dropbox (Partners HealthCare)/Microdosimetry Project/ChemMicrodosimetry/nucleusSims/'
    #particles = ['xray', 'proton', 'proton', 'proton', 'proton', 'proton', 'proton', 'proton', 'alpha', 'alpha', 'alpha', 'alpha', 'alpha']
    #energies = ['250keV', '100MeV', '20MeV', '10MeV', '5MeV', '2MeV', '1MeV', '0.8MeV', '8MeV', '6MeV', '4MeV', '2MeV', '1MeV']
    #c = ['black', 'gray', 'slategray', 'skyblue', 'blue', 'green', 'olive', 'orange', 'peru', 'brown', 'salmon', 'red', 'darkred']

    particles = ['alpha', 'alpha', 'alpha', 'alpha', 'alpha']
    energies = ['8MeV', '6MeV', '4MeV', '2MeV', '1MeV']
    c = ['black', 'gray', 'red', 'green', 'blue']
    particles.reverse()
    energies.reverse()
    c.reverse()

    data = {}
    # Loop over all particles and energies specified above
    for ip, e in enumerate(energies):
        identifier = particles[ip] + '/' + e
        data[identifier] = {}
        basepath = base + particles[ip] + '/sims/' + e + '.txt/'
        nfiles = len(os.listdir(basepath))
        #Section to get what directories actually contains both dose and SDD. Disregard others!
        listOfAvailableDirs = []
        for j in range(nfiles):
            newpath = basepath + str(j) + '/'
            if str(j) in os.listdir(basepath):
                files = os.listdir(newpath)
                if 'DNADamage_sdd.txt' in files and 'DNADamage.phsp' in files:
                    if os.path.getsize(newpath + 'DNADamage_sdd.txt') > 0: #only those with actual data
                        listOfAvailableDirs.append(j)
        print(identifier, str(int(ip/len(energies)*100)) + '%...')
        # Preparing arrays
        Dose = np.linspace(0, maxdose, 100)
        DSB = np.zeros(Dose.shape)
        alldsbs = np.zeros([len(Dose), nboot])
        allssbs = np.zeros([len(Dose), nboot])
        allbds = np.zeros([len(Dose), nboot])
        SSB = np.zeros(Dose.shape)
        BD = np.zeros(Dose.shape)
        ntries = np.zeros(Dose.shape)
        # Bootstrapping
        for i in range(nboot):
            damage = DamageToDNA()
            # Getting new random sequence along the available directories each time
            neworder = random.sample(listOfAvailableDirs, len(listOfAvailableDirs))
            for j in neworder:
                path = basepath + str(j) + '/'
                damage.readSDDAndDose(path, version='1.0')

            damage.populateDamages(stopAtDose=maxdose)
            damage.computeStrandBreaks()
            #damage.printDamageCount()
            dose, dsb = damage.getDoseResponseCurve(plot=False, q='dsb')
            dose, ssb = damage.getDoseResponseCurve(plot=False, q='ssb')
            dose, bd = damage.getDoseResponseCurve(plot=False, q='bd')
            dose = np.array(dose)
            dsb = np.array(dsb)
            ssb = np.array(dsb)
            bd = np.array(dsb)
            dose = np.insert(dose, 0, 0)
            dsb = np.insert(dsb, 0, 0)
            ssb = np.insert(ssb, 0, 0)
            bd = np.insert(bd, 0, 0)
            fdsb = interpolate.interp1d(dose, dsb)
            fssb = interpolate.interp1d(dose, ssb)
            fbd = interpolate.interp1d(dose, bd)
            for id, d in enumerate(Dose):
                if d >= 0 and d < np.max(dose):
                    DSB[id] += fdsb(d)
                    alldsbs[id, i] = fdsb(d)
                    SSB[id] += fssb(d)
                    allssbs[id, i] = fssb(d)
                    BD[id] += fbd(d)
                    allbds[id, i] = fbd(d)
                    ntries[id] += 1

        vardsb = np.zeros(Dose.shape)
        varssb = np.zeros(Dose.shape)
        varbd = np.zeros(Dose.shape)
        for j in range(alldsbs.shape[0]):
            vardsb[j] = np.var(alldsbs[j, :])
            varssb[j] = np.var(allssbs[j, :])
            varbd[j] = np.var(allbds[j, :])

        stddsb = np.sqrt(vardsb)
        stdssb = np.sqrt(varssb)
        stdbd = np.sqrt(varbd)

        DSB = DSB / ntries
        SSB = SSB / ntries
        BD = BD / ntries

        data[identifier]['Dose'] = Dose
        data[identifier]['DSB'] = DSB
        data[identifier]['SSB'] = SSB
        data[identifier]['BD'] = BD
        data[identifier]['stdDSB'] = stddsb
        data[identifier]['stdSSB'] = stdssb
        data[identifier]['stdBD'] = stdbd

        ax.plot(Dose, DSB, '-', label=particles[ip]+'-'+e, color=c[ip])
        ax.fill_between(Dose, DSB-stddsb, DSB+stddsb, color=c[ip], alpha=0.15)
        ax.set_xlabel('Dose (Gy)')
        ax.set_ylabel('DSB')
        fig.tight_layout()

    with open(picklefile, 'wb') as handle:
        pickle.dump(data, handle)

def readFromPickle(file):
    with open(file, 'rb') as f:
        data = pickle.load(f)
    for identifier in data:
        Dose = data[identifier]['Dose']
        maxdose = np.max(Dose)
        DSB = data[identifier]['DSB']
        SSB = data[identifier]['SSB']
        BD = data[identifier]['BD']
        stddsb = data[identifier]['stdDSB']
        stdssb = data[identifier]['stdSSB']
        stdbd = data[identifier]['stdBD']
        ax.plot(Dose, DSB, '-', label=identifier)
        ax.fill_between(Dose, DSB-stddsb, DSB+stddsb, alpha=0.15)
        ax.set_xlabel('Dose (Gy)')
        ax.set_ylabel('DSB')
        fig.tight_layout()
    return maxdose

if readfromdisk:
    maxdose = readFromPickle(picklefile)


ax.grid()
ax.set_xlim(0, 1.01*maxdose)
ax.set_ylim(0, None)
ax.legend()
plt.show()