#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/8/22 2:01 PM

@author: alejandrobertolet
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class AverageTimeCurveOverRuns:
    def __init__(self, nRuns=0, refdata=None):
        self._nRuns = nRuns
        self.runlist = []
        self.refdata = refdata

    def AddTimeCurveForSingleRun(self, tc):
        self.runlist.append(tc)

    def DoStatistics(self, scaledToInitialValue=True):
        self.times = self.runlist[0].times
        self.avgyvalues = np.zeros(len(self.times))
        self.varyvalues = np.zeros(len(self.times))
        for j in range(len(self.times)):
            yvaluesfortimej = np.array([])
            for tc in self.runlist:
                if tc.yvalues[0] != 0:
                    if scaledToInitialValue:
                        yvaluesfortimej = np.append(yvaluesfortimej, tc.yvalues[j]/tc.yvalues[0])
                    else:
                        yvaluesfortimej = np.append(yvaluesfortimej, tc.yvalues[j])
            self.avgyvalues[j] = np.mean(yvaluesfortimej)
            self.varyvalues[j] = np.var(yvaluesfortimej)

    def Plot(self, fsize=(10,6)):
        if self.runlist[0].timeunit == 'h':
            xlabel = 'Time (h)'
        elif self.timeunit == 'min':
            xlabel = 'Time (min)'
        else:
            xlabel = 'Time (s)'
        ylabel = self.runlist[0].ylabel
        fig = plt.figure(figsize=fsize, tight_layout=True)
        ax = fig.add_subplot(111)
        ax.plot(self.times, self.avgyvalues, label='Model')
        ax.fill_between(self.times, self.avgyvalues-np.sqrt(self.varyvalues), self.avgyvalues+np.sqrt(self.varyvalues), alpha=0.5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xlim([0, np.max(self.times)*1.01])
        ax.set_ylim([0, np.max(self.avgyvalues + np.sqrt(self.varyvalues))])
        plt.show()
        if self.refdata is not None:
            ax.scatter(self.refdata.x, self.refdata.y, label='Experimental points')
            ax.legend()
        ax.grid()
        return fig, ax

    def WriteCSV(self, resultsFolder='./'):
        df = pd.DataFrame({"Time" : self.times, "Avg-" + self.runlist[0].ylabel : self.avgyvalues, "Var-" + self.runlist[0].ylabel : self.varyvalues})
        df.to_csv(resultsFolder + self.runlist[0].ylabel + '.csv', index=False)

    @property
    def NRuns(self):
        if len(self.runlist) > 0:
            return len(self.runlist)
        else:
            return self._nRuns
    @NRuns.setter
    def NRuns(self, n):
        self._nRuns = n

class TimeCurveForSingleRun:
    def __init__(self, ylabel='', timeunit='h', refdata=None):
        self.ylabel = ylabel
        self.timeunit = timeunit
        self.refdata = refdata
        self.times = np.array([])
        self.yvalues = np.array([])

    def AddTimePoint(self, t, y):
        if self.timeunit == 'h':
            t = t / 3600
        elif self.timeunit == 'min':
            t = t / 60
        self.times = np.append(self.times, t)
        self.yvalues = np.append(self.yvalues, y)

    def Plot(self, fsize=(10,6), scaledToInitialValue=True):
        if self.timeunit == 'h':
            xlabel = 'Time (h)'
        elif self.timeunit == 'min':
            xlabel = 'Time (min)'
        else:
            xlabel = 'Time (s)'
        if scaledToInitialValue:
            if self.yvalues[0] != 0:
                self.yvalues = self.yvalues / self.yvalues[0]
        fig = plt.figure(figsize=fsize, tight_layout=True)
        ax = fig.add_subplot(111)
        ax.plot(self.times, self.yvalues, label='Model')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(self.ylabel)
        ax.set_xlim([0, np.max(self.times)*1.01])
        ax.set_ylim([0, np.max(self.yvalues)])
        plt.show()
        if self.refdata is not None:
            ax.scatter(self.refdata.x, self.refdata.y, label='Experimental points')
            ax.legend()
        ax.grid()
        return fig, ax

    def WriteCSV(self, resultsFolder='./', scaledToInitialValue=True):
        if scaledToInitialValue:
            df = pd.DataFrame({"Time" : self.times, self.ylabel : self.yvalues / self.yvalues[0], "Initial number" : self.yvalues[0]})
        else:
            df = pd.DataFrame({"Time": self.times, self.ylabel: self.yvalues})
        df.to_csv(resultsFolder + self.ylabel + '.csv', index=False)