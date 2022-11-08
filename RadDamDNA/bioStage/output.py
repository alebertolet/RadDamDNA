#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/8/22 2:01 PM

@author: alejandrobertolet
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class TimeCurveForSingleRun:
    def __init__(self, ylabel='', timeunit='h', refdata=None):
        self.ylabel = ylabel
        self.timeunit = timeunit
        self.refdata = refdata
        self.times = np.array([])
        self.yvalues = np.array([])

    def AddTimePoint(self, t, y):
        self.times = np.append(self.times, t)
        self.yvalues = np.append(self.yvalues, y)

    def Plot(self, fsize=(10,6), scaledToInitialValue=True):
        if self.timeunit == 'h':
            self.times = self.times / 3600
            xlabel = 'Time (h)'
        elif self.timeunit == 'min':
            self.times = self.times / 60
            xlabel = 'Time (min)'
        else:
            xlabel = 'Time (s)'
        if scaledToInitialValue:
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

    def WriteCSV(self, resultsFolder='./'):
        df = pd.DataFrame({"Time" : self.times, self.ylabel : self.yvalues})
        df.to_csv(resultsFolder + self.ylabel + '.csv', index=False)