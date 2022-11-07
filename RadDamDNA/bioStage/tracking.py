#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 11/6/22 2:30 PM

@author: alejandrobertolet
"""

class BeStep:
    def __init__(self, initialPosition, initialTime = 0, complexity = 0):
        self.Position = initialPosition
        self.Time = initialTime
        self.Complexity = complexity
        self.Status = 'Damaged'
        self.RepairProteinStatus = 'NotPresent'

    @property
    def Time(self):
        return self._time
    @Time.setter
    def Time(self, t):
        self._time = t

    @property
    def Position(self):
        return self._position
    @Position.setter
    def Position(self, p):
        self._position = p

    @property
    def Complexity(self):
        return self._complexity
    @Complexity.setter
    def Complexity(self, c):
        self._complexity = c

    @property
    def Status(self):
        return self._status
    @Status.setter
    def Status(self, s):
        if s == 'Damaged' or s == 1:
            self._status = 1
        if s == 'Repaired' or s == 2:
            self._status = 2
        if s == 'Misrepaired' or s == 3:
            self._status = 3

    @property
    def RepairProteinStatus(self):
        return self._repprotstatus
    @RepairProteinStatus.setter
    def RepairProteinStatus(self, s):
        if s == 'Present' or s == 1 or s:
            self._repprotstatus = 1
        if s == 'NotPresent' or s == 0 or s is False:
            self._repprotstatus = 0

class DamageTrack:
    def __init__(self, trackid = -1):
        self.TrackID = trackid
        self.Status = 'Damaged'

    def GetTrackID(self):
        return self.TrackID

    @property
    def TrackID(self):
        return self._trackid

    @TrackID.setter
    def TrackID(self, id):
        self._trackid = id

    @property
    def ChromosomeID(self):
        return self._chrId

    @ChromosomeID.setter
    def ChromosomeID(self, c):
        self._chrId = c

    @property
    def BasePairID(self):
        return self._bpId

    @BasePairID.setter
    def BasePairID(self, b):
        self._bpId = b

    @property
    def StrandID(self):
        return self._strandid

    @StrandID.setter
    def StrandID(self, s):
        self._strandid = s

    @property
    def Time(self):
        return self._time
    @Time.setter
    def Time(self, t):
        self._time = t

    @property
    def Position(self):
        return self._position
    @Position.setter
    def Position(self, p):
        self._position = p

    @property
    def Complexity(self):
        return self._complexity
    @Complexity.setter
    def Complexity(self, c):
        self._complexity = c

    @property
    def Status(self):
        return self._status
    @Status.setter
    def Status(self, s):
        if s == 'Damaged' or s == 1:
            self._status = 1
        if s == 'Repaired' or s == 2:
            self._status = 2
        if s == 'Misrepaired' or s == 3:
            self._status = 3

class BeTrack(DamageTrack):
    def __init__(self, trackid = -1):
        super().__init__(trackid)
        self.Steps = []

    def AddNewStep(self, step):
        self.Steps.append(step)

    def GetLastStep(self):
        if len(self.Steps) > 0:
            return self.Steps[-1]
        else:
            print("No steps were found for this break end track.")
            return -1

    def GetStepAtIndex(self, i):
        if i >=0 and i < len(self.Steps):
            return self.Steps[i]
        else:
            print("Error at accessing at step ", str(i), " in this break end track.")
            return -1

    @property
    def OriginTime(self):
        if len(self.Steps) > 0:
            return self.Steps[0].Time
        else:
            print ("Error at getting origin time, this break end has no steps.")
            return -1

    @property
    def IsRepaired(self):
        if len(self.Steps) > 0:
            for s in self.Steps:
                if s.Status == 2:
                    return True
            return False
        else:
            return -1

    @property
    def IsMisrepaired(self):
        if len(self.Steps) > 0:
            for s in self.Steps:
                if s.Status == 3:
                    return True
            return False
        else:
            return -1

    @property
    def RepairTime(self):
        if len(self.Steps) > 0 and self.IsRepaired:
            for s in self.Steps:
                if s.Status == 2:
                    return s.Time
        else:
            return -1

    @property
    def MisrepairTime(self):
        if len(self.Steps) > 0 and self.IsMisrepaired:
            for s in self.Steps:
                if s.Status == 3:
                    return s.Time
        else:
            return -1
