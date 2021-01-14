#!/usr/bin/env python

### Simple trajectory format.
### coded by MSS

import numpy as np

import base


#Trajectory format is as follows:
#  NAtom
#  AtomName1
#  AtomName2
#    ...
#  AtomNameN
#  Time [PEnergy [LogWeight]]
#  x1 y1 z1
#  x2 y2 z2
#    ...
#  xN yN zN
#  Time [PEnergy [LogWeight]]
#  x1 y1 z1
#    ...



class SimpleWrite(base.TrajWriteClass):
    """Simple trajectory format for writing."""

    def __init__(self, FileName, **kwargs):
        """Opens a simple trajectory format for writing."""
        base.TrajWriteClass.__init__(self, FileName, Mode = "wb", **kwargs)
        self.__WroteHead = False
        self.__Time = 0.

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the file."""
        if not "Time" in FrameData:
            FrameData["Time"] = self.__Time
        if self.__WroteHead == False:
            sHead = ""
            if "NAtom" in HeadVars:
                sHead += "%d\n" % HeadVars["NAtom"]
                self.__NAtom = HeadVars["NAtom"]
            else:
                raise ValueError("Cannot find 'NAtom' in HeadVars.")
            if "AtomNames" in HeadVars:
                for i in range(self.__NAtom):
                    sHead += "%s\n" % HeadVars["AtomNames"][i]
            else:
                raise ValueError("Cannot find 'NAtom' in HeadVars.")
            self.fobj.write(sHead)
            self.__WroteHead = True
        else:
            if not len(Pos) == self.__NAtom:
                raise ValueError("Number of atoms changed.")
        Time = FrameData.get("Time", self.__Time)
        self.__Time = Time + 1.
        PEnergy = FrameData.get("PEnergy", 0.)
        LogWeight = FrameData.get("LogWeight", 0.)
        sFrame = "%14.7e %14.7e %14.7e\n" % (Time, PEnergy, LogWeight)
        for i in range(self.__NAtom):
            sFrame += " ".join(["%12.5e" % x for x in Pos[i]]) + "\n"
        self.fobj.write(sFrame)


class Simple(base.TrajClass):
    """Simple trajectory format reader; can handle gzipped or normal."""

    def __init__(self, FileName, Sys = None, Mask = None):
        """Simple trajectory format reader; can handle gzipped or normal."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.FileName = FileName
        self.__ParseInit()
        self.SetSlice()
        self.__fobj = None

    def __ParseInit(self):
        f = base.FileOpen(self.FileName)
        try:
            self.NAtom = int(f.readline())
        except ValueError:
            raise ValueError("Cannot parse trajectory file.")
        self.AtomNames = []
        for i in range(self.NAtom):
            self.AtomNames.append(f.readline().strip())
        self.__BytesHead = f.tell()
        for i in range(self.NAtom + 1):
            f.readline()
        self.__BytesFrame = f.tell() - self.__BytesHead
        f.close()
        self.__Bytes = base.FileLen(self.FileName)
        if not ((self.__Bytes - self.__BytesHead) % self.__BytesFrame) == 0:
            raise ValueError("Trajectory does not have consistent frame sizes or bytes per frame.")
        self.NFrameAbs = (self.__Bytes - self.__BytesHead) / self.__BytesFrame   

    def Close(self):
        if not self.__fobj is None:
            self.__fobj.close()
            self.__fobj = None

    def ParseFrame(self, IndexAbs):
        if self.__fobj is None:
            self.__fobj = base.FileOpen(self.FileName)
        fpos = self.__BytesHead + IndexAbs * self.__BytesFrame
        if not self.__fobj.tell() == fpos:
            self.__fobj.seek(fpos)
        line = self.__fobj.readline()
        vals = line.split()
        self.FrameData = {}
        try:
            self.FrameData["Time"] = float(vals[0])
            if len(vals) > 1:
                self.FrameData["PEnergy"] = float(vals[1])
            if len(vals) > 2:
                self.FrameData["LogWeight"] = float(vals[2])
        except ValueError:
            raise ValueError("Cannot parse trajectory frame data.")
        Data = self.__fobj.read(self.__BytesFrame - len(line))
        try:
            Pos = np.fromstring(Data, float, sep = ' ').reshape((-1,3))
        except ValueError:
            raise ValueError("Cannot parse trajectory coordinates.")        
        return Pos

    
        