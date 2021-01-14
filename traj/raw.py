#!/usr/bin/env python

### Raw trajectory format.
### coded by MSS

import numpy as np
import cPickle

import base

#Trajectory format is as follows:
#  BytesHeadVars BytesField BytesCoord NField
#  FieldName1
#  FieldName2
#     ...
#  FieldNameM
#  CPICKLED DICTIONARY OF HEAD FIELDS
#  CPICKLED DICTIONARY OF FRAME DATA
#  RAW COORDS
#  CPICKLED DICTIONARY OF FRAME DATA
#  RAW COORDS
#    ...


class RawWrite(base.TrajWriteClass):
    """Raw trajectory format for writing."""

    def __init__(self, FileName, **kwargs):
        """Opens a raw trajectory format for writing."""
        base.TrajWriteClass.__init__(self, FileName, Mode = "wb", **kwargs)
        self.__WroteHead = False
        self.__Dim = None
        self.__NAtom = None
        self.__BytesHeadVars = None
        self.__BytesCoord = None
        self.__BytesField = None
        self.__Time = 0.

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the file."""
        if not "Time" in FrameData:
            FrameData["Time"] = self.__Time
        if self.__WroteHead == False:
            sHeadVars = cPickle.dumps(HeadVars, protocol = 2)
        sField = cPickle.dumps(FrameData, protocol = 2)
        sPos = Pos.tostring()
        if self.__WroteHead == False:
            self.__BytesHeadVars = len(sHeadVars)
            self.__BytesField = len(sField)
            self.__BytesCoord = len(sPos)
            self.fobj.write("%d %d %d\n" % (self.__BytesHeadVars,
                                            self.__BytesField,
                                            self.__BytesCoord))
            self.fobj.write(sHeadVars)
            self.__WroteHead = True
        else:
            if not len(sField) == self.__BytesField:
                raise ValueError("Field byte size changed.")
            if not len(sPos) == self.__BytesCoord:
                raise ValueError("Coord byte size changed.")
        self.fobj.write(sField + sPos)
        Time = FrameData.get("Time", self.__Time)
        self.__Time = Time + 1.        

class Raw(base.TrajClass):
    """Raw trajectory format reader."""

    def __init__(self, FileName, Sys = None, Mask = None):
        """Raw trajectory format reader."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.FileName = FileName
        self.__ParseInit()
        self.SetSlice()
        self.__fobj = None

    def __ParseInit(self):
        f = base.FileOpen(self.FileName)
        try:
            l = f.readline()
            s1, s2, s3 = l.split()[:3]
            self.__BytesHeadVars = int(s1)
            self.__BytesField = int(s2)
            self.__BytesCoord = int(s3)
        except ValueError:
            print "Byte allocation line:"
            print l
            raise ValueError("Cannot parse byte allocations.")
        #read head vars
        try:
            TrajAttrDict = cPickle.loads(f.read(self.__BytesHeadVars))
        except ValueError:
            raise ValueError("Cannot parse header fields.")
        for (k,v) in TrajAttrDict.items():
            setattr(self, k, v)
        self.TrajAttr.update(TrajAttrDict.keys())
        #read byte sizes
        self.__BytesFrame = self.__BytesField + self.__BytesCoord
        self.__BytesHead = f.tell()
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
        s = self.__fobj.read(self.__BytesField)
        try:
            self.FrameData = cPickle.loads(s)
        except ValueError:
            raise ValueError("Cannot parse trajectory frame fields.")
        s = self.__fobj.read(self.__BytesCoord)
        try:
            Pos = np.fromstring(s, dtype=float, count=self.NAtom*self.Dim,
                                sep = '').reshape((-1,self.Dim))
        except ValueError:
            raise ValueError("Cannot parse trajectory coordinates.")
        return Pos
        
    def Delete(self):
        """Deletes this trajectory file."""
        import os
        self.Close()
        if os.path.isfile(self.FileName):
            os.remove(self.FileName)

    
        