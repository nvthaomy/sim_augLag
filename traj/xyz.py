#!/usr/bin/env python

### xyz trajectory format reader.
### coded by MSS

import numpy as np

import base

DEBUG = False



class XYZWrite(base.TrajWriteClass):
    """XYZ format for writing."""

    def __init__(self, FileName, MinImage = False, ScalePos = False,
                 Precision = 6, **kwargs):
        """Opens a xyz trajectory format for writing."""
        base.TrajWriteClass.__init__(self, FileName, Mode = "wb", **kwargs)
        self.__FrameNum = 0
        self.MinImage = MinImage
        self.ScalePos = ScalePos
        self.Precision = Precision

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the file."""
        self.__FrameNum += 1
        TimeStep = FrameData.get("TimeStep", self.__FrameNum - 1)
        if "AtomNames" in HeadVars:
            AtomNames = HeadVars["AtomNames"]
        elif "AtomTypes" in FrameData:
            AtomNames = FrameData["AtomNames"]
        else:
            AtomNames = ["C"] * len(Pos)
        Pos = Pos.copy()    
        BoxL = FrameData["BoxL"]
        BoxLo = FrameData.get("BoxLo", -BoxL)
        BoxHi = FrameData.get("BoxHi", +BoxL)
        Mask = BoxL > 0
        if self.MinImage and self.ScalePos:
            Pos1 = (Pos - BoxLo) / BoxL
            Pos1 = Pos1 - np.around(Pos1)
            Pos[:,Mask] = Pos1[:,Mask]
        elif self.MinImage:
            Pos1 = Pos - BoxL * np.around(Pos / BoxL)
            Pos[:,Mask] = Pos1[:,Mask]
        elif self.ScalePos:
            Pos1 = Pos / BoxL
            Pos[:,Mask] = Pos1[:,Mask]
        s = """%d
frame: %d time:%14.7e
""" % (len(Pos), self.__FrameNum, TimeStep)
        for (ind, (x, y, z)) in enumerate(Pos):
            s += "%s %.*f %.*f %.*f\n" % (AtomNames[ind], 
                                          self.Precision, x, 
                                          self.Precision, y,
                                          self.Precision, z)
        self.fobj.write(s)



class XYZ(base.TrajClass):
    """XYZ trajectory format reader; can handle compressed or normal."""
    
    def __init__(self, TrjFile, Sys = None, Mask = None):
        """XYZ trajectory format reader; can handle compressed or normal."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.TrjFile = TrjFile
        self.__fobj = None
        self.__ParseInit()
        self.SetSlice()
        self.Reset()
          
    def __ParseInit(self):
        from sim.utility import ProgressBar
        #get byte sizes for the file
        f = base.FileOpen(self.TrjFile)    
        self.__FrameByteIndices = [0]
        self.AtomNames = []
        #parse through all frames to get bytes in each frame
        #first line should be a label
        pb = ProgressBar("Getting trajectory frame byte indices.", BarLen = 0)
        #parse the first frame
        self.NAtom = int(f.readline())
        f.readline()  #read the frame time line
        for i in xrange(self.NAtom):
            self.AtomNames.append(f.readline().split()[0])
        self.__FrameByteIndices.append(f.tell())
        NFrame = 1
        #now read other frames
        while True:
            FirstLine = f.readline()
            if len(FirstLine) == 0:
                break
            NFrame += 1
            if not int(FirstLine) == self.NAtom:
                raise ValueError("Different number of atoms found in frame %d" % NFrame)
            f.readline()  #read the frame time line
            for i in xrange(self.NAtom):
                f.readline()  #read the atom xyz line
            self.__FrameByteIndices.append(f.tell())
            pb.Update(len(self.__FrameByteIndices))
        f.close()
        pb.Clear()
        #get frame numbers
        self.NFrameAbs = len(self.__FrameByteIndices) - 1 
        
    def Close(self):
        if not self.__fobj is None:
            self.__fobj.close()
            self.__fobj = None

    def ParseFrame(self, IndexAbs):
        if self.__fobj is None:
            self.__fobj = base.FileOpen(self.TrjFile)
        fpos = self.__FrameByteIndices[IndexAbs]
        BytesFrame = self.__FrameByteIndices[IndexAbs+1] - fpos
        if not self.__fobj.tell() == fpos:
            self.__fobj.seek(fpos)
        s = self.__fobj.read(BytesFrame)
        s = s.rstrip()
        #parse the frame data
        Pos = []
        for line in s.split("\n")[2:]:
            Pos.append(line.split()[1:])
        Pos = np.array(Pos, dtype=float)
        self.FrameData = {}
        self.FrameData["Pos"] = Pos
        return Pos

    
        