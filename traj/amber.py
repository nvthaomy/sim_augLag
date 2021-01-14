#!/usr/bin/env python

### Amber trajectory format.
### coded by MSS

import numpy as np

import base
                
                    

class AmberWrite(base.TrajWriteClass):
    """Raw trajectory format writer."""

    def __init__(self, TrjFile, **kwargs):
        """Opens an AMBER coordinate trajectory format for writing.
This will not write a PrmtopFile."""
        base.TrajWriteClass.__init__(self, TrjFile, Mode = "wb", **kwargs)
        self.__WroteHead = False
        self.__Dim = None
        self.__NAtom = None

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the coordinate file."""
        if self.__WroteHead == False:
            self.fobj.write("ACE".ljust(80) + "\n")
            self.__Dim, self.__NAtom = Pos.shape
        else:
            if not (self.__Dim, self.__NAtom) == Pos.shape:
                raise ValueError("Found inconsistent dimensions in position array.")
        #make strings of the values
        if self.Units:
            Scale = 1. / self.Units.Angstrom
        else:
            Scale = 1.
        Vals = ["%8.3f" % x for x in Pos.flat * Scale]
        #we need to add trailing spaces so that there are ten fields per line
        r = len(Vals) % 10
        if r > 0: 
            Vals.extend([" "*10] * (10 - r))
        for i in xrange(0, len(Vals), 10):
            self.fobj.write("".join(Vals[i:i+10]) + "\n")
            

class Amber(base.TrajClass):
    """Amber trajectory format reader; can handle compressed or normal."""
    TrajAttrDflt = ["Seq", "AtomRes", "HasBox"]
    EneFieldsConvert = {"E_pot":"PEnergy", "EKinetic":"KEnergy", "time(ps)":"Time"}
    
    def __init__(self, TrjFile, PrmtopFile,
                 EneFile = None, LogWeightFile = None,
                 Sys = None, Mask = None):
        """Amber trajectory format reader; can handle compressed or normal."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.TrjFile = TrjFile
        self.PrmtopFile = PrmtopFile
        self.EneFile = EneFile
        self.LogWeightFile = LogWeightFile
        self.__fobj = None
        self.__Enefobj = None
        self.__ParseInit()
        self.SetSlice()
        self.Reset()
        self.TrajAttr.update(Amber.TrajAttrDflt)

    def __ParseInit(self):
        def GetPrmtopBlock(Prmtop, Flag, DoStrip = True):
            "Gets Prmtop data for a specified Flag."
            import re
            p = re.compile(r"\%FLAG "
                         + re.escape(Flag)
                         + r"\s*\%FORMAT\([\d\.]*[\D]*([\d\.]*)\)\s*$([^\%]*)",
                         re.MULTILINE)
            m = p.search(Prmtop)
            if m is None:
                raise ValueError("Could not parse block %s in prmtop file." % Flag)
            BlockLen = int(m.group(1))
            Dat = m.group(2).replace("\n","")
            Blocks = [Dat[i:i+BlockLen] for i in range(0,len(Dat)-1,BlockLen)]
            if DoStrip:
                Blocks = [x.strip() for x in Blocks]
            return Blocks
        #read prmtop
        Prmtop = base.FileOpen(self.PrmtopFile, "rU").read()
        #get box information
        Pointers = GetPrmtopBlock(Prmtop, "POINTERS")
        IfBox = int(Pointers[27])
        if IfBox == 0:
            self.HasBox = False
        elif IfBox == 1:
            self.HasBox = True
        else:
            raise TypeError("Cannot parse AMBER trajectory with IFBOX = %d" % IfBox)
        #extract names, etc
        self.AtomNames = GetPrmtopBlock(Prmtop, "ATOM_NAME")
        self.NAtom = len(self.AtomNames)
        self.Seq = GetPrmtopBlock(Prmtop, "RESIDUE_LABEL")
        ResPtr = GetPrmtopBlock(Prmtop, "RESIDUE_POINTER")
        ResPtr = [int(x) for x in ResPtr] + [self.NAtom + 1]
        self.AtomRes = []
        for i in range(0, len(ResPtr) - 1):
            self.AtomRes.extend([i] * (ResPtr[i+1] - ResPtr[i]))
        #get byte sizes for the file
        self.__Bytes = base.FileLen(self.TrjFile)
        f = base.FileOpen(self.TrjFile)
        LinesPerFrame = (3 * self.NAtom + 9) / 10
        self.__BytesHead = len(f.readline())
        self.__BytesFrame = 0
        for i in range(LinesPerFrame):
            self.__BytesFrame += len(f.readline())
        if self.HasBox:
            self.__BytesBox = len(f.readline())
            self.__BytesFrame += self.__BytesBox
        else:
            self.__BytesBox = 0
        f.close()
        if not ((self.__Bytes - self.__BytesHead) % self.__BytesFrame) == 0:
            print "Total bytes: %d" % self.__Bytes
            print "Head bytes : %d" % self.__BytesHead
            print "Frame bytes: %d" % self.__BytesFrame
            raise ValueError("Trajectory does not have consistent frame sizes or bytes per frame.")
        self.NFrameAbs = (self.__Bytes - self.__BytesHead) / self.__BytesFrame
        #check if an ene file is specified
        if not self.EneFile is None:
            LVal = -1
            self.__EneFields = []
            self.__EneBytesHead = 0
            self.__EneBytesFrame = 0
            self.__EneBytes = base.FileLen(self.EneFile)
            f = base.FileOpen(self.EneFile)
            LinesPerFrame = 0
            while True:
                s = f.readline()
                NewLVal = int(s.split()[0][1:])
                if NewLVal < LVal:
                    break
                else:
                    LVal = NewLVal
                self.__EneBytesHead += len(s)
                LinesPerFrame += 1
                self.__EneFields.extend(s.split()[1:])
            self.__EneBytesFrame += len(s)
            for i in range(LinesPerFrame - 1):
                self.__EneBytesFrame += len(f.readline())
            f.close()
            if not ((self.__EneBytes - self.__EneBytesHead) % self.__EneBytesFrame) == 0:
                raise ValueError("Ene file does not have consistent frame sizes or bytes per frame.")
            if not (self.__EneBytes - self.__EneBytesHead) / self.__EneBytesFrame == self.NFrameAbs:
                raise ValueError("Number of energy fields does not equal number of frames.")
        #check if a Log Weight file is specified
        if not self.LogWeightFile is None:
            f = base.FileOpen(self.LogWeightFile)
            s = f.read()
            f.close()
            try:
                self.LogWeight = [float(x) for x in s.strip().split()]
                self.LogWeight = np.array(self.LogWeight, float)
            except ValueError:
                raise ValueError("Could not process log weight file.")
            if not len(self.LogWeight) == self.NFrameAbs:
                raise ValueError("Number of log weights does not equal number of frames.")
            

    def Close(self):
        if not self.__fobj is None:
            self.__fobj.close()
            self.__fobj = None
        if not self.__Enefobj is None:
            self.__Enefobj.close()
            self.__Enefobj = None

    def ParseFrame(self, IndexAbs):
        if self.Units:
            Scale = self.Units.Angstrom
        else:
            Scale = 1.
        if self.__fobj is None:
            self.__fobj = base.FileOpen(self.TrjFile)
        fpos = self.__BytesHead + IndexAbs * self.__BytesFrame
        if not self.__fobj.tell() == fpos:
            self.__fobj.seek(fpos)
        s = self.__fobj.read(self.__BytesFrame)
        self.FrameData = {}
        if self.HasBox:
            sbox = s[-self.__BytesBox:]
            s = s[:-self.__BytesBox]
            #get the box lengths
            try:
                sbox = sbox.replace("\n","")
                BoxL = np.fromstring(sbox, dtype=float, sep=" ")
                self.FrameData["BoxL"] = BoxL * Scale
            except ValueError:
                raise ValueError("Cannot parse box coordinates in trj file.")             
        try:
            s = s.replace("\n","")
            s = " ".join([s[i:i+8] for i in xrange(0, len(s), 8)])
            Pos = np.fromstring(s, dtype=float, count=3*self.NAtom, sep=" ")
            Pos = Pos.reshape(self.NAtom,3) * Scale
        except ValueError:
            raise ValueError("Cannot parse frame coordinates in trj file.")
        if not self.EneFile is None:
            if self.__Enefobj is None:
                self.__Enefobj = base.FileOpen(self.EneFile)
            fpos = self.__EneBytesHead + IndexAbs * self.__EneBytesFrame
            if not self.__Enefobj.tell() == fpos:
                self.__Enefobj.seek(fpos)
            s = self.__Enefobj.read(self.__EneBytesFrame)
            FieldInd = 0
            try:
                for line in s.split("\n"):
                    for val in line.split()[1:]:
                        self.FrameData[self.__EneFields[FieldInd]] = float(val)
                        FieldInd += 1
            except ValueError:
                raise ValueError("Cannot parse frame data in ene file.")
            for (k,v) in Amber.EneFieldsConvert.items():
                if k in self.FrameData:
                    self.FrameData[v] = self.FrameData[k]
        if not self.LogWeightFile is None:
            self.FrameData["LogWeight"] = self.LogWeight[IndexAbs]
        return Pos

    
        