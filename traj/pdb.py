#!/usr/bin/env python

### Pdb file / trajectory format.
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



class PdbWrite(base.TrajWriteClass):
    """Pdb format for writing."""

    def __init__(self, FileName, **kwargs):
        """Opens a pdb trajectory format for writing."""
        base.TrajWriteClass.__init__(self, FileName, Mode = "wb", **kwargs)
        self.__NAtom = None
        self.__Time = 0.
        self.__ModelNum = 0

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the file."""
        if self.Units:
            Scale = 1. / self.Units.Angstrom
        else:
            Scale = 1.        
        if not "Time" in FrameData:
            FrameData["Time"] = self.__Time
        if self.__NAtom is None:
            sHead = ""
            if "NAtom" in HeadVars:
                sHead += "%d\n" % HeadVars["NAtom"]
                self.__NAtom = HeadVars["NAtom"]
            else:
                raise ValueError("Cannot find 'NAtom' in HeadVars.")
        else:
            if not len(Pos) == self.__NAtom:
                raise ValueError("Number of atoms changed.")
        if "AtomNames" in HeadVars:
            AtomNames = HeadVars["AtomNames"]
        else:
            raise ValueError("Cannot find 'AtomNames' in HeadVars.")
        AtomIsProtein = HeadVars.get("AtomIsProtein", [False]*self.__NAtom)
        AtomRes = HeadVars.get("AtomRes", ["SYS"]*self.__NAtom)
        AtomResNum = HeadVars.get("AtomResNum", range(self.__NAtom))
        AtomChain = HeadVars.get("AtomChain", [" "]*self.__NAtom)
        AtomOccupancy = HeadVars.get("AtomOccupancy", [0.]*self.__NAtom)
        AtomBFactor = HeadVars.get("AtomBFactor", [0.]*self.__NAtom)
        if "AtomElement" in HeadVars:
            AtomElement = HeadVars["AtomElement"]
        else:
            AtomElement = [x.strip("0123456789") for x in AtomNames]
        AtomCharge = HeadVars.get("AtomCharge", [0]*self.__NAtom)
        if "Angstrom" in HeadVars:
            Pos = Pos.copy() * Scale       
        #header part
        self.__ModelNum += 1
        s = "MODEL     %-4i\n" % self.__ModelNum
        if "REMARK" in FrameData:
            s += "REMARK %s\n" % FrameData["REMARK"]
        #add the coordinates to the pdb string
        for i in range(self.__NAtom):
            x, y, z = Pos[i]
            an = AtomNames[i].strip().ljust(4)[:4]
            ar = AtomRes[i]
            arn = AtomResNum[i]
            ac = AtomChain[i]
            ao = AtomOccupancy[i]
            ab = AtomBFactor[i]
            ae = AtomElement[i].rjust(2)[:2]
            ach = AtomCharge[i]
            if AtomIsProtein[i]:
                s += "ATOM  "
            else:
                s += "HETATM"
            s += "%5d %4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           %2s%2d\n" % (
                 i+1, an, ar, ac, arn+1, x, y, z, ao, ab, ae, ach)               
        #make the footer
        s += "TER\nENDMDL\n"
        self.fobj.write(s)


class Pdb(base.TrajClass):
    """Pdb file and trajectory format reader; can handle compressed or normal."""
    TrajAttrDflt = ["NAtom", "AtomTypes", "AtomNames", "AtomRes", "AtomResNum",
                    "AtomChain", "AtomOccupancy", "AtomBFactor", "AtomElement",
                    "AtomCharge"]
    
    def __init__(self, PdbFile, Sys = None, Mask = None):
        """Pdb file and trajectory format reader; can handle compressed or normal."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask)
        self.PdbFile = PdbFile
        self.__fobj = None
        self.__ParseInit()
        self.SetSlice()
        self.Reset()
        self.TrajAttr.update(Pdb.TrajAttrDflt)
        
    
    def __ParsePdb(self, FrameString):
        Pos = []
        AtomNames = []
        AtomRes = []
        AtomResNum = []
        AtomChain = []
        AtomOccupancy = []
        AtomBFactor = []
        AtomElement = []
        AtomCharge = []
        Seq = []
        #sort through lines
        CurResSeq = ""
        CurResNum = -1
        for line in FrameString.split("\n"):
            if not line[:6] in ["ATOM  ", "HETATM"]: continue
            AtomNames.append(line[12:16].strip())
            AtomRes.append(line[17:20].strip())
            AtomChain.append(line[21:22])
            ResSeq = line[22:26].strip()
            if not CurResSeq == ResSeq:
                CurResNum += 1
                CurResSeq = ResSeq
                Seq.append(AtomRes[-1])
            AtomResNum.append(CurResNum)
            x, y, z = line[30:38], line[38:46], line[46:54]
            try:
                x, y, z = float(x), float(y), float(z)
            except ValueError:
                raise ValueError("Could not convert x, y, z values in pdb file to floats.")
            Pos.append([x, y, z])
            AtomOccupancy.append(line[54:60].strip())
            BFactor = line[60:66].strip()
            if len(BFactor):
                try:
                    BFactor = float(BFactor)
                except ValueError:
                    raise ValueError("Could not convert B-factor in pdb file to float.")
                AtomBFactor.append(BFactor)
            else:
                AtomBFactor.append(0.)
            Element = line[76:78].strip()
            if len(Element):
                Element = Element.strip("0123456789")
            AtomElement.append(Element)
            Charge = line[78:80].strip()
            if len(Charge):
                if Charge[-1] in "+-":
                    Charge = Charge[::-1]
                try:
                    Charge = int(Charge)
                except ValueError:
                    raise ValueError("Could not convert charge in pdb file to int.")
                AtomCharge.append(Charge)
            else:
                AtomCharge.append(0)
        #convert to arrays
        Pos = np.array(Pos, dtype = float)
        AtomResNum = np.array(AtomResNum, dtype = int)
        AtomBFactor = np.array(AtomBFactor, dtype = float)
        AtomCharge = np.array(AtomCharge, dtype = int)        
        return (Pos, AtomNames, AtomRes, AtomResNum, AtomChain, AtomOccupancy,
                AtomBFactor, AtomElement, AtomCharge, Seq)
                
                
    def __ParseInit(self):
        from sim.utility import ProgressBar
        #find out how many frames
        f = base.FileOpen(self.PdbFile)    
        self.__FrameByteIndices = [0]
        #parse through all frames to get bytes in each frame
        #first line should be a label
        pb = ProgressBar("Getting trajectory frame byte indices.", BarLen = 0)
        CurPos = 0
        FoundModel = False
        line = f.readline()
        FirstFrame = ""
        #read through lines and count models and get indices
        while True:
            if line.startswith("MODEL"):
                if FoundModel == False:
                    FoundModel = True
                else:
                    #this is a new model, add the index
                    self.__FrameByteIndices.append(CurPos)
            if len(self.__FrameByteIndices) == 1:
                FirstFrame += line
            CurPos = f.tell()
            pb.Update(len(self.__FrameByteIndices))
            line = f.readline()
            if len(line) == 0:
                self.__FrameByteIndices.append(CurPos)
                break
        f.close()
        pb.Clear()        
        #parse the first frame to get info
        (self.Pos, self.AtomNames, self.AtomRes, self.AtomResNum, 
         self.AtomChain, self.AtomOccupancy, self.AtomBFactor, self.AtomElement, 
         self.AtomCharge, self.Seq) = self.__ParsePdb(FirstFrame)
        self.NAtom = len(self.Pos)        
        self.AtomTypes = self.AtomNames
        #get frame numbers
        self.NFrameAbs = len(self.__FrameByteIndices) - 1          

    def Close(self):
        if not self.__fobj is None:
            self.__fobj.close()
            self.__fobj = None

    def ParseFrame(self, IndexAbs):     
        if self.__fobj is None:
            self.__fobj = base.FileOpen(self.PdbFile)
        fpos = self.__FrameByteIndices[IndexAbs]
        BytesFrame = self.__FrameByteIndices[IndexAbs+1] - fpos
        if not self.__fobj.tell() == fpos:
            self.__fobj.seek(fpos)
        s = self.__fobj.read(BytesFrame)
        #parse the frame data
        (self.Pos, self.AtomNames, self.AtomRes, self.AtomResNum, 
         self.AtomChain, self.AtomOccupancy, self.AtomBFactor, self.AtomElement, 
         self.AtomCharge, self.Seq) = self.__ParsePdb(s)
        self.AtomTypes = self.AtomNames
        if self.Units:
            self.Pos *= self.Units.Angstrom        
        if not len(self.Pos) == self.NAtom:
            raise ValueError("Different number of atoms found in this frame.")
        return self.Pos

