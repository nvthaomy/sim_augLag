#!/usr/bin/env python

### LAMMPS trajectory format reader.
### coded by MSS

import numpy as np

import base

DEBUG = False



class LammpsWrite(base.TrajWriteClass):
    """Lammps format for writing."""

    def __init__(self, FileName, MinImage = False, ScalePos = False,
                 Precision = 6, **kwargs):
        """Opens a lammps trajectory format for writing."""
        base.TrajWriteClass.__init__(self, FileName, Mode = "wb", **kwargs)
        self.__Time = 0.
        self.MinImage = MinImage
        self.ScalePos = ScalePos
        self.Precision = Precision
        self.SysAttr.append("BoxL")

    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        """Writes a frame to the file."""
        if self.Units:
            Scale = 1. / self.Units.Angstrom
        else:
            Scale = 1.
        TimeStep = FrameData.get("TimeStep", self.__Time)
        if "AtomTypes" in HeadVars:
            AtomTypes = HeadVars["AtomTypes"]
        elif "AtomTypes" in FrameData:
            AtomTypes = FrameData["AtomTypes"]
        else:
            AtomTypes = np.zeros(len(Pos), dtype=int)
        if "AtomNames" in HeadVars:
            AtomNames = HeadVars["AtomNames"]
        elif "AtomNames" in FrameData:
            AtomNames = FrameData["AtomNames"]
        else:
            AtomNames = ["C"] * len(Pos)
        BoxL = FrameData["BoxL"]
        BoxLo = FrameData.get("BoxLo", -0.5 * BoxL)
        BoxHi = FrameData.get("BoxHi", 0.5 * BoxL)
        Pos = Pos.copy()    
        Mask = BoxL > 0
        if self.MinImage and self.ScalePos:
            Opt = "s"
            Pos1 = Pos / BoxL
            Pos1 = Pos1 - np.around(Pos1)
            Pos[:,Mask] = Pos1[:,Mask]
        elif self.MinImage:
            Opt = ""
            Pos1 = Pos - BoxL * np.around(Pos / BoxL)
            Pos[:,Mask] = Pos1[:,Mask]
        elif self.ScalePos:
            Opt = "su"
            Pos1 = Pos / BoxL
            Pos[:,Mask] = Pos1[:,Mask]
        else:
            Opt = "u"
        #scale positions
        Pos = Pos * Scale
        BoxL = BoxL * Scale
        BoxLo = BoxLo * Scale
        BoxHi = BoxHi * Scale
        self.__Time += 1.
        s = """ITEM: TIMESTEP
%d
ITEM: NUMBER OF ATOMS
%d
ITEM: BOX BOUNDS pp pp pp
%f %f
%f %f
%f %f
""" % (TimeStep, len(Pos), BoxLo[0], BoxHi[0],
       BoxLo[1], BoxHi[1], BoxLo[2], BoxHi[2])
        s += "ITEM: ATOMS id type x%s y%s z%s element\n" % (Opt, Opt, Opt)
        for (ind, (x, y, z)) in enumerate(Pos):
            s += "%d %d %.*f %.*f %.*f %s\n" % (ind+1, AtomTypes[ind], 
                                                self.Precision, x, 
                                                self.Precision, y,
                                                self.Precision, z,
                                                AtomNames[ind])
        self.fobj.write(s)


def ParseLog(LogFile, Token = None):
    """Parses log file results to extract thermo quantities."""
    Log = file(LogFile, "r").read()
    if not Token is None:
        ind = Log.rfind(Token)
        if ind < 0:
            raise ValueError("Could not find token %s in logfile %s" % (Token, LogFile))
        else:
            Log = Log[ind:]
    d = {}
    NFrameLog = 0
    for s in Log.split("---------------- Step")[1:]:
        s = s[s.find("\n")+1:]
        l = s.split()
        Found = False
        while len(l) >= 3:
            if l[1] == '=':
                key = l[0]
                val = float(l[2])
                if key in d:
                    d[key].append(val)
                else:
                    d[key] = [val]
                l = l[3:]
                if not Found:
                    NFrameLog += 1
                    Found = True
            elif Found:
                #done with this block
                break
            else:
                #shift one over
                l = l[1:]
    for key in d.keys():
        d[key] = np.array(d[key])
    return d, NFrameLog    

def ParseOMMLog(LogFile):
    with open(LogFile,'r') as file:
        header = file.readline()

    d = {}
    data = np.loadtxt(LogFile)

    for ind,column in enumerate(header.split("\t")):
        if "Potential" in column:
            PEindex = ind
            d.update({"PotEng":data[:,PEindex]})
        if "Temp" in column:
            Tindex = ind
            d.update({'Temp':data[:,Tindex]})
    
    NFrameLog = data.shape[0]

    return d, NFrameLog

class Lammps(base.TrajClass):
    """Lammps trajectory format reader; can handle compressed or normal."""
    TrajAttrDflt = ["NAtom", "AtomTypes", "AtomNames", "Init_BoxL"]
    EneFieldsConvert = {"PotEng":"PEnergy", "KinEng":"KEnergy", "TotEng":"TEnergy"}
    
    def __init__(self, TrjFile, LogFile = None, LogFileToken = None, OMMLog=False, 
                 Sys = None, Mask = None, **kwargs):
        """Lammps trajectory format reader; can handle compressed or normal."""
        base.TrajClass.__init__(self, Sys=Sys, Mask=Mask, **kwargs)
        self.TrjFile = TrjFile
        self.LogFile = LogFile
        self.LogFileToken = LogFileToken
        self.OMMLog = OMMLog
        self.__fobj = None
        self.__ParseInit()
        self.SetSlice()
        self.Reset()
        self.TrajAttr.update(Lammps.TrajAttrDflt)
        
    
    def __ParseBlocks(self, FrameString, IndexAbs):
        Blocks = FrameString.split("ITEM: ")
        ret = {}
        DoScale = False
        
        for Block in Blocks:
            if len(Block) == 0: continue
            try:
                Label, Contents = Block.split("\n", 1)
            except ValueError:
                print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                raise ValueError("Error parsing lammps frame.")
            Label = Label.strip()
            Contents = Contents.strip()
            
            if Label.startswith("TimeStep"):
                try:
                    ret["TimeStep"] = float(Contents.strip())
                except ValueError:
                    print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                    raise ValueError("Error parsing TimeStep in frame.")
                    
            elif Label.startswith("NUMBER OF ATOMS"):
                try:
                    ret["NAtom"] = int(Contents)
                except ValueError:
                    print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                    raise ValueError("Error parsing number of atoms in frame.")                
                    
            elif Label.startswith("BOX BOUNDS"):
                is_periodic = ['p' in i for i in Label.split('BOX BOUNDS')[1].split()]
                ret["BoxLo"] = []
                ret["BoxHi"] = []
                for (i, line) in enumerate(Contents.split("\n")):
                    line = line.strip()
                    #set non periodic box dimensions to zero
                    if not is_periodic[i]:
                        lo, hi = 0., 0.
                    else:
                        try:
                            lo, hi = line.split(" ")
                            lo, hi = float(lo), float(hi)
                        except ValueError:
                            print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                            raise ValueError("Error parsing box limits in frame.")
                    ret["BoxLo"].append(lo)
                    ret["BoxHi"].append(hi)
                ret["BoxLo"] = np.array(ret["BoxLo"], float) 
                ret["BoxHi"] = np.array(ret["BoxHi"], float)
                ret["BoxL"] = ret["BoxHi"] - ret["BoxLo"]
                
            elif Label.startswith("ATOMS "):
                #get the column headers
                Keys = Label.split()[1:]
                #number of columns
                NColumns = len(Keys)
                #checks on format
                if not Keys[0] == "id":
                    print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                    print Keys
                    raise ValueError("First key in ATOMS section is not id.")
                if not (Keys[1] == "type" or Keys[2] == "type"):
                    print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                    print Keys
                    raise ValueError("Second key in ATOMS section is not type.")
                # Check to see if the type is over one in the LAMMPS trajectory,
                #   if so, increment the index by 1. 
                if Keys[1] == "type":
                    increment = 0
                elif Keys[2] == "type":
                    increment = 1
                if Keys[(2+increment)] == 'xs' and Keys[(3+increment)] == 'ys' and Keys[(4+increment)] == 'zs':
                    DoScale == True
                elif Keys[(2+increment)] == 'x' and Keys[(3+increment)] == 'y' and Keys[(4+increment)] == 'z':
                    DoScale == False
                elif Keys[(2+increment)] == 'xsu' and Keys[(3+increment)] == 'ysu' and Keys[(4+increment)] == 'zsu':
                    DoScale == True
                elif Keys[(2+increment)] == 'xu' and Keys[(3+increment)] == 'yu' and Keys[(4+increment)] == 'zu':
                    DoScale = False
                else:
                    print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                    print Keys
                    raise ValueError("2nd-4th key in ATOMS section are not appropriate coordinates.")
                if 'element' in Keys:
                    ElementInd = Keys.index('element')
                else:
                    ElementInd = 1
                #get the data
                if NColumns == 5:
                    #if we are all numerics, we can do this fast
                    try:
                        a = np.fromstring(Contents, sep = ' ', dtype = float)
                    except ValueError:
                        print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                        raise ValueError("Error parsing atoms in frame block.")
                else:
                    #there is probably element data, need to do longer parse
                    try:
                        a = np.array(Contents.split(), dtype=str)
                    except ValueError:
                        print "ABSOLUTE FRAME INDEX %d" % IndexAbs
                        raise ValueError("Error parsing atoms in frame block.")
                #reshape data array
                a = a.reshape((-1,NColumns))
                AtomIDs = np.array(a[:,0], dtype=int)
                AtomIDs = AtomIDs - AtomIDs.min()
                NAtom = len(AtomIDs)
                if AtomIDs.max() >= NAtom:
                  raise ValueError("Found missing or nonconsecutive atom indices in frame block.")
                Posx = np.zeros(NAtom, dtype=float)
                Posy = np.zeros(NAtom, dtype=float)
                Posz = np.zeros(NAtom, dtype=float)
                Posx[AtomIDs] = np.array(a[:,2+increment], dtype=float)
                Posy[AtomIDs] = np.array(a[:,3+increment], dtype=float)
                Posz[AtomIDs] = np.array(a[:,4+increment], dtype=float)
                ret["Pos"] = np.concatenate((Posx[:,np.newaxis], Posy[:,np.newaxis], Posz[:,np.newaxis]), axis=1)
                AtomTypes = np.zeros(NAtom, dtype=int)
                AtomTypes[AtomIDs] = np.array(a[:,1], dtype=int)
                ret["AtomTypes"] = AtomTypes
                AtomNames = np.zeros(NAtom, dtype=a.dtype)
                AtomNames[AtomIDs] = a[:,ElementInd]
                ret["AtomNames"] = list(AtomNames.astype(str))
        
        if not ("Pos" in ret and "BoxLo" in ret and "BoxHi" in ret):
            if DEBUG: print FrameString
            raise ValueError("Could not find all of position and box size blocks in frame.")
                
        #scale by units
        if self.Units:
            Scale = self.Units.AngstromScale
        else:
            Scale = 1.
        ret["BoxLo"] *= Scale
        ret["BoxHi"] *= Scale
        ret["BoxL"] *= Scale
        
        #scale the atomic coordinates
        if DoScale:
            ret["Pos"] = ret["BoxL"][np.newaxis,:] * ret["Pos"] + ret["BoxLo"][np.newaxis,:]
        else:
            ret["Pos"] *= Scale
                
        if not ret["NAtom"] == len(ret["Pos"]):
            if DEBUG: print FrameString
            raise ValueError("Indicated number of atoms does not match coordinates length.")   
        
        return ret
                       
          
    def __ParseInit(self):
        from sim.utility import ProgressBar
        #get byte sizes for the file
        f = base.FileOpen(self.TrjFile)    
        self.__FrameByteIndices = [0]
        #parse through all frames to get bytes in each frame
        #first line should be a label
        pb = ProgressBar("Getting trajectory frame byte indices.", BarLen = 0)
        FirstLine = f.readline()
        FirstFrame = FirstLine
        FirstFrameActive = True
        CurPos = f.tell()
        NFrame = 1
        while True:
            NextLine = f.readline()
            if len(NextLine) == 0:
                break
            if NextLine == FirstLine:
                FirstFrameActive = False
                self.__FrameByteIndices.append(CurPos)
                NFrame += 1
            if FirstFrameActive:
                FirstFrame = FirstFrame + NextLine
            CurPos = f.tell()
            pb.Update(len(self.__FrameByteIndices))
        self.__FrameByteIndices.append(f.tell())
        f.close()
        pb.Clear()
        #parse the first frame to get info
        ret = self.__ParseBlocks(FirstFrame, 0)
        self.NAtom = ret["NAtom"]
        self.AtomTypes = ret["AtomTypes"]
        #extract names, etc
        self.AtomNames = ret["AtomNames"]
        #get frame numbers
        self.NFrameAbs = len(self.__FrameByteIndices) - 1
        # get initial box side dimension
        self.Init_BoxL = ret["BoxL"]
        #parse the log file
        if not self.LogFile is None:
            if self.OMMLog:
                self.ThermoDict,NFrameLog = ParseOMMLog(self.LogFile)
            else:
                self.ThermoDict, NFrameLog = ParseLog(self.LogFile, Token = self.LogFileToken)
            
            for key in self.ThermoDict.keys():
                if key in Lammps.EneFieldsConvert:
                    newkey = Lammps.EneFieldsConvert[key]
                    self.ThermoDict[newkey] = self.ThermoDict[key]
            if not NFrameLog == self.NFrameAbs:
                raise ValueError("LAMMPS log file has %d thermo outputs but trajectory has %d frames." % (NFrameLog, self.NFrameAbs))


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
        #parse the frame data
        self.FrameData = self.__ParseBlocks(s, IndexAbs)
        if not self.FrameData["NAtom"] == self.NAtom:
            raise ValueError("Different number of atoms found in this frame.")
        Pos = self.FrameData["Pos"]
        del self.FrameData["Pos"]
        del self.FrameData["AtomTypes"]
        del self.FrameData["AtomNames"]
        #add energies to frame data if any
        if hasattr(self, "ThermoDict"):
            for key, val in self.ThermoDict.items():
                self.FrameData[key] = val[IndexAbs]
        return Pos

    
        
