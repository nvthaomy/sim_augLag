#/usr/bin/env python


### Trajectory routines objects in SIM suite.
### coded by MSS

import numpy as np
import copy
import os
import gzip, bz2

import sim.utility as utility


class TrajWriteClass(object):
    SysAttrDflt = ["PEnergy", "KEnergy", "TEnergy", "BoxL"]
    
    def __init__(self, FileName, Mode = "wb", Units = None):
        self.fobj = FileOpen(FileName, Mode)
        self.HeadVars = None
        self.SysAttr = copy.copy(TrajWriteClass.SysAttrDflt)
        self.Units = Units
        
    def WriteFrame(self, Pos, HeadVars = {}, FrameData = {}):
        pass
    
    def Close(self):
        self.fobj.close()

    def WriteFrameSys(self, Sys):
        self.Units = Sys.Units
        if self.HeadVars is None:
            self.HeadVars = {"NAtom"     : len(Sys.Atom),
                             "Dim"       : Sys.Dim,
                             "AtomNames" : [a.Name for a in Sys.Atom],
                             "TempSet"   : Sys.TempSet}
        FrameData = dict([(Attr, getattr(Sys, Attr))
                          for Attr in self.SysAttr])
        self.WriteFrame(Sys.Pos.copy(), HeadVars = self.HeadVars,
                        FrameData = FrameData)

    def AddAction(self, IntegrateObj, StepFreq = None, CycleFreq = None):
        """Adds an action to an integrator to output trajectory."""
        if hasattr(self, "Action"):
            raise ValueError("Action already added to an integrator.")
        self.Action = IntegrateObj.AddAction(Fn = self.WriteFrameSys,
                                             StepFreq = StepFreq,
                                             CycleFreq = CycleFreq,
                                             Name = "TrajWrite")
        return self.Action

    def DelAction(self, IntegrateObj):
        """Removes an action from an integrator."""
        IntegrateObj.DelAction(self.Action)
        del self.Action        


def Convert(InTraj, OutTrajClass, FileName, 
            Verbose = False, FrameDataFields = None, NewFrameData = None,
            FrameIndices = None, BlurSigma = None, BlurNSample = 1, *args, **kwargs):
    """Converts a given trajectory to another trajectory format."""
    if not isinstance(InTraj, TrajClass):
        raise TypeError("InTraj is not a TraClass instance.")
    if not issubclass(OutTrajClass, TrajWriteClass):
        raise TypeError("OutTrajClass is not a TrajWriteClass object.")
    if not NewFrameData is None:
        for (k,v) in NewFrameData.iteritems():
            if not len(v) == len(InTraj):
                raise ValueError("%s in NewFrameData is not the same length as InTraj." % k)
    OutTraj = OutTrajClass(FileName, *args, **kwargs)
    if Verbose:
        p = utility.ProgressBar("Converting trajectory", len(InTraj))
    HeadVars = InTraj.GetHeadVars()
    for (i, Pos) in enumerate(InTraj):
        FrameData = InTraj.FrameData
        if not FrameDataFields is None:
            FrameData = dict([(k, FrameData.get(k, 0.))
                              for k in FrameDataFields])
        if not NewFrameData is None:
            for (k,v) in NewFrameData.iteritems():
                FrameData[k] = NewFrameData[k][i]
        if BlurSigma is None:                
            OutTraj.WriteFrame(Pos, HeadVars, FrameData)
        else:
            for j in range(BlurNSample):
                BlurPos = Pos + np.random.normal(scale=BlurSigma, size=Pos.shape)
                OutTraj.WriteFrame(BlurPos, HeadVars, FrameData)
        if Verbose:
            p.Update(i)
    if Verbose:
        p.Clear()
    OutTraj.Close()
    return OutTraj
    
    
def ConvertList(InTrajList, OutTrajClass, FileName, 
                Verbose = False, FrameDataFields = None, 
                FrameIndices = None, *args, **kwargs):
    """Converts a given trajectory to another trajectory format.
InTraj is a list of trajectory objects.  FrameIndices is a list of tuples (TrajNum, FrameNum)."""
    if not all([isinstance(t, TrajClass) for t in InTrajList]):
        raise TypeError("InTraj is not a TraClass instance.")
    if not issubclass(OutTrajClass, TrajWriteClass):
        raise TypeError("OutTrajClass is not a TrajWriteClass object.")
    if FrameIndices is None:
        FrameIndices = []
        for (i, InTraj) in enumerate(InTrajList):
            FrameIndices.extend([(i,j) for j in range(len(InTraj))])
    OutTraj = OutTrajClass(FileName, *args, **kwargs)
    if Verbose:
        p = utility.ProgressBar("Converting trajectory list", len(FrameIndices))
    HeadVars = InTrajList[0].GetHeadVars()
    for (k, (i,j)) in enumerate(FrameIndices):
        InTraj = InTrajList[i]
        Pos = InTraj[j]
        FrameData = InTraj.FrameData
        if not FrameDataFields is None:
            FrameData = dict([(k, FrameData.get(k, 0.))
                              for k in FrameDataFields])
        OutTraj.WriteFrame(Pos, HeadVars, FrameData)
        if Verbose:
            p.Update(k)
    if Verbose:
        p.Clear()
    OutTraj.Close()
    return OutTraj    
    
    


class TrajClass(object):
    TrajAttrDflt = ["NAtom", "Dim", "AtomNames"]

    def __init__(self, Sys = None, Mask = None, Dim = 3, Units = None):
        """Base class for trajectory files."""
        #In derived classes, the following must be defined after calling:
        #  NAtom, AtomNames, NFrameAbs
        #Also, additional variables can be marked for use in copying
        #  and in headers by adding their names to the TrajAttr set
        #assign input variables if specified
        self.Sys = Sys
        self.Mask = Mask
        self.Dim = Dim
        self.Units = Units
        self.FrameData = {}
        self.Index = -1
        self.IndexAbs = -1
        self.Pos = None
        #see if a list or set of attributes for inclusion in the
        #header have been created; add default values 
        self.TrajAttr = set(TrajClass.TrajAttrDflt)
        #calls to SetSlice and Reset should be made after initialization
        # in derived classes

    def SetSlice(self, Start = None, Stop = None, Stride = None):
        if Start is None:
            self.sStart = 0
        else:
            self.sStart = Start
        if Stop is None:
            self.sStop = self.NFrameAbs
        else:
            self.sStop = Stop
        if Stride is None:
            self.sStride = 1
        else:
            self.sStride = Stride
        self.NFrame = max((self.sStop - self.sStart) / self.sStride, 0)
        if (self.sStop - self.sStart) % self.sStride > 0:
            self.NFrame += 1
            
    def GetHeadVars(self):
        """Returns a dictionary of variable names and values for header."""
        return dict([(f, getattr(self,f)) for f in sorted(self.TrajAttr)])

    def Copy(self):
        self.Reset()
        ret = copy.copy(self)
        for Attr in self.TrajAttr:
            setattr(ret, Attr, copy.deepcopy(getattr(self, Attr)))
        ret.Mask = copy.deepcopy(self.Mask)
        return ret

    def copy(self):
        """Alias for Copy()"""
        return self.Copy()

    def Close(self):
        #code here to safely close
        pass

    def close(self):
        """Alias for Close()"""
        return self.Close()

    def Reset(self):
        """Closes files and resets trajectory to beginning."""
        self.Close()
        self.FrameData = {}
        self.Index = -1
        self.IndexAbs = -1
        self.Pos = None

    def ParseFrame(self, IndexAbs):
        #code here to return Pos and set FrameData
        pass
       
    def Get(self, Index, Mask = None,
            UseMask = True, UpdateSys = True):
        """Returns a frame configuration at specified (relative) index."""
        if Index < 0: Index += self.NFrame
        if Index < 0 or Index >= self.NFrame:
            raise IndexError("Frame index out of range.")
        self.Index = Index
        self.IndexAbs = self.sStart + Index * self.sStride
        self.FrameData = {}
        Pos = self.ParseFrame(self.IndexAbs)
        if Mask is None: Mask = self.Mask
        if UseMask and not Mask is None:
            Pos = Pos[Mask]
        self.Pos = Pos
        if UpdateSys and not self.Sys is None:
            if not UseMask or Mask is None:
                self.Sys.Pos[:,:] = Pos
            else:
                self.Sys.Pos[Mask,:] = Pos
        return Pos

    def __getitem__(self, Index):
        """Returns a position array if Index is an integer;
a FrameData element if a string;
otherwise returns a new trajectory object with given slice."""
        if type(Index) is slice:
            Start, Stop, Stride = Index.start, Index.stop, Index.step
            if Start is None: Start = 0
            if Stop is None: Stop = self.NFrame
            if Stride is None: Stride = 1
            if Start < 0: Start += self.NFrame
            if Stop < 0: Stop += self.NFrame
            Start = max(0, Start)
            Stop = min(self.NFrame, Stop)
            if Stop < Start:
                Stop = Start
            #make a copy to return
            ret = self.Copy()
            #get absolute values
            ret.SetSlice(Start = self.sStart + Start * self.sStride,
                         Stop = self.sStart + Stop * self.sStride,
                         Stride = self.sStride * Stride)
            return ret
        elif type(Index) is str:
            return self.FrameData[Index]
        else:
            return self.Get(Index)
        
    def _setattr__(self, n, v):
        if n == "Sys":
            self.Units = v.Units
        object.__setattr__(n, v)

    def __len__(self):
        return self.NFrame

    def __iter__(self):
        self.Reset()
        return self

    def next(self):
        self.Index += 1
        if self.Index < self.NFrame:
          return self.Get(self.Index)
        else:
          self.Reset()
          raise StopIteration

    def GetIndices(self):
        "Returns the relative configuration indices."
        return range(self.NFrame)

    def GetIndicesAbs(self):
        """Returns the absolute configuration indices."""
        return range(self.sStart, self.sStop, self.sStride)
        


class Multi(TrajClass):
    
    def __init__(self, TrajList, Sys = None, Mask = None, NewFrameData = None):
        """Returns a multi-trajectory class."""
        TrajClass.__init__(self, Sys = Sys, Mask = Mask)
        self.TrajList = list(TrajList)
        if not len(self.TrajList):
            raise IndexError("TrajList is length zero.")
        if not all([isinstance(t, TrajClass) for t in self.TrajList]):
            raise TypeError("TrajList")
        #get any other settings from first object
        t = self.TrajList[0]
        for Attr in t.TrajAttr:
            setattr(self, Attr, copy.deepcopy(getattr(t, Attr)))
        #parse through trajectory list
        for Traj in self.TrajList:
            if self.AtomNames is None:
                self.AtomNames = copy.copy(Traj.AtomNames)
                self.NAtom = Traj.NAtom
            else:
                if not self.NAtom == Traj.NAtom:
                    raise TypeError("Not all trajectories have the same number of atoms.")
            if self.Dim is None:
                self.Dim = Traj.Dim
            else:
                if not self.Dim == Traj.Dim:
                    raise TypeError("Not all trajectories have the same dimensions.")                
        #check weights
        self.NewFrameData = None
        if not NewFrameData is None:
            for (k,v) in NewFrameData.iteritems():
                if not len(v) == len(self.TrajList):
                    raise TypeError("%s in NewFrameData is not the same length as TrajList." % k)
                elif any([len(t) != len(l) for (t,l) in zip(self.TrajList, NewFrameData[k])]):
                    raise TypeError("Not all trajectory lengths match %s length." % k)
            self.NewFrameData = NewFrameData
        self.__Traj = None
        self.__TrajNFrame = [t.NFrame for t in self.TrajList]
        self.__StartInd = [sum(self.__TrajNFrame[:i])
                           for i in range(len(self.TrajList) + 1)]
        self.NFrameAbs = sum(self.__TrajNFrame)
        self.SetSlice()
        self.FrameData = {}
        self.Index = -1
        self.IndexAbs = -1
        self.Pos = None    

    def Copy(self):
        self.Reset()
        ret = copy.copy(self)
        for Attr in self.TrajAttr:
            setattr(ret, Attr, copy.deepcopy(getattr(self, Attr)))
        ret.TrajList = copy.copy(self.TrajList)
        ret.Mask = copy.deepcopy(self.Mask)
        ret.Sys = None
        return ret

    def Close(self):
        if not self.__Traj is None:
            self.__Traj.Close()
            self.__Traj = None

    def Get(self, Index, Mask = None, UseMask = True, UpdateSys = True):
        if Index < 0: Index += self.NFrame
        if Index < 0 or Index >= self.NFrame:
            raise IndexError("Frame index out of range.")
        self.Index = Index
        self.IndexAbs = self.sStart + Index * self.sStride
        for (tNum, StartInd) in enumerate(self.__StartInd):
            if self.IndexAbs < self.__StartInd[tNum+1]:
                break
        Traj = self.TrajList[tNum]
        if self.__Traj is None:
            self.__Traj = Traj
        elif not Traj is self.__Traj:
            self.__Traj.Close()
            self.__Traj = Traj
        j = self.IndexAbs - StartInd
        Pos = self.__Traj.Get(j, UseMask = False, UpdateSys = False)
        self.FrameData = self.__Traj.FrameData.copy()
        if not self.NewFrameData is None:
            for (k,v) in self.NewFrameData.iteritems():
                self.FrameData[k] = self.NewFrameData[k][tNum][j]
        if Mask is None: Mask = self.Mask
        if UseMask and not Mask is None:
            Pos = Pos[Mask]
        self.Pos = Pos
        if UpdateSys and not self.Sys is None:
            if not UseMask or Mask is None:
                self.Sys.Pos[:,:] = Pos
            else:
                self.Sys.Pos[Mask,:] = Pos
        return Pos
            

def FileOpen(FileName, mode = "rb"):
    """Returns a file object with the correct class."""
    if FileName.endswith(".gz"):
        return gzip.GzipFile(FileName, mode)
    elif FileName.endswith(".bz2"):
        return bz2.BZ2File(FileName, mode)
    else:
        return file(FileName, mode)

def FileLen(FileName):
    """Returns the length of a file."""
    if FileName.endswith(".gz"):
        #for gzip, need to fake it by scanning thru file
        fobj = gzip.GzipFile(FileName, "rb")
        fobj.seek(0)
        n = 0
        inc = 1024
        while n == fobj.tell():
            n += inc
            fobj.seek(n)
        n = max(fobj.tell(), 0)
        fobj.close()
        return n
    elif FileName.endswith(".bz2"):
        fobj = bz2.BZ2File(FileName, "rb")
        fobj.seek(0, os.SEEK_END)
        n = fobj.tell()
        fobj.close()
        return n
    else:
        return os.path.getsize(FileName)

def FileCountLines(FileName):
    """Returns the number of lines in a file."""
    fobj = FileOpen(FileName)
    n = 0
    for line in fobj:
        n += 1
    fobj.close()
    return n
        
     
        

        
    
            
    