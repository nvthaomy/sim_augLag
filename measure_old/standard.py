#/usr/bin/env python


### Energy measuring objects in SIM suite.
### coded by MSS

import numpy as np

import base



class KEnergy(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring kinetic energy."""
        Source = """
>>> init
Val = KEnergy
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "KEnergy", StepFreq, CycleFreq, Active,
                                      Source = Source)


class PEnergy(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring potential energy."""
        Source =  """
>>> init
Val = PEnergy
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "PEnergy", StepFreq, CycleFreq, Active,
                                      Source = Source)

class TEnergy(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring total energy."""
        Source =  """
>>> init
Val = TEnergy
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "TEnergy", StepFreq, CycleFreq, Active,
                                      Source = Source)


class PETerm(base.MeasureAvgClass):
    def __init__(self, Sys, Name, Ind, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring potential energy from one forcefield term."""
        Source = """
>>> init
Val = Terms(%d)
[DOMEASURE]""" % Ind
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, Active,
                                      Source = Source)


class Vol(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring number density."""
        Source = """
>>> init
Val = product(BoxL)
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "Vol", StepFreq, CycleFreq, Active, Source = Source)
        

class NMol(base.MeasureAvgClass):
    def __init__(self, Sys, MolID, MolName, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring number of molecules of a given type."""
        Source = """
>>> init
Val = float(NActiveMID(%d))
[DOMEASURE]""" % MolID
        base.MeasureAvgClass.__init__(self, Sys, "N_"+MolName, StepFreq, CycleFreq, Active,
                                      Source = Source)
        
class Rho(base.MeasureAvgClass):
    def __init__(self, Sys, MolID, MolName, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring number density."""
        Source = """
>>> init
Val = float(NActiveMID(%d)) / product(BoxL)
[DOMEASURE]""" % MolID
        base.MeasureAvgClass.__init__(self, Sys, "r_"+MolName, StepFreq, CycleFreq, Active,
                                      Source = Source)
        
class MolFrac(base.MeasureAvgClass):
    def __init__(self, Sys, MolID, MolName, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring mole fraction."""
        Source = """
>>> init
Val = float(NActiveMID(%d)) / sum(NActiveMID)
[DOMEASURE]""" % MolID
        base.MeasureAvgClass.__init__(self, Sys, "x_"+MolName, StepFreq, CycleFreq, Active,
                                      Source = Source)        


class KTemp(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring temperature from kinetic energy."""
        Source = """
>>> inputmaps
bool ConservesMomentum = Sys.Flags.ConservesMomentum
>>> init
if (merge(NDOF-Dim, NDOF, ConservesMomentum) <= 0) then
    Val = 0.
else
    Val = 2.d0 * KEnergy / (kB * merge(NDOF-Dim, NDOF, ConservesMomentum))
endif
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "KTemp", StepFreq, CycleFreq, Active, 
                                      Source = Source)


class Pressure(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring pressure."""
        Source = """
>>> init
Val = (kB*NDOF*TempSet - Virial) / (Dim * product(BoxL))
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "Pressure", StepFreq, CycleFreq, Active,
                                      Source = Source)
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcVirial = True
            

class Virial(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring virial."""
        Source = """
>>> init
Val = Virial
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "Virial", StepFreq, CycleFreq, Active,
                                      Source = Source)
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcVirial = True        


class DUParam(base.MeasureAvgClass):
    
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring DUParam."""
        object.__setattr__(self, "LibVars", {})
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.Name = "DUParam"
        self.Sys = Sys

    def PreLoad(self):
        """Run before library making."""
        N1, N2 = self.Sys.ForceField.ArraySizes()
        Source = """
>>> init
Val = DUParam
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, self.Sys, "DUParam", self.StepFreq,
                                      self.CycleFreq, self.Active, NVal = N1,
                                      UseCovar = True, Report = False,
                                      Source = Source)
        base.MeasureAvgClass.PreLoad(self)
        
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcDUParam = True        

    def ResetAvgs(self):
        """Resets average counters to zero."""
        base.MeasureAvgClass.ResetAvgs(self)


class DDUParam(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring DDUParam."""
        object.__setattr__(self, "LibVars", {})
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.Name = "DDUParam"
        self.__MakeMatrix = None
        self.Sys = Sys

    def PreLoad(self):
        """Run before library making."""
        N1, N2 = self.Sys.ForceField.ArraySizes()
        Source = """
>>> init
Val = DDUParam
[DOMEASURE]"""        
        base.MeasureAvgClass.__init__(self, self.Sys, "DDUParam", self.StepFreq,
                                      self.CycleFreq, self.Active, NVal = N2,
                                      Report = False, Source = Source)
        base.MeasureAvgClass.PreLoad(self)

    def PostLoad(self):
        """Run after library making."""
        base.MeasureAvgClass.PostLoad(self)
        #copy the function over so we can convert array to matrix
        self.__MakeMatrix = self.Sys.ForceField.Param.DDU.MakeMatrix
        
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcDUParam = True           

    def __getattr__(self, name):
        if name == "MatrixAvg":
            if self.__MakeMatrix is None:
                raise AttributeError("Can only access Matrix after compilation.")
            MatSum = self.__MakeMatrix(self.ValSum)
            return MatSum / (self.Count + 1.e-300)
        else:
            return base.MeasureAvgClass.__getattr__(self, name)



class DWParam(base.MeasureAvgClass):
    
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring DWParam."""
        object.__setattr__(self, "LibVars", {})
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.Name = "DWParam"
        self.Sys = Sys

    def PreLoad(self):
        """Run before library making."""
        N1, N2 = self.Sys.ForceField.ArraySizes()
        Source = """
>>> init
Val = DWParam
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, self.Sys, "DWParam", self.StepFreq,
                                      self.CycleFreq, self.Active, NVal = N1,
                                      UseCovar = False, Report = False,
                                      Source = Source)
        base.MeasureAvgClass.PreLoad(self)
        
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcDWParam = True           

    def ResetAvgs(self):
        """Resets average counters to zero."""
        base.MeasureAvgClass.ResetAvgs(self)


class DUParamDWParam(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring DUParam*DWParam."""
        object.__setattr__(self, "LibVars", {})
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.Name = "DUParamDWParam"
        self.__MakeMatrix = None
        self.__N1 = None
        self.Sys = Sys

    def PreLoad(self):
        """Run before library making."""
        N1, N2 = self.Sys.ForceField.ArraySizes()
        self.__N1 = N1
        Source = """
>>> init
do i = 0, NDParam - 1
    do j = 0, NDParam - 1
        Val(i*NDParam + j) = DUParam(i) * DWParam(j)
    enddo
enddo
[DOMEASURE]
"""        
        base.MeasureAvgClass.__init__(self, self.Sys, "DUParamDWParam", self.StepFreq,
                                      self.CycleFreq, self.Active, NVal = N1*N1,
                                      Report = False, Source = Source)
        base.MeasureAvgClass.PreLoad(self)
        
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: 
            self.Sys.Flags.CalcDUParam = True   
            self.Sys.Flags.CalcDWParam = True 

    def __getattr__(self, name):
        if name == "MatrixAvg":
            MatSum = self.ValSum.reshape((self.__N1, self.__N1), order = 'C')
            return MatSum / (self.Count + 1.e-300)
        else:
            return base.MeasureAvgClass.__getattr__(self, name)
        
        
class FluctTerm(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring the fluctuation energy term."""
        Source = """
>>> init
Val = FluctTerm
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "FluctTerm", StepFreq, CycleFreq, Active,
                                      Source = Source)        
        
class FluctE0(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring the mean fluctuation energy."""
        Source = """
>>> init
Val = FluctE0
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "FluctE0", StepFreq, CycleFreq, Active,
                                      Source = Source) 
        
class FluctA0(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring the variance of the fluctuation energy."""
        Source = """
>>> init
Val = FluctA0
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, "FluctA0", StepFreq, CycleFreq, Active,
                                      Source = Source)           
            
            
class DDWParam(base.MeasureAvgClass):
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = False):
        """Class for measuring DDWParam."""
        object.__setattr__(self, "LibVars", {})
        self.StepFreq = StepFreq
        self.CycleFreq = CycleFreq
        self.Active = Active
        self.Name = "DDWParam"
        self.__MakeMatrix = None
        self.Sys = Sys

    def PreLoad(self):
        """Run before library making."""
        N1, N2 = self.Sys.ForceField.ArraySizes()
        Source = """
>>> init
Val = DDWParam
[DOMEASURE]"""        
        base.MeasureAvgClass.__init__(self, self.Sys, "DDWParam", self.StepFreq,
                                      self.CycleFreq, self.Active, NVal = N2,
                                      Report = False, Source = Source)
        base.MeasureAvgClass.PreLoad(self)

    def PostLoad(self):
        """Run after library making."""
        base.MeasureAvgClass.PostLoad(self)
        #copy the function over so we can convert array to matrix
        self.__MakeMatrix = self.Sys.ForceField.Param.DDW.MakeMatrix
        
    def Init(self):
        """Initializes measure, used before integration."""
        if self.Active: self.Sys.Flags.CalcDWParam = True      

    def __getattr__(self, name):
        if name == "MatrixAvg":
            if self.__MakeMatrix is None:
                raise AttributeError("Can only access Matrix after compilation.")
            MatSum = self.__MakeMatrix(self.ValSum)
            return MatSum / (self.Count + 1.e-300)
        else:
            return base.MeasureAvgClass.__getattr__(self, name)
            

class EinsteinSDC(base.MeasureClass):
    
    ModuleVars = ["PosOrigin", "Sums", "Counts", "MinStep", "MaxStep", "OriginDeltaStep"]
    ModuleConsts = ["NOrigins"]
    
    def __init__(self, Sys, StepFreq = 100, Active = False,
                 MinStep = 0, MaxStep = 10000,
                 NOrigins = 1, OriginDeltaStep = 1000):
        """Initializes a new self-diffusion coefficient measurement using the Einstein approach.
StepFreq: frequency of updates to SDC in steps
MinStep: minimum step number for computing SDC
MaxStep: maximum step number for computing SDC
NOrigins: number of time origins
OriginDeltaStep: steps between time origins"""
        Source = """
>>> defs
int SDCOr
int SDCInd
float SDCrsq
>>> final
do SDCOr = 0, NOrigins - 1
    if (StepNum == OriginDeltaStep * SDCOr + StepFreq) then
        SDCOrigin(SDCOr, :, :) = Pos
    endif
    SDCInd = int((StepNum - MinStep - SDCOr * OriginDeltaStep) / StepFreq - 1) 
    if (SDCInd < 0 .or. SDCInd >= SDCNBin) cycle
    do i = 0, NAtom - 1
        if (MolActive(MInd(i)) /= 1) cycle
        SIDi = SID(i)
        SDCrsq = sum((Pos(i,:) - PosOrigin(SDCOr, i, :))**2)
        Sums(SIDi, SDCInd) = Sums(SIDi, SDCInd) + SDCrsq
        Counts(SIDi, SDCInd) = Counts(SIDi, SDCInd) + 1.d0
    enddo       
enddo
"""
        NVal = Sys.World.NSID 
        ValName = ["sdc:%s" % s.Label for s in Sys.World.SiteTypes]
        base.MeasureClass.__init__(self, Sys, "SDC", StepFreq = StepFreq,
                                   Active = Active, NVal = NVal, ValName = ValName,
                                   Source = Source)
        self.Setup(StepFreq = StepFreq, MinStep = MinStep, MaxStep = MaxStep, 
                   NOrigins = NOrigins, OriginDeltaStep = OriginDeltaStep) 
        
    def Setup(self, StepFreq = 100, MinStep = 0, MaxStep = 10000, 
              NOrigins = 1, OriginDeltaStep = 1000):
        """Sets up self-diffusion coefficient parameters.
StepFreq: frequency of SDC updates
MinStep: minimum step number for computing SDC
MaxStep: maximum step number for computing SDC
NOrigins: number of time origins
OriginDeltaStep: steps between time origins""" 
        if self.Sys.Loaded:
            raise ValueError("Cannot change SDC setup after loading.")
        StepFreq = int(StepFreq)
        self.StepFreq = StepFreq
        self.MinStep = np.array([int(MinStep / StepFreq) * StepFreq], dtype=int)
        self.MaxStep = np.array([int(MaxStep / StepFreq) * StepFreq], dtype=int)
        self.NOrigins = NOrigins
        self.OriginDeltaStep = np.array([int(OriginDeltaStep / StepFreq) * StepFreq], dtype=int)
        
    def PreLoad(self):
        """Run before library making."""
        NSID = self.Sys.World.NSID
        NBin = int((self.MaxStep - self.MinStep) / self.StepFreq)
        NAtom, Dim = self.Sys.Pos.shape
        self.PosOrigin = np.zeros((self.NOrigins, NAtom, Dim), dtype=float, order='F')
        self.Sums = np.zeros((NSID, NBin), dtype=float, order='F')
        self.Counts = np.zeros((NSID, NBin), dtype=float, order='F')    
        base.MeasureAvgClass.PreLoad(self)
        
    def Init(self):
        if not self.Active: return
        self.Avg.fill(0.)
        self.Val.fill(0.)
        self.Sums.fill(0.)
        self.Counts.fill(0.)
        if hasattr(self.Sys.Int.Method, "TimeStep"):
            self.TimeStep = self.Sys.Int.Method.TimeStep
        else:
            self.TimeStep = 1
    
    def Finalize(self):
        if not self.Active: return
        NSID = self.Sys.World.NSID
        MSD = self.Sums / (self.Counts + 1.e-300)
        Times = np.arange(self.MinStep, self.MaxStep, self.StepFreq) * self.TimeStep
        for i in xrange(NSID):
            Mask = (self.Counts[i] > 0)
            if Mask.sum() < 2:
                SDC = 0
            else:
                Slope, Int = np.polyfit(Times[Mask], MSD[i,Mask], 1)
                SDC = Slope / (2. * self.Sys.Dim)
            self.Avg[i] = SDC
     
    def SDCString(self):
        """Returns a string of SDC results."""
        Labels = [s.Label for s in self.Sys.World.SiteTypes]
        NSID = self.Sys.World.NSID
        NTime = self.Sums.shape[1]
        Counts = self.Counts
        MSD = self.Sums / (self.Counts + 1.e-300)
        Times = np.arange(self.MinStep, self.MaxStep, self.StepFreq)
        TimeLabel = "step"
        if not self.TimeStep is None:
            Times *= self.TimeStep
            TimeLabel = "time"
        s = "%-20s %11s\n" % ("site", "sdc")
        for i in xrange(NSID):
            s += "%-20s %11.4e\n" % (Labels[i], self.Avg[i])
        s += "\n"
        s += "%-11s" % TimeLabel + " " + " ".join(["%-20s %-20s" % ("msd:"+l, "n:"+l) for l in Labels]) + "\n"
        for j in xrange(NTime):
            s += "%-11.4e" % Times[j] + " " + " ".join(["%-20.4e %-20.4e" % (MSD[i,j], Counts[i,j]) 
                                                        for i in xrange(NSID)]) + "\n"
        return s
        


def Add(Measures):
    """Adds standard measures to a measures class."""
    Sys = Measures.Sys
    MList = [KEnergy(Sys), PEnergy(Sys), TEnergy(Sys),
             KTemp(Sys), Vol(Sys), Pressure(Sys), Virial(Sys)]
    for MolID in range(Sys.World.NMID):
        M = NMol(Sys, MolID, Sys.World[MolID].Name)
        MList.append(M)    
    for MolID in range(Sys.World.NMID):
        M = Rho(Sys, MolID, Sys.World[MolID].Name, Active = False)
        MList.append(M)
    for MolID in range(Sys.World.NMID):
        M = MolFrac(Sys, MolID, Sys.World[MolID].Name, Active = False)
        MList.append(M)        
    for (Ind, P) in enumerate(Sys.ForceField):
        M = PETerm(Sys, P.Name, Ind, StepFreq = 0, CycleFreq = 1, Active = False)
        MList.append(M)
    MList.append(DUParam(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(DDUParam(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(DWParam(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(DDWParam(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(DUParamDWParam(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(FluctTerm(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(FluctE0(Sys, StepFreq = 0, CycleFreq = 1, Active = False))
    MList.append(FluctA0(Sys, StepFreq = 0, CycleFreq = 1, Active = False))    
    #MList.append(EinsteinSDC(Sys, Active = False))
    Measures.extend(MList)

