#/usr/bin/env python


### Distance and angle measuring objects in SIM suite.
### coded by MSS


import base
import numpy as np



class Distance(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, Name = "Distance",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring a distance."""
        if not Filter.N == 2:
            raise TypeError("Filter must be for atom pairs.")
        Source = """
>>> mainlooppair
VAL = DIJ
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)


class Angle(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, Name = "Angle",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring an angle."""
        if not Filter.N == 3 and Filter.Intra:
            raise TypeError("Filter must be for atom triples.")
        Source = """
>>> angleloop
VAL = THETA
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)
  
  
class Dihedral(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, Name = "Dihedral",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring a dihedral angle."""
        if not Filter.N == 4 and Filter.Intra:
            raise TypeError("Filter must be for atom quartets.")
        Source = """
>>> torsionloop
VAL = PHI
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)
                                      
class CoordNum(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, Cut, Name = "CoordNum",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring mean coordination number."""
        if not Filter.N == 2:
            raise TypeError("Filter must be for atom pairs.")
        Source = """
>>> defs
int CoordCount__
int CoordSum__
>>> init
CoordCount__ = 0
CoordSum__ = 0
>>> mainloopbeforepair
CoordCount__ = CoordCount__ + 1
>>> mainlooppair
if (DIJSQ < CutSq) CoordSum__ = CoordSum__ + 2
>>> final
VAL = float(CoordSum__) / (float(CoordCount__) + 1.e-10)
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)    
        self.CutSq = Cut * Cut
        self.ModuleVars.append("CutSq")

        
class CoordNumAuto(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, MinDistFrac, Name = "CoordNum",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring a distancemean coordination number."""
        if not Filter.N == 2:
            raise TypeError("Filter must be for atom pairs.")
        Source = """
>>> defs
float MinPairDistSq__(NAtom)
float CutSq__
int CoordCount__
int CoordSum__
>>> init
MinPairDistSq__ = 1.d300
CoordCount__ = 0
CoordSum__ = 0
>>> prelooppair
MinPairDistSq__(I) = min(MinPairDistSq__(I), DIJSQ)
MinPairDistSq__(J) = min(MinPairDistSq__(J), DIJSQ)
>>> beforemainloop
CutSq__ = (sum(MinPairDistSq__) / NAtom) * MinDistFrac
>>> mainloopbeforepair
CoordCount__ = CoordCount__ + 1
>>> mainlooppair
if (DIJSQ < CutSq__) CoordSum__ = CoordSum__ + 2
>>> final
VAL = float(CoordSum__) / (float(CoordCount__) + 1.e-10)
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)    
        self.MinDistFrac = MinDistFrac * np.ones(1, dtype=float)
        self.ModuleVars.append("MinDistFrac")     
        
        
class ClusterSize(base.MeasureAvgClass):
    def __init__(self, Sys, Filter, Cut, Name = "ClusterSize",
                 StepFreq = 0, CycleFreq = 1, Active = True):
        """Class for measuring mean cluster size"""
        if not Filter.N == 2:
            raise TypeError("Filter must be for atom pairs.")
        Source = """
>>> defs
int ClustID__(NAtom)
int NClust__
int NClustAtoms__
int ClustI
int ClustJ
>>> init
ClustID__ = 0
NClust__ = 0
NClustAtoms__ = 0
>>> preloopbeforepair
NClust__ = NClust__ + 1
NClustAtoms__ = NClustAtoms__ + 1
ClustID__(I) = NClust__
>>> mainlooppair
if (DIJSQ < CutSq) then
    ClustI = ClustID__(I)
    ClustJ = ClustID__(J)
    if (.not. ClustI == ClustJ) then
        do K = 1, NATOM
            if (ClustID__(K) == ClustJ) ClustID__(K) = ClustI
        enddo
        NClust__ = NClust__ - 1
    endif
endif
>>> final
VAL = float(NClustAtoms__) / (float(NClust__) + 1.e-10)
[DOMEASURE]"""
        base.MeasureAvgClass.__init__(self, Sys, Name, StepFreq, CycleFreq, 
                                      Active, Source = Source, Filter = Filter)    
        self.CutSq = Cut * Cut 
        self.ModuleVars.append("CutSq")      
        