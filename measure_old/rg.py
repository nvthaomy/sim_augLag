#/usr/bin/env python

### Radius of gyration measuring object in SIM suite.
### coded by SPC

import base

class Rg(base.MeasureAvgClass):
    
    ModuleVars = base.MeasureAvgClass.ModuleVars + ["MolInd"]

    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True, MolInd = 0):
        """Class for computing radius of gyration."""   
        Source = """
>>> defs
float PosSum(Dim)
float PosSumSq(Dim)
float RgN
float LastPos(Dim)
float ThisPos(Dim)
float dPos(Dim)
int ind

>>> init
!check molInd
if (MolInd < 0 .or. MolInd >= NMol) then
    print *, "Invalid molecule index specified in RG measure"
    stop
endif
!initialize averages
PosSum = 0.d0
PosSumSq = 0.d0
LastPos = 0.d0
DO ind = MolRange(MolInd), MolRange(MolInd+1)-1
    ThisPos = Pos(ind,:)
    !minimum image, if periodic
    if (DoMinImage) then
        !compute distance from last particle
        dPos = ThisPos - LastPos
        !minimum-image
        dPos = dPos - BoxL * dnint(dPos * iBoxL)
        !compute minimum-imaged position
        ThisPos = LastPos + dPos
    endif
    !add to running sums
    PosSum = PosSum + ThisPos
    PosSumSq = PosSumSq + ThisPos*ThisPos
    !update last position
    LastPos = ThisPos
END DO
!compute Rg
RgN = MolRange(MolInd+1) - MolRange(MolInd)
Val = sqrt( sum(PosSumSq / RgN - (PosSum*PosSum) / RgN**2))
"""
        base.MeasureAvgClass.__init__(self, Sys, "Rg", StepFreq, CycleFreq, Active,
                                      Source = Source)
        self.MolInd = MolInd
class RgEnsemble(base.MeasureAvgClass): 
    
    ModuleVars = base.MeasureAvgClass.ModuleVars + ["MolIndices"] + ["NMolsRg"]
    
    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True, MolIndices = [None],idx=0):
        """Class for computing radius of gyration."""   
        prev_idx = []
        for m in Sys.Measures:
            if isinstance(m,RgEnsemble):
                prev_idx.append(m.idx)
        if len(prev_idx)>0:
            idx = max(prev_idx)+1 
        Source = """
>>> defs
float PosSum(Dim)
float PosSumSq(Dim)
float RgN
float LastPos(Dim)
float ThisPos(Dim)
float dPos(Dim)
int ind
int jind 
int mRgind
int NMolsRg{i} = {}
int MolIndices{i}(NMolsRg{i}) = {}


>>> init
!print *, "Molecule Indices"
!print *, MolIndices{i}
DO jind = 0, NMolsRg{i} - 1
    mRgind = MolIndices{i}(jind)
    !print *, "MolIndex"
    !print *, mRgind
    !check mRgind
    if (mRgind < 0 .or. mRgind >= NMol) then
        print *, "Invalid molecule index specified in RGEnsemble measure"
        stop
    endif
    !initialize averages
    PosSum = 0.d0
    PosSumSq = 0.d0
    LastPos = 0.d0
    DO ind = MolRange(mRgind), MolRange(mRgind+1)-1
        ThisPos = Pos(ind,:)
        !minimum image, if periodic
        if (DoMinImage) then
            !compute distance from last particle
            dPos = ThisPos - LastPos
            !minimum-image
            dPos = dPos - BoxL * dnint(dPos * iBoxL)
            !compute minimum-imaged position
            ThisPos = LastPos + dPos
        endif
        !add to running sums
        PosSum = PosSum + ThisPos
        PosSumSq = PosSumSq + ThisPos*ThisPos
        !update last position
        LastPos = ThisPos
    END DO
    !compute Rg
    RgN = MolRange(mRgind+1) - MolRange(mRgind)
    Val = Val + sqrt( sum(PosSumSq / RgN - (PosSum*PosSum) / RgN**2))
END DO
Val = Val / NMolsRg{i}
[DOMEASURE]
""".format(len(MolIndices),MolIndices,i=idx)
        base.MeasureAvgClass.__init__(self, Sys, "RgEnsemble", StepFreq, CycleFreq, Active,
                                      Source = Source)
        self.MolIndices = MolIndices
        self.NMolsRg = len(MolIndices)
        self.idx = idx

