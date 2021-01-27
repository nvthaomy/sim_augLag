#/usr/bin/env python

### Radius of gyration measuring object in SIM suite.
### coded by SPC

import numpy as np

import base


class Rg(base.MeasureAvgClass):
    
    ModuleVars = base.MeasureAvgClass.ModuleVars + ["MolIndices"]

    def __init__(self, Sys, StepFreq = 0, CycleFreq = 1, Active = True, MolIndices = []):
        """Class for computing radius of gyration for specific molecule indices."""   
        Source = """
>>> defs
float PosSum(Dim)
float PosSumSq(Dim)
float RgN
float LastPos(Dim)
float ThisPos(Dim)
float dPos(Dim)
float ThisRg
int ind
int j
int MolInd
int NMol

>>> init
NMol = size(MolIndices)
Val = 0.
do j = 0, NMol-1
    MolInd = MolIndices(j)
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
    ThisRg = sqrt( sum(PosSumSq / RgN - (PosSum*PosSum) / RgN**2))
    Val = Val + ThisRg
enddo
Val = Val / dble(NMol)
[DOMEASURE]
"""
        base.MeasureAvgClass.__init__(self, Sys, "Rg", StepFreq, CycleFreq, Active,
                                      Source = Source)
        self.MolIndices = MolIndices
        
    
    def PreLoad(self):
        """Run before library making."""
        self.MolIndices = np.array(self.MolIndices, dtype=int)
        base.MeasureAvgClass.PreLoad(self)

