#/usr/bin/env python


### Monte Carlo volume scaling routines in SIM suite.
### coded by MSS

import numpy as np
from mcmoves import MCMoveClass


class ScaleVolume(MCMoveClass):
    
    Source = """
>>> defs
int i
int istart
int istop
int ind
int j
int m
int NActiveAxes
float DeltaPos(Dim)
float CentroidPos(Dim)
float DeltaVol
float ScaleFactor(Dim)
float OldVol
float NewVol
float OldE
float NewE
float r
float lnP
float OldBoxL(Dim)
bool Acc
bool AxisMask(Dim)

>>> externals
ran2

>>> Main
!update the attempt
NAtt = NAtt + 1

NActiveAxes = count(UseAxis)
OldBoxL = BoxL
OldVol = product(BoxL)
OldE = PEnergy
[external:saveenergystate(Mode = 0)]

!check if we need to scale independently
if (.not. DoIsotropic) then
    !find a random active axis
    call ran2int(NActiveAxes, i)
    j = -1
    do ind = 0, Dim - 1
        if (UseAxis(ind)) j = j + 1
        if (j == i) exit
    enddo
    AxisMask = .false.
    AxisMask(ind) = .true.
else
    AxisMask = UseAxis
endif

!choose a random volume change
call ran2(r)
DeltaVol = Delta * (2.d0 * r - 1.d0)
NewVol = OldVol + DeltaVol

if (NewVol > 0. .and. NewVol >= MinVol .and. NewVol <= MaxVol) then
   
    ScaleFactor = merge((NewVol / OldVol)**(1.d0 / dble(count(AxisMask))), 1.d0, AxisMask)
    BoxL = BoxL * ScaleFactor
    OldPos = Pos
        
    !now scale the molecule centers of mass
    do m = 0, NMol - 1
        !skip frozen and inactive mols
        if (MolActive(m) < 1) cycle
        
        !find the current centroid
        istart = MolRange(m)
        istop = MolRange(m+1) - 1
        CentroidPos = sum(Pos(istart:istop,:), dim=1) / AtomsPerMol(MolID(m))
        
        !find displacement
        DeltaPos = CentroidPos * (ScaleFactor - 1.d0)
        
        !update atom positions
        do i = istart, istop
            Pos(i,:) = Pos(i,:) + DeltaPos   
        enddo
    enddo
           
    !update energy 
    [external:calcenergyforces(Mode = 0, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
    NewE = PEnergy
    
    lnP = Beta * (OldE - NewE) - Beta * PresSet * DeltaVol
    lnP = lnP + NActiveMol * log(NewVol / OldVol)
    if (lnP >= 0) then
        Acc = .true.
    else
        call ran2(r)
        Acc = (exp(lnP) > r)
    endif
    
    if (Acc) then
        NAcc = NAcc + 1
    else
        BoxL = OldBoxL
        Pos = OldPos
        [external:revertenergystate(Mode = 0)]
    endif     
endif
"""    

    ModuleVars = MCMoveClass.ModuleVars + ["UseAxis", "DoIsotropic", "Delta", 
                                           "MaxVol", "MinVol"]
    Name = "constant temperature & pressure linear volume scaling moves"
    
    def __init__(self, Sys):
        """Initializes a class for particle additions and deletions at const T."""
        MCMoveClass.__init__(self, Sys)
        self.UseAxis = np.ones(self.Sys.Dim, dtype=bool)
        self.DoIsotropic = True
        self.Weight = 0.
        self.Delta = 10. * Sys.Units.Angstrom**self.Sys.Dim
        self.MinDelta = 0.
        self.MaxDelta = 1.e300
        self.TargetAcc = 0.7
        self.MaxVol = 1.e300
        self.MinVol = 0.
        
    def Init(self):
        """Runs at the start of a cycle."""
        #set the minimum volume to be the case when the minimum box length
        #is twice the biggest cutoff
        BoxL = [x for x in self.Sys.BoxL if x > 0.]
        minBoxL = min(BoxL)
        maxCut = 0.
        for P in self.Sys.ForceField:
            if not P.Cut is None:
                maxCut = max(P.Cut, maxCut)
        CurVol = np.prod(self.Sys.BoxL)
        self.MinVol = 1.0000001 * CurVol * (2. * maxCut / minBoxL)**self.Sys.World.Dim
        
    def Cleanup(self):
        """Runs when object is eliminated."""
        self.Sys = None

    def UpdateDelta(self):
        """Updates the maximum move to achieve target acceptance."""
        if self.NAtt > 0:
            self.Delta *= (self.NAcc / self.NAtt + 0.1) / (self.TargetAcc + 0.1)
            self.Delta = min(self.MaxDelta, max(self.MinDelta, self.Delta))   
            
    def AccFrac(self):
        """Returns the acceptance fraction(s)."""
        acc1 = self.NAcc / (self.NAtt +1.e-8)
        return [acc1]
        

            
        
    
        