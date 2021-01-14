#/usr/bin/env python


### Monte Carlo semigrand mutation (as addition/deletion) routines in SIM suite.
### coded by MSS

import numpy as np
from mcmoves import MCMoveClass


class MutateMolecules(MCMoveClass):
    
    Source = """
>>> defs
int i
int istart
int istop
int ind
int DeleteMID
int DeleteMol
int DeleteAtomsPerMol
int InsertMol
int InsertMID
int InsertAtomsPerMol
int Delta1
int N1
float NewPos(Dim)
float DeleteCentroidPos(Dim)
float InsertCentroidPos(Dim)
float OldE
float NewE
float r
float lnP
float Pacc0
float RotMat(Dim,Dim)
bool Acc

>>> externals
ran2
ran2int
randomrotmat3D

>>> Main
!update the attempt
NAtt = NAtt + 1

!pick a random molecule to mutate that's active and type MID1 or MID2
call ran2int(NActiveMID(MID1) + NActiveMID(MID2), ind)
ActiveInd = -1
DeleteMol = -1
do while (ActiveInd < ind)
    DeleteMol = DeleteMol + 1
    DeleteMID = MolID(DeleteMol)
    if (MolActive(DeleteMol) == 1 .and. (DeleteMID == MID1 .or. DeleteMID == MID2)) then
        ActiveInd = ActiveInd + 1
    endif
enddo
DeleteAtomsPerMol = AtomsPerMol(DeleteMID)

!find the molecule type to replace with
if (DeleteMID == MID1 .and. NInactiveMID(MID2) > 0) then
    InsertMID = MID2
    Delta1 = -1
elseif (DeleteMID == MID2 .and. NInactiveMID(MID1) > 0) then
    InsertMID = MID1
    Delta1 = 1
else
    InsertMID = -1
endif

!order parameter for weights and TM
N1 = NActiveMID(MID1)

if (InsertMID >= 0) then
    
    !find a replacement molecule
    call ran2int(NInactiveMID(InsertMID), ind)
    ActiveInd = -1
    InsertMol = -1
    do while (ActiveInd < ind)
        InsertMol = InsertMol + 1
        if (MolActive(InsertMol) == -1 .and. MolID(InsertMol) == InsertMID) then
            ActiveInd = ActiveInd + 1
        endif
    enddo 
    InsertAtomsPerMol = AtomsPerMol(InsertMID)
    
    !save current energy
    OldE = PEnergy
    [external:saveenergystate(Mode = 0)]
    
    !get the center of mass of the original molecule
    istart = MolRange(DeleteMol)
    istop = MolRange(DeleteMol+1) - 1
    DeleteCentroidPos = sum(Pos(istart:istop,:), dim=1) / DeleteAtomsPerMol
        
    !update energy for deleting the original
    TargetMol = DeleteMol
    MolActive(DeleteMol) = -1
    [external:calcenergyforces(Mode = -2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)] 

    !now add the new molecule at the same location
    istart = MolRange(InsertMol)
    istop = MolRange(InsertMol+1) - 1
    InsertCentroidPos = sum(Pos(istart:istop,:), dim=1) / InsertAtomsPerMol

    if (InsertAtomsPerMol > 1) then
        !pick a random rotation
        call RandomRotMat3D(3.1415926535897931D0, RotMat)
        do i = istart, istop
            Pos(i,:) = matmul(RotMat, Pos(i,:) - InsertCentroidPos) + DeleteCentroidPos
        enddo
    else
        Pos(istart,:) = DeleteCentroidPos
    endif
    
    !update energy for adding this molecule
    TargetMol = InsertMol
    MolActive(InsertMol) = 1
    [external:calcenergyforces(Mode = +2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
    
    !get the new energy
    NewE = PEnergy
    
    lnP = Beta * (OldE - NewE) + Beta * (MuSet(InsertMID) - MuSet(DeleteMID))
    Pacc0 = exp(min(0.d0, lnP))        
    lnP = lnP + BoltzWeights(N1 + Delta1) - BoltzWeights(N1)
    
    !update transition matrix
    TM(N1, 1) = TM(N1, 1) + (1.d0 - Pacc0)
    TM(N1, 1+Delta1) = TM(N1, 1+Delta1) + Pacc0
    
    if (lnP >= 0) then
        Acc = .true.
    else
        call ran2(r)
        Acc = (exp(lnP) > r)
    endif

    if (Acc) then
        NAcc = NAcc + 1
        NActiveMID(MID1) = NActiveMID(MID1) + Delta1
        NActiveMID(MID2) = NActiveMID(MID2) - Delta1
        NInactiveMID(MID1) = NInactiveMID(MID1) - Delta1
        NInactiveMID(MID2) = NInactiveMID(MID2) + Delta1
        N1 = N1 + Delta1
        BoltzWeights(N1) = BoltzWeights(N1) + MFactor
    else
        MolActive(InsertMol) = -1
        MolActive(DeleteMol) = 1
        BoltzWeights(N1) = BoltzWeights(N1) + MFactor
        [external:revertenergystate(Mode = 0)]
    endif   
    
else

    !update transition matrix
    TM(N1, 1) = TM(N1, 1) + 1.d0
    
endif
"""    

    ModuleVars = MCMoveClass.ModuleVars + ["MID1", "MID2", "BoltzWeights", "TM", "MFactor"]
    Name = "constant temperature & chemical potential addition and deletion of molecules"
    
    def __init__(self, Sys):
        """Initializes a class for particle mutations at const T."""
        MCMoveClass.__init__(self, Sys)
        self.Weight = 0.
        self.MID1 = 0
        self.MID2 = 0
        self.MFactor = 0.
        
    def Cleanup(self):
        """Runs when object is eliminated."""
        self.Sys = None
        
    def Reset(self):
        """Resets counters."""
        self.NAcc = 0.
        self.NAtt = 0.      
        
    def PreLoad(self):
        """Runs before fortran compilation."""
        n = self.Sys.NActiveMID[self.MID1] + self.Sys.NInactiveMID[self.MID1]
        self.BoltzWeights = np.zeros((n+1), dtype = float)
        self.TM = np.zeros((n+1, 3), dtype=float)
        MCMoveClass.PreLoad(self)
    
    def Init(self):
        """Runs at the start of a cycle."""
        #check
        for MID in [self.MID1, self.MID2]:
            MolType = self.Sys.World[MID]
            if len(MolType) == 1 or MolType.Rigid: continue
            raise ValueError("This move cannot treat nonrigid polyatomic molecules (MID=%d)." % MID) 
        #check
        NAct = self.Sys.NActiveMID
        NIna = self.Sys.NInactiveMID
        if self.MID1 == self.MID2:
            raise ValueError("Both species in semigrand move are the same (MID = %d)" % self.MID1)
        if NAct[self.MID1] + NIna[self.MID1] != NAct[self.MID2] + NIna[self.MID2]:
            raise ValueError("Semigrand move for MID = %d and %d have different total num of molecules." % (self.MID1, self.MID2))
        if NAct[self.MID1] != NIna[self.MID2] or NAct[self.MID2] != NIna[self.MID1]:
            raise ValueError("Semigrand move requires equal and opposite active/inactive for MID = %d and %d." % (self.MID1, self.MID2))    
            
            
    def ResetTM(self):
        """Resets the transition matrix."""
        self.TM.fill(0.)
        
    def CalcTMWeights(self):
        """Uses the transition matrix to estimate weights."""
        n = len(self.TM)
        Weights = np.zeros(n, dtype=float)
        TM = self.TM
        TMSum = self.TM.sum(axis=1)
        for i in range(1, n):
            if (TM[i,0] <= 0. or TM[i-1,2] <= 0):
                Weights[i] = Weights[i-1]
            else:
                Weights[i] = Weights[i-1] + np.log(TM[i-1,2] * TMSum[i] / (TM[i,0] * TMSum[i-1]))
        return Weights
    
    def AccFrac(self):
        """Returns the acceptance fraction."""
        acc1 = self.NAcc / (self.NAtt +1.e-8)
        return [acc1]
              
        

            
        
    
        