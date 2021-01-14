#/usr/bin/env python


### Monte Carlo addition/deletion routines in SIM suite.
### coded by MSS

import numpy as np
from mcmoves import MCMoveClass


class AddDeleteMolecules(MCMoveClass):
    
    Source = """
>>> defs
int i
int istart
int istop
int ind
int ThisMID
int ThisAtomsPerMol
float NewPos(Dim)
float CentroidPos(Dim)
float OldE
float NewE
float r
float lnP
float Pacc0
float RotMat(Dim,Dim)
float ProbIns
float ProbDel
bool Acc
bool DoInsert


>>> externals
ran2
ran2array
ran2int
randomrotmat3D

>>> Main
!choose whether to insert or delete
call ran2(r)
DoInsert = (r > 0.5d0)

!update the attempt
NAtt = NAtt + 1
   
!do an insertion
if (DoInsert) then

    !make sure there are available molecules to insert
    if (NInsertMols > 0) then
       
        call ran2int(NInsertMols, ind)
        TargetMol = MolMoves(ind)
        ThisMID = MolID(TargetMol)
        ThisAtomsPerMol = AtomsPerMol(ThisMID)
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1
        
        NAttMID(ThisMID) = NAttMID(ThisMID) + 1
        OldE = PEnergy
        [external:saveenergystate(Mode = 2)]
   
        !pick a random location 
        call ran2array(Dim, NewPos)
        NewPos = (NewPos - 0.5d0) * BoxL
        if (ThisAtomsPerMol > 1) then
            !pick a random rotation
            CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
            call RandomRotMat3D(3.1415926535897931D0, RotMat)
            do i = istart, istop
                Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + NewPos
            enddo
        else
            Pos(istart,:) = NewPos
        endif
        
        !update energy for adding this molecule
        [external:calcenergyforces(Mode = +2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
        NewE = PEnergy
        
        lnP = Beta * (OldE - NewE) + Beta * MuSet(ThisMID)
        ProbIns = dble(NInactiveMID(ThisMID)) / dble(NInsertMols) / product(BoxL)
        ProbDel = 1.d0 / dble(NDeleteMols + 1)
        lnP = lnP + log(ProbDel / ProbIns)
        Pacc0 = exp(min(0.d0, lnP))
        lnP = lnP + BoltzWeights(ThisMID, NActiveMID(ThisMID)+1) - BoltzWeights(ThisMID, NActiveMID(ThisMID))
        
        !update transition matrix
        TM(NActiveMol, 1) = TM(NActiveMol, 1) + (1.d0 - Pacc0)
        TM(NActiveMol, 2) = TM(NActiveMol, 2) + Pacc0
        
        if (lnP >= 0) then
            Acc = .true.
        else
            call ran2(r)
            Acc = (exp(lnP) > r)
        endif

        if (Acc) then
            NAcc = NAcc + 1
            NActiveMol = NActiveMol + 1
            NActiveMID(ThisMID) = NActiveMID(ThisMID) + 1
            NInactiveMol = NInactiveMol - 1
            NInactiveMID(ThisMID) = NInactiveMID(ThisMID) - 1
            MolActive(TargetMol) = 1
            NAccMID(ThisMID) = NAccMID(ThisMID) + 1
            !swap the indices of the inserted atom the last possible inserted
            MolMoves(ind) = MolMoves(NInsertMols - 1)
            MolMoves(NInsertMols - 1) = TargetMol
            NInsertMols = NInsertMols - 1
            NDeleteMols = NDeleteMols + 1
        else
            MolActive(TargetMol) = -1
            [external:revertenergystate(Mode = 2)]
        endif 
               
    else
    
        !update transition matrix
        TM(NActiveMol, 1) = TM(NActiveMol, 1) + 1.d0
        
    endif
    
!do deletion
else

    !make sure there are available molecules to delete
    if (NDeleteMols > 0) then
    
        call ran2int(NDeleteMols, ind)
        ind = ind + NInsertMols
        TargetMol = MolMoves(ind)
        ThisMID = MolID(TargetMol)
        istart = MolRange(TargetMol)
        istop = MolRange(TargetMol+1) - 1           
        
        NAttMID(ThisMID) = NAttMID(ThisMID) + 1
        OldE = PEnergy
        [external:saveenergystate(Mode = 2)]
        
        !update energy for deleting this molecule
        [external:calcenergyforces(Mode = -2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]           
        NewE = PEnergy
        
        lnP = Beta * (OldE - NewE) - Beta * MuSet(ThisMID)
        ProbDel = 1.d0 / dble(NDeleteMols)
        ProbIns = dble(NInactiveMID(ThisMID) + 1) / dble(NInsertMols+1) / product(BoxL)
        lnP = lnP + log(ProbIns / ProbDel)
        Pacc0 = exp(min(0.d0, lnP))
        lnP = lnP + BoltzWeights(ThisMID, NActiveMID(ThisMID)-1) - BoltzWeights(ThisMID, NActiveMID(ThisMID))
        
        !update transition matrix
        TM(NActiveMol, 1) = TM(NActiveMol, 1) + (1.d0 - Pacc0)
        TM(NActiveMol, 0) = TM(NActiveMol, 0) + Pacc0        

        if (lnP >= 0) then
            Acc = .true.
        else
            call ran2(r)
            Acc = (exp(lnP) > r)
        endif

        if (Acc) then
            NAcc = NAcc + 1
            NActiveMol = NActiveMol - 1
            NActiveMID(ThisMID) = NActiveMID(ThisMID) - 1
            NInactiveMol = NInactiveMol + 1
            NInactiveMID(ThisMID) = NInactiveMID(ThisMID) + 1
            MolActive(TargetMol) = -1
            NAccMID(ThisMID) = NAccMID(ThisMID) + 1
            !swap the indices of the deleted atom and the first possible deleted
            MolMoves(ind) = MolMoves(NInsertMols)
            MolMoves(NInsertMols) = TargetMol
            NInsertMols = NInsertMols + 1
            NDeleteMols = NDeleteMols - 1
        else
            [external:revertenergystate(Mode = 2)]
        endif   
        
    else
    
        !update transition matrix
        TM(NActiveMol, 1) = TM(NActiveMol, 1) + 1.d0
        
    endif
    
endif

!update weights
do i = 0, NMID-1
    BoltzWeights(i, NActiveMID(i)) = BoltzWeights(i, NActiveMID(i)) + MFactor
enddo    
"""    

    ModuleVars = MCMoveClass.ModuleVars + ["NInsertMols", "NDeleteMols", "MolMoves", 
                                           "BoltzWeights", "NAccMID", "NAttMID", "TM",
                                           "MFactor"]
    Name = "constant temperature & chemical potential addition and deletion of molecules"
    
    def __init__(self, Sys):
        """Initializes a class for particle additions and deletions at const T."""
        MCMoveClass.__init__(self, Sys)
        self.UseMID = np.array([True for MolType in Sys.World])
        self.Weight = 0.
        self.MFactor = 0.
        
    def Cleanup(self):
        """Runs when object is eliminated."""
        self.Sys = None
        
    def Reset(self):
        """Resets counters."""
        self.NAcc = 0.
        self.NAtt = 0.     
        #these are for rotational updates
        self.NAccMID = np.zeros_like(self.Sys.NActiveMID)
        self.NAttMID = np.zeros_like(self.Sys.NActiveMID)     
        
    def PreLoad(self):
        """Runs before fortran compilation."""
        self.MolMoves = np.zeros(self.Sys.NMol, dtype = int)
        self.NDeleteMols = 0
        self.NInsertMols = 0
        self.MolMoves = np.zeros(self.Sys.NMol, dtype = int)
        self.BoltzWeights = np.zeros((self.Sys.World.NMID, self.Sys.NMol+1), dtype = float)
        self.TM = np.zeros((self.Sys.NMol+1, 3), dtype=float)
        MCMoveClass.PreLoad(self)
    
    def Init(self):
        """Runs at the start of a cycle."""
        #check
        for (ThisUseMID, MolType) in zip(self.UseMID, self.Sys.World):
            if len(MolType) == 1 or MolType.Rigid or not ThisUseMID: continue
            raise ValueError("This move cannot treat nonrigid polyatomic molecules (MID=%d)." % MolType.MID) 
        #make a list of insertable/deletable molecules
        InsMolMoves = [m.MInd for m in self.Sys.Mol if self.UseMID[m.MID] and self.Sys.MolActive[m.MInd]==-1]
        DelMolMoves = [m.MInd for m in self.Sys.Mol if self.UseMID[m.MID] and self.Sys.MolActive[m.MInd]==+1]
        m = np.array(InsMolMoves + DelMolMoves)
        self.MolMoves[:len(m)] = m
        self.NDeleteMols = len(DelMolMoves)
        self.NInsertMols = len(InsMolMoves)
        if not len(self.MolMoves):
            raise ValueError("Did not find any molecules to add or delete.")

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

            
        
    
        