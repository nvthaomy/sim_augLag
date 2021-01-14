#/usr/bin/env python


### Monte Carlo integration routine objects in SIM suite.
### coded by MSS

import numpy as np
import sim.fortran as fortran



class MCMoveClass(object):
    
    ModuleVars = ["NAtt", "NAcc"]
    ModuleConsts = []
    
    Source = ""

    def __init__(self, Sys, Name = ""):
        """Initializes a MC move; to be modified for derived classes."""
        object.__setattr__(self, "LibVars", {})
        self.Weight = 1.
        self.Sys = Sys
        self.Name = Name
        self.FortCode = None
        #reset counters
        self.Reset()
        
    def Reset(self):
        """Resets counters."""
        self.NAcc = 0.
        self.NAtt = 0.        

    def PreLoad(self):
        """Run before fortran compilation."""
        self.Sys.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts, self.FortPrefix)
        for (Var, ModuleVar) in self.LibVars0.items():
            self.Source = fortran.ReplaceToken(self.Source, Var, ModuleVar)
        self.FortCode = fortran.FortCode(Source = self.Source, Comment = self.Name)
        
    def PostLoad(self):
        """Run after fortran compilation."""
        self.Sys.Lib.VarPostLoad(self)     
    
    def Cleanup(self):
        """Runs when object is eliminated."""
        self.Sys = None   
        
    def Init(self):
        """Runs at the start of a cycle."""
        pass
    
    def Finalize(self):
        """Runs at the end of a cycle."""
        pass    
    
    def __setattr__(self, name, val):
        Var = self.LibVars.get(name, None)
        if Var is None:
            v = self.__dict__.get(name, None)
            if isinstance(v, np.ndarray):
                v[:] = val
            else:
                object.__setattr__(self, name, val)   
        else:
            if Var.shape:
                Var[:] = val
            else:
                Var.itemset(val)
                    
    def __getattr__(self, name):
        Var = self.LibVars.get(name, None)
        if Var is None:
            raise AttributeError(name) 
        else:
            if not Var.shape: Var = Var.item()
            return Var
    

class CanonicalDisplacements(MCMoveClass):
    
    Source = """   
>>> defs
int i
int istart
int istop
int m
int ActiveInd
int ind
int ThisAtomsPerMol
int Step
float OldAtomPos(Dim)
float dPos(Dim)
float CentroidPos(Dim)
float OldE
float NewE
float r
float lnP
float RotMat(Dim,Dim)
bool Acc
bool Rotate 

>>> externals
ran2
ran2array
ran2int
randomrotmatxyz

>>> Main
!check if active
if (NActiveMol == 0) then
    NAtt = NAtt + 1
    return
endif

!pick a random molecule that's active
call ran2int(NActiveMol, ind)
ActiveInd = -1
TargetMol = -1
do while (ActiveInd < ind)
    TargetMol = TargetMol + 1
    if (MolActive(TargetMol) == 1) ActiveInd = ActiveInd + 1
enddo
ThisAtomsPerMol = AtomsPerMol(MolID(TargetMol))
 

!check if this is a rigid molecule or not
if (MolIsRigid(TargetMol)) then

    !rigid molecule; find start and stop atoms
    istart = MolRange(TargetMol)
    istop = MolRange(TargetMol+1) - 1
    
    !save the old information
    OldE = PEnergy
    OldPos(istart:istop,:) = Pos(istart:istop,:)
    [external:saveenergystate(Mode = 2)]
    
    !calculate energy without molecule
    [external:calcenergyforces(Mode = -2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
    
    !decide whether to move or rotate
    call ran2(r)
    Rotate = (r > 0.5d0)
    
    if (Rotate) then
            
        NAtt2 = NAtt2 + 1
        !rotate about a random axis 
        call RandomRotMatxyz(Delta2, RotMat)
        CentroidPos = sum(Pos(istart:istop,:), dim=1) / ThisAtomsPerMol
        do i = istart, istop
            Pos(i,:) = matmul(RotMat, Pos(i,:) - CentroidPos) + CentroidPos
        enddo
        
    else
            
        NAtt = NAtt + 1        
        !translation
        call ran2array(Dim, dPos)
        dPos = Delta * (2.d0 * dPos - 1.d0)
        do i = istart, istop
            Pos(i,:) = Pos(i,:) + dPos
        enddo
        
    endif
    
    !calculate energy with molecule
    [external:calcenergyforces(Mode = +2, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
    NewE = PEnergy

    lnP = Beta * (OldE - NewE)
    if (lnP >= 0) then
        Acc = .true.
    else
        call ran2(r)
        Acc = (exp(lnP) > r)
    endif
    if (Acc) then
        if (Rotate) then
            NAcc2 = NAcc2 + 1
        else
            NAcc = NAcc + 1
        endif
    else
        !revert to old state
        Pos(istart:istop,:) = OldPos(istart:istop,:)
        [external:revertenergystate(Mode = 2)]
    endif     
        
else

    !this is a non-rigid molecule; pick a random atom to displace
    call ran2int(ThisAtomsPerMol, i)
    TargetAtom = i + MolRange(TargetMol)

    !save the old information
    OldE = PEnergy
    OldAtomPos = Pos(TargetAtom,:)
    [external:saveenergystate(Mode = 1)]
    
    !calculate energy without atom
    [external:calcenergyforces(Mode = -1, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]

    NAtt = NAtt + 1
    
    !random displacement
    call ran2array(Dim, dPos)
    dPos = Delta * (2.d0 * dPos - 1.d0)
    Pos(TargetAtom,:) = Pos(TargetAtom,:) + dPos
    
    !calculate energy with atom
    [external:calcenergyforces(Mode = +1, CalcForce = .false., CalcVirial = .false., CalcDUParam = .false., CalcDWParam = .false.)]
    NewE = PEnergy

    lnP = Beta * (OldE - NewE)
    if (lnP >= 0) then
        Acc = .true.
    else
        call ran2(r)
        Acc = (exp(lnP) > r)
    endif
    if (Acc) then
        NAcc = NAcc + 1
    else
        !revert to old state
        Pos(TargetAtom,:) = OldAtomPos
        [external:revertenergystate(Mode = 1)]
    endif
    
endif
"""    

    ModuleVars = MCMoveClass.ModuleVars + ["NAtt2", "NAcc2", "Delta", "Delta2"]
    Name = "canonical molecule displacements and rotations."
    
    def __init__(self, Sys):
        """Initializes a class for single-particle displacements at const T."""
        MCMoveClass.__init__(self, Sys)
        self.Delta = Sys.Units.Angstrom * 0.1
        self.Delta2 = np.pi / 18.
        self.MinDelta = 0.
        self.MaxDelta = 1.e300
        self.TargetAcc = 0.3
        self.MinDelta2 = 0.
        self.MaxDelta2 = np.pi
        self.TargetAcc2 = 0.3  

    def PreLoad(self):
        """Run before fortran compilation."""
        if not np.any(abs(self.Sys.MolActive) == 1):
            raise ValueError("Did not find any molecules to move.")  

        MCMoveClass.PreLoad(self)
        
    def PostLoad(self):
        """Run after fortran compilation."""
        MCMoveClass.PostLoad(self)
        
    def Cleanup(self):
        """Runs when object is eliminated."""
        self.Sub0 = None
        self.Sys = None
        
    def Reset(self):
        """Resets counters."""
        self.NAcc = 0.
        self.NAtt = 0.     
        #these are for rotational updates
        self.NAcc2 = 0.
        self.NAtt2 = 0.   

    def UpdateDelta(self):
        """Updates the maximum move to achieve target acceptance."""
        if self.NAtt > 0:
            acc1 = (self.NAcc / (self.NAtt + 1.e-8)) 
            self.Delta *= (acc1 + 5.) / (self.TargetAcc + 5.)
            self.Delta = min(self.MaxDelta, max(self.MinDelta, self.Delta))
        if self.NAtt2 > 0:
            acc2 = (self.NAcc2 / (self.NAtt2 + 1.e-8))
            self.Delta2 *= (acc2 + 5.) / (self.TargetAcc2 + 5.)
            self.Delta2 = min(self.MaxDelta2, max(self.MinDelta2, self.Delta2))       
            
    def AccFrac(self):
        """Returns the acceptance fractions."""
        acc1 = self.NAcc / (self.NAtt +1.e-8)
        acc2 = self.NAcc2 / (self.NAtt2 + 1.e-8)
        return [acc1, acc2]

        

        

            
        
    
        