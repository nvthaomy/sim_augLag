#/usr/bin/env python


### Class definitions for systems in SIM suite.
### coded by MSS

import numpy as np


import flags
import positions
import velocities

#must import visualize manually
#import visualize

import sim.chem as chem
import sim.geom as geom



ForceCompile = True




#define an iterator over center of masses
def RigidUnits(Sys):
    """An iterator over rigid units."""
    for Mol in Sys.Mol:
        if Mol.HasRigid:
            yield Mol 
        else:
            for Atom in Mol:
                yield Atom



class System(object):
       
    ModuleVars = ["BoxL", "AIDCount", "MolRange", "MolID", "NDOF", "COMPos",
                  "Pos", "Vel", "Force", "Mass", "iMass", "sqrtMass", "MInd", 
                  "NActiveMol", "NInactiveMol", "NActiveMID", "NInactiveMID", 
                  "AID", "SID", "MID", "MolActive", "PEnergy", "KEnergy", "TEnergy",
                  "Virial", "TempSet", "PresSet", "MuSet", "TargetAtom", "TargetMol"]

    ModuleConsts = ["NAtom", "NMol", "NDOFMID", "AtomsPerMol", "MaxAtomsPerMol", "MaxRigidAtoms", 
                    "MolIsRigid", "RBondLengthSq", "RBondInd", "RBondRange"]
    
    UpdateActiveSource = """
>>> defs
int i
int j
int Start
int Stop
int ActiveState
int ThisMolID
int ThisAID

>>> Main
NDOF = 0
AIDCount = 0
NActiveMol = 0
NInactiveMol = 0
NActiveMID = 0
NInactiveMID = 0
do i = 0, NMol-1
    Start = MolRange(i) 
    Stop = MolRange(i+1) - 1
    if (MolActive(i) == 1) then
        NActiveMol = NActiveMol + 1
        ThisMolID = MolID(i)
        NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
        NDOF = NDOF + NDOFMID(ThisMolID)
        do j = Start, Stop
            ThisAID = AID(j)
            AIDCount(ThisAID) = AIDCount(ThisAID) + 1
        enddo
    elseif (MolActive(i) == -1) then
        NInactiveMol = NInactiveMol + 1
        ThisMolID = MolID(i)
        NInactiveMID(ThisMolID) = NInactiveMId(ThisMolID) + 1
    endif
enddo
"""    

    HideMolSource = """
>>> inputmaps
int MolInd = MolInd

>>> defs
int i
int ThisMolID
int ThisAID

>>> Main
if (MolActive(MolInd) == 1) then
    ThisMolID = MID(MolInd)
    MolActive(MolInd) = -1
    NActiveMol = NActiveMol - 1
    NInactiveMol = NInactiveMol + 1
    NActiveMID(ThisMolID) = NActiveMID(ThisMolID) - 1
    NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) + 1
    NDOF = NDOF - NDOFMID(ThisMolID)
    do i = MolRange(MolInd), MolRange(MolInd + 1) - 1
        ThisAID = AID(i)
        AIDCount(ThisAID) = AIDCount(ThisAID) - 1
    enddo
endif
"""

    ShowMolSource = """
>>> inputmaps
int MolInd = MolInd

>>> defs
int i
int ThisMolID
int ThisAID

>>> Main
if (MolActive(MolInd) == -1) then
    ThisMolID = MID(MolInd)
    MolActive(MolInd) = 1
    NActiveMol = NActiveMol + 1
    NInactiveMol = NInactiveMol - 1
    NActiveMID(ThisMolID) = NActiveMID(ThisMolID) + 1
    NInactiveMID(ThisMolID) = NInactiveMID(ThisMolID) - 1   
    NDOF = NDOF + NDOFMID(ThisMolID)
    do i = MolRange(MolInd), MolRange(MolInd+1) - 1
        ThisAID = AID(i)
        AIDCount(ThisAID) = AIDCount(ThisAID) + 1
    enddo
endif
"""   
        
    def __init__(self, World, Name = "sys", BoxL = None):
        """Creates a new system.
Arguments:
    World: a chemical world from which this system is built
    Name: string name of system
"""
        import sim.fortran as fortran
        import sim.measure as measure
        import sim.integrate as integrate
        import sim.potential as potential
        
        #add this dictionary first since it affects setattr routine
        object.__setattr__(self, "LibVars", {})
        
        self.Loaded = False
        
        #load the world if not already loaded
        if not World.Locked:
            World.Load()
            
        #make references
        self.World = World
        self.Name = Name

        #check dimensionality
        if any([MolType.HasRigid for MolType in World]) and World.Dim != 3:
            raise ValueError("Rigid bonds currently only supported for Dim=3.")
              
        #flags (options)
        self.Flags = flags.FlagsClass()

        #units
        self.Units = World.Units  
       
        #dimensionality
        self.Dim = World.Dim
        
        #box side lengths
        #have to add this variable using the parent class b/c setattr is overrided below
        object.__setattr__(self, "BoxL", np.zeros(self.Dim, float))
        if not BoxL is None: self.BoxL = BoxL
        
        #fortran module
        self.Lib = fortran.Module(Name, KeepSource = True, UseFortModule = True)
        #by default, force compilation of the system library each time
        #it is created
        self.Lib.ForceCompile = ForceCompile
        #add base fortran
        fortran.commoncode.AddCommonFortran(self.Lib)
        #add subroutine from this object
        fc = fortran.FortCode(Source = System.UpdateActiveSource)
        self.ActiveSub = fortran.Subroutine("UpdateActive", FortCode = fc, 
                                            SubStr = "self.ActiveSub")
        self.Lib.append(self.ActiveSub)   
        fc = fortran.FortCode(Source = System.HideMolSource)
        self.HideMolSub = fortran.Subroutine("HideMol", FortCode = fc, 
                                            SubStr = "self.HideMolSub")
        self.Lib.append(self.HideMolSub)  
        fc = fortran.FortCode(Source = System.ShowMolSource)
        self.ShowMolSub = fortran.Subroutine("ShowMol", FortCode = fc, 
                                            SubStr = "self.ShowMolSub")
        self.Lib.append(self.ShowMolSub)     
        
        #energies and thermo
        self.PEnergy = 0.
        self.KEnergy = 0.
        self.TEnergy = 0.
        self.Virial = 0.
        self.TempSet = self.Units.TempDefault
        self.PresSet = 0.
        self.MuSet = np.zeros((World.NMID), int)

        #chemistry information
        self.AtomsPerMol = np.array([len(MolType) for MolType in World], dtype=int)
        self.MaxAtomsPerMol = np.max(self.AtomsPerMol)
        MaxN = max([len(MolType) for MolType in World])
        self.COMPos = np.zeros((len(World), MaxN, World.Dim))
        self.MaxRigidAtoms = 0
        for (i, MolType) in enumerate(World):
            if MolType.COMPos is None: continue
            self.COMPos[i,:,:] = MolType.COMPos
            if MolType.Rigid:
                self.MaxRigidAtoms = max(self.MaxRigidAtoms, len(MolType))  
        #lists of rigid bonds
        self.RBondRange = [0]
        self.RBondLengthSq = []
        self.RBondInd = []
        for MolType in World:
            for Bond in MolType.Bonds:
                if Bond.RLength > 0: 
                    self.RBondInd.append([Bond.SType1.AInd, Bond.SType2.AInd])
                    self.RBondLengthSq.append(Bond.RLength**2)
            self.RBondRange.append(len(self.RBondLengthSq))
        #pad the arrays at the end b/c fortran can't do zero-length array constants
        self.RBondInd.append([-1, -1])
        self.RBondLengthSq.append(-1)
        self.RBondRange = np.array(self.RBondRange, int)
        self.RBondLengthSq = np.array(self.RBondLengthSq, float)
        self.RBondInd = np.array(self.RBondInd, int)   
        self.NDOFMID = np.array([chem.NDOF(MolType, World.Dim) for MolType in World], dtype=int)
                        
        #molecule (and atom) lists
        self.Mol = []
        self.Atom = []

        #number of atoms and degrees of freedom
        self.NAtom = 0
        self.NMol = 0
        self.NDOF = 0    

        #counts
        self.NActiveMol = 0
        self.NActiveMID = np.zeros(World.NMID, int)
        self.NInactiveMol = 0
        self.NInactiveMID = np.zeros(World.NMID, int)
        self.AIDCount = np.zeros(World.NAID, int)

        #initialize molecule arrays        
        self.__MolArrayVars = ["MolRange", "MolID", "MolActive", "MolIsRigid"]
        self.MolRange = np.zeros(1, int)
        self.MolID = np.zeros(0, int)
        self.MolActive = np.array([], dtype=int)
        self.MolIsRigid = np.array([], dtype=bool)
        self.NActiveMol = 0
        
        #initialize atom arrays
        self.__AtomArrayVars = []
        Dim = World.Dim
        #initialize atom arrays
        self.AddAtomArray("Pos", Dim, dtype=float)
        self.AddAtomArray("Vel", Dim, dtype=float)
        self.AddAtomArray("Force", Dim, dtype=float)
        self.AddAtomArray('Mass', 1, dtype=float)
        self.AddAtomArray('iMass', 1, dtype=float)
        self.AddAtomArray('sqrtMass', 1, dtype=float)
        self.AddAtomArray('MInd', 1, dtype=int)
        self.AddAtomArray('AID', 1, dtype=int)
        self.AddAtomArray('SID', 1, dtype=int)
        self.AddAtomArray('MID', 1, dtype=int)
        
        #target atoms and molecules for individual energy calcs
        self.TargetAtom = 0
        self.TargetMol = 0

        #force field
        self.ForceField = potential.ForceField(self)
        
        #measures
        self.Measures = measure.MeasuresClass(self)
        measure.standard.Add(self.Measures)
        
        #integrator
        self.Int = integrate.Integrator(self)
        integrate.standard.Add(self.Int)
        
        
    def Load(self):
        """Compiles all Fortran functions in this system."""
        #constants from other objects
        for Parent in [self.World.Units, self.World]:
            self.Lib.VarPreLoad(Parent, Consts = Parent.ModuleConsts)
        #module variables
        self.Lib.VarPreLoad(self, self.ModuleVars, self.ModuleConsts)
        self.Lib.Load()
        self.Lib.VarPostLoad(self)
        self.Loaded = True
        self.UpdateActive()
        
        
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
                    
    def __len__(self):
        return len(self.Mol)


    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        for Mol in self.Mol:
            Mol.Cleanup()
        del self.Mol[:]
        del self.Atom[:]
        self.LibVars.clear()
        self.ForceField.Cleanup()
        self.ForceField = None
        self.Measures.Cleanup()
        self.Measures = None
        self.Int.Cleanup()
        self.Int = None
        self.Flags = None
        self.World = None
        self.Units = None
        self.Lib = None


    def __repr__(self):
        return repr(self.Mol)

        
    def AddAtomArray(self, VarName, Dim = 1, dtype = float, AutoNumber = False):
        """Adds a new variable to the atom arrays."""
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        if AutoNumber:
            i = 0
            while hasattr(self, VarName + "%d" % i):
                i += 1
            VarName = VarName + "%d" % i
        elif hasattr(self, VarName):
            raise ValueError("Atom array %s already exists." % VarName)
        if Dim == 1:
            Var = np.zeros(self.NAtom, dtype=dtype)
        else:
            Var = np.zeros((self.NAtom, Dim), dtype=dtype)
            Var = np.asfortranarray(Var)
        self.__AtomArrayVars.append(VarName)
        object.__setattr__(self, VarName, Var)
        return VarName        
        
        
    def ResizeAtomArrays(self, N):
        """Resizes all atom arrays."""
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        for varname in self.__AtomArrayVars:
            var = getattr(self, varname)
            newshape = (N,) + var.shape[1:]
            newvar = np.resize(var, newshape)
            delattr(self, varname)
            setattr(self, varname, newvar)
        self.NAtom = N


    def AddMol(self, Mol):
        """Adds a molecule."""
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        #check type
        if not self.World.HasMolType(Mol):
            raise TypeError("World does not include moltype %s" % repr(Mol))
        
        M = len(Mol)
    
        #resize atom arrays
        StartInd = len(self.Atom)
        self.ResizeAtomArrays(StartInd + M)
                            
        #add the molecule list
        Mol.MInd = len(self.Mol)
        Mol.Parent = self
        self.Mol.append(Mol)

        #update the atom list
        self.Atom.extend(Mol)
        for Ind in range(StartInd, StartInd + M):
            Atom = self.Atom[Ind]
            Atom.Ind = Ind
            #set id and properties
            self.SID[Ind] = Atom.SID
            self.AID[Ind] = Atom.AID
            self.MID[Ind] = Mol.MID
            self.Mass[Ind] = Atom.Mass
            self.iMass[Ind] = 1. / Atom.Mass
            self.sqrtMass[Ind] = np.sqrt(Atom.Mass)
            self.MInd[Ind] = Mol.MInd

        #update the molecule and atom counts
        self.NActiveMID[Mol.MID] += 1
        self.NActiveMol += 1
        for Atom in Mol:
            self.AIDCount[Atom.AID] += 1

        #update the degrees of freedom
        self.NDOF += chem.NDOF(Mol, self.Dim)
        self.NMol += 1

        #update the molecule arrays
        for varname in self.__MolArrayVars:
            var = getattr(self, varname)
            m = len(var)
            newshape = (m + 1,) + var.shape[1:]
            newvar = np.resize(var, newshape)
            delattr(self, varname)
            setattr(self, varname, newvar)        
        self.MolRange[-1] = StartInd + M
        self.MolID[-1] = Mol.MID
        self.MolActive[-1] = 1
        self.MolIsRigid[-1] = Mol.Rigid


    def __add__(self, Mol):
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        self.AddMol(Mol)
        return self

    def append(self, Mol):
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        self.AddMol(Mol)

    def extend(self, MolList):
        if self.Loaded: raise ValueError("Cannot run this after loading.")
        for Mol in MolList:
            self.AddMol(Mol)                          
        
    def UpdateActive(self):
        """Updates all variables that depend on active atoms or not."""
        exec(self.ActiveSub.CallObj)
            
    def SetActiveState(self, Mol, ActiveState):
        """Changes the active status of a molecule."""
        if type(Mol) is int:
            Mol = self.Mol[Mol]
        self.MolActive[Mol.Ind] = ActiveState
        self.UpdateActive(self)   


    def SetMolPos(self, Mol, Pos, Orient = None):
        """Sets the position of a rigid molecule."""
        if type(Mol) is int:
            Mol = self.Mol[Mol]
        Mol.Pos = Pos.copy()
        if not Orient is None:
            Mol.Orient = Orient.copy()
        MInd = Mol.MInd
        rm = geom.RotMatEuler(Mol.Orient[0], Mol.Orient[1], Mol.Orient[2])
        Pos = np.dot(Mol.MType.COMPos, rm) + Mol.Pos
        self.Pos[self.MolRange[MInd]:self.MolRange[MInd+1],:] = Pos


    def SetPos(self, Unit, Pos, Orient = None):
        """Sets the position of a unit."""
        if isinstance(Unit, chem.Mol):
            self.SetMolPos(Unit, Pos, Orient)
        else:
            Unit.Pos = Pos       


    def ResetPos(self):
        """Resets atoms to original box."""
        BoxL = np.where(self.BoxL == 0, 1.e300, self.BoxL)
        for Unit in RigidUnits(self):
            Pos = Unit.Pos - BoxL * np.round_(Unit.Pos / BoxL)
            self.SetPos(Unit, Pos)


    def RigidUnitsIter(self):
        """Returns an iterator over rigid units."""
        return RigidUnits(self)
        
        
    def ScaleBox(self, BoxL, UpdatePos = True):
        """Scales the box size."""
        if not type(BoxL) in [list, np.ndarray]:
            BoxL = np.array([BoxL]*self.Dim, dtype=float)
        if np.all(BoxL == self.BoxL): 
            return 
        if UpdatePos:
            Scale = BoxL / self.BoxL
            Scale[self.BoxL <= 0.] = 1.
            for Unit in RigidUnits(self):
                NewPos = Unit.Pos * Scale
                self.SetPos(Unit, NewPos)
        self.BoxL = np.where(BoxL <= 0., 0., BoxL)  
        
        
    def MinImage(self, Vec):
        """Returns a minimum imaged vector based on simulation box."""
        Mask = np.nonzero(self.BoxL > 0.)[0]
        ret = np.array(Vec, dtype=float)
        ret[...,Mask] = ret[...,Mask] - self.BoxL[Mask] * np.round_(ret[...,Mask] / self.BoxL[Mask])
        return ret
        
        
    def ApplySettings(self, Settings):
        """Applies settings from Settings dictionary."""
        if not type(Settings) is dict:
            raise TypeError("Settings must be a dictionary.")
        for (k,v) in Settings.items():
            if k == "TempSet":
                self.TempSet = v
            elif k == "PresSet":
                self.PresSet = v
            elif k == "BoxL":
                self.ScaleBox(v)
            else:
                raise ValueError("I do not recognize setting %s" % k)
                
    def Check(self):
        """Checks system to ensure nothing is wrong."""
        self.ForceField.Check()
        
        

