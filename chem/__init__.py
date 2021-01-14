#/usr/bin/env python


### Class definitions for atoms, molecules, and world in SIM suite.
### coded by MSS

import numpy as np
import sim.units as units


#maximum bond order to tabulate
BondOrdLimit = 4

#CONVENTIONS:
#AInd - index of an atom within a molecule 
#MInd - index of a molecule within a system
#Ind - absolute index of atom
#AID - type index of atom
#SID - type index of an atom site within a molecule
#MID - type index of molecule



def NDOF(Mol1, Dim):
    """Returns the number of translational degrees of freedom for a molecule."""
    if isinstance(Mol1, Mol):
        MolType1 = Mol1.MType
    elif isinstance(Mol1, MolType):
        MolType1 = Mol1
    else:
        raise TypeError("Expected Mol or MolType")
    n = len(MolType1)
    if n == 1:
        RigidDOF = Dim
    elif Dim == 1:
        RigidDOF = 1
    elif Dim == 2:
        RigidDOF = 3
    elif Dim == 3:
        if n == 2:
            RigidDOF = 5
        else:
            RigidDOF = 6
    else:
        raise ValueError("Don't know how to count DOFs for Dim = %d" % Dim)
    if MolType1.Rigid:
        return RigidDOF
    else:
        NConstr = len([b for b in MolType1.Bonds if b.Rigid])
        return max(RigidDOF, Dim * n - NConstr)


def PairID(ID1, ID2):
    """Returns a unique pair ID."""
    if ID1 < ID2:
        ID1, ID2 = ID2, ID1
    return ID1 * (ID1 + 1) / 2 + ID2


def GetPair(ID):
    """Returns a pair number from an ID"""
    ID1 = int(0.5*(-1+np.sqrt(1+8*ID)) + 0.00001)
    ID2 = ID - ID1*(ID1+1)/2
    return ID1, ID2


def NPairID(NID):
    """Returns total number of unique pair IDs."""
    return NID * (NID + 1) / 2


def NMultiID(NID, Rank):
    """Returns a unique ID for multiple atom IDs."""
    if Rank == 1:
        return NID
    else:
        return NID**(Rank - 1) * (NID + 1) / 2


def MultiID(NID, Ind):
    """Returns a unique ID for multiple atom IDs."""
    if Ind[0] > Ind[-1]:
        x = Ind[0] * (Ind[0] + 1) / 2 + Ind[-1]
    else:
        x = Ind[-1] * (Ind[-1] + 1) / 2 + Ind[0]
    if len(Ind) > 2:
        y = np.product([i * NID**n for (n,i) in enumerate(Ind[1:-1])])
    else:
        y = 0
    return NPairID(NID) * y + x
    

def ReduceBondOrdMatrix(BondOrdMatrix, BondOrdLimit = BondOrdLimit):
    """Reduces a bond order matrix to compact form.  
Returns BondOrdData, BondOrdStart, BondOrdShift.
BondOrdData[BondOrdStart[SIDi] + SIDj - BondOrdShift[SIDi]] gives the bond order between i and j."""
    BondOrdData = []
    BondOrdStart = []
    BondOrdShift = []
    for (SIDi, Row) in enumerate(BondOrdMatrix):
        #filter out too high bond orders
        Row = Row.copy()
        Row[np.abs(Row) > BondOrdLimit] = BondOrdLimit + 1
        #find where normal bond orders start and stop
        UseInd = np.where(Row <= BondOrdLimit)[0]
        SIDStart, SIDStop = UseInd[0], UseInd[-1]
        #update the bond order data
        BondOrdStart.append(len(BondOrdData))
        BondOrdShift.append(SIDStart)
        BondOrdData.extend(Row[SIDStart:SIDStop+1])
    BondOrdStart.append(len(BondOrdData))
    #convert to arrays
    BondOrdData = np.array(BondOrdData, dtype=int)
    BondOrdStart = np.array(BondOrdStart, dtype=int)
    BondOrdShift = np.array(BondOrdShift, dtype=int)
    return BondOrdData, BondOrdStart, BondOrdShift
    



class AtomType(object):
    LockedVars = ["Mass", "AID"]    
    
    def __init__(self, Name = "", Element = None, Mass = 1.,
                 Charge = 0.0, FCharge = 0.0,
                 Color = (0,1,0), Radius = 1.0, Opacity=1.):
        """Creates a new atom type.
Arguments:
    Name: string name of atom
    Element: string name of element
    Mass: mass, float
    Charge: partial charge, float
    FCharge: formal charge, float
    Color: (R,G,B) tuple for visualized color
    Radius: radius for visualization
"""
        #variable to indicate that atom type is locked; automatically set
        self.Locked = False        
        self.Name = Name
        self.Element = Element
        self.Mass = Mass
        self.Charge = Charge
        self.FCharge = FCharge
        #unique type identifier for this atom
        self.AID = 0
        #visualization variables
        self.Color = Color
        self.Radius = Radius
        self.Opacity = Opacity

    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("Atom type object already locked and loaded.")
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Locked = False        

    def New(self, Parent = None):
        """Returns a new atom instance of self."""
        return Atom(self, Parent = Parent)

    def __repr__(self):
        return self.Name   
        
    def __setattr__(self, name, val):
        if name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify atom type property %s because it is locked." % name)
        else:
            object.__setattr__(self, name, val)


class SiteType(object):
    ALinkVars = ["Name", "Element", "Mass", "Charge", "FCharge", "AID"]
    MLinkVars = ["MID"]
    LockedVars = ["AType", "MType", "Parent", "AInd", "SID"]
    
    def __init__(self, AType, MType):
        """Creates a new atom-within-a-molecule type.
Arguments:
    AType: atom type
    MType: molecule type
"""
        #variable to indicate that site type is locked; automatically set
        self.Locked = False        
        self.AType = AType
        self.MType = MType
        self.Parent = MType
        #index within molecule
        self.AInd = 0
        #unique type identifier for this atom
        self.SID = 0
        #visualization variables
        self.Color = AType.Color
        self.Radius = AType.Radius
        self.Opacity = AType.Opacity
        
    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("Site type object already locked and loaded.")
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        if not self.AType is None: 
            self.AType.Cleanup()
        self.Locked = False
        self.AType = None
        self.MType = None
        self.Parent = None

    def New(self, Parent = None):
        """Returns a new atom instance of self."""
        return Atom(self, Parent = Parent)

    def AtomBonds(self, Rigid = None):
        """Returns a list of bonded site types."""
        return self.MType.AtomBonds(self, Rigid = Rigid)

    def __getattr__(self, name):
        if name in self.ALinkVars:
            return getattr(self.AType, name)
        elif name in self.MLinkVars:
            return getattr(self.MType, name)
        elif name == "Label":
            if self.MType is None:
                raise ValueError("Cannot access Label until site has a parent molecule.")
            else:
                return "%s:%s:%d" % (self.MType.Name, self.AType.Name, self.AInd)
        else:
            raise AttributeError(name)

    def __setattr__(self, name, val):
        if name in self.ALinkVars:
            raise AttributeError("Cannot set %s from this object" % name)
        elif name in self.MLinkVars:
            raise AttributeError("Cannot set %s from this object" % name)
        elif name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify atom site type property %s because it is locked." % name)
        else:
            object.__setattr__(self, name, val)

    def __repr__(self):
        return "%s%d" % (self.AType.Name, self.AInd)
    


class Atom(object):
    ArrayVars = ["Pos", "Vel", "Force"]
    SLinkVars = ["AType", "MType", "AInd", "SID"]
    ALinkVars = SiteType.ALinkVars
    MLinkVars = ["MID"]
    LockedVars = ["SType", "Parent"]
    
    def __init__(self, SType, Parent = None):
        """Creates a new atom based on SType."""
        #variable to indicate that atom is locked; automatically set
        self.Locked = False
        self.SType = SType
        #index in master system arrays
        self.Ind = 0
        #parent
        self.Parent = Parent
        #visualization variables
        self.Color = SType.Color
        self.Radius = SType.Radius
        self.Opacity = SType.Opacity
        #variable to indicate that atom is locked; set automatically
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Locked = False
        self.Parent = None
        self.SType = None        

    def AtomBonds(self, Rigid = None):
        """Returns a list of bonded atoms."""
        return self.Parent.AtomBonds(self, Rigid = Rigid)

    def __repr__(self):
        return repr(self.SType)

    def __getattr__(self, name):
        if name in Atom.ArrayVars:
            v = getattr(self.Parent.Parent, name)[self.Ind]
            try:
                v = getattr(self.Parent.Parent, name)[self.Ind]
            except AttributeError:
                raise AttributeError("Cannot access Atom.%s until linked to a System" % name)
            return v
        elif name in self.SLinkVars:
            return getattr(self.SType, name)
        elif name in self.ALinkVars:
            return getattr(self.SType.AType, name)
        elif name in self.MLinkVars:
            return getattr(self.SType.MType, name)
        else:
            raise AttributeError(name)

    def __setattr__(self, name, val):
        if name in Atom.ArrayVars:
            try:
                v = getattr(self.Parent.Parent, name)
            except AttributeError:
                raise AttributeError("Cannot access Atom.%s until linked to a System" % name)
            v[self.Ind] = val
        elif name in self.SLinkVars or name in self.ALinkVars or name in self.MLinkVars:
            raise AttributeError("Cannot set %s from this object" % name)
        elif name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify atom property %s because it is locked." % name)
        else:
            object.__setattr__(self, name, val)    

 

class Bond(object):
    LockedVars = ["SType1", "SType2"]
    
    def __init__(self, SType1, SType2, RLength = None):
        """Creates a new bond.
Arguments:
    SType1, SType2: SiteType instances
    RLength: length of bond if rigid, float
"""
        #variable to indicate that bond type is locked; automatically set
        self.Locked = False
        self.SType1 = SType1
        self.SType2 = SType2
        if not self.SType1.MType == self.SType2.MType:
            raise TypeError("Site types must be in the same molecule.")
        if self.SType1 == self.SType2:
            raise ValueError("Cannot bond sites to themselves.")
        self.RLength = RLength

    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("Bond object already locked and loaded.")
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Locked = False
        self.SType1 = None
        self.SType2 = None

    def GetAtom(self, SType, Rigid = None):
        """Returns the other bonded atom if SType matches one."""
        if not Rigid is None and Rigid == (self.RLength > 0):
            return None
        if self.SType1 == SType:
            return self.SType2
        elif self.SType2 == SType:
            return self.SType1
        else:
            return None

    def IsBond(self, SType1, SType2, Rigid = None):
        """Returns True if matches bond."""
        if not Rigid is None and Rigid == (self.RLength > 0):
            return False
        return ((self.SType1 == SType1 and self.SType2 == SType2)
            or (self.SType2 == SType1 and self.SType1 == SType2))    

    def __eq__(self, other):
        """Returns True if bond atoms are the same."""
        if isinstance(other, Bond):
            s1, s2 = self.SType1, self.SType2
            o1, o2 = other.SType1, other.SType2
            return (s1==o1 and s2==o2) or (s2==o1 and s1==o2)
        else:
            raise TypeError("Did not find a Bond type.")

    def __ne__(self, other):
        """Returns True if bond atoms are not the same."""
        return not self.__eq__(other)      

    def __repr__(self):
        return "%s-%s" % (self.SType1, self.SType2)
        
    def __getattr__(self, name):
        if name == "Rigid":
            return self.RLength > 0
        else:
            raise AttributeError(name)
        
    def __setattr__(self, name, val):
        if name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify bond property %s because it is locked." % name)
        else:
            object.__setattr__(self, name, val)
        


class BondList(list):
    LockedVars = ["HasRigid", "NRigid"]
    
    def __init__(self):
        """Creates a new bond list."""
        list.__init__(self)
        #variable to indicate if bond list is locked; set automatically
        self.Locked = False
        #whether or not has a rigid bond
        self.HasRigid = False
        #number of rigid bonds
        self.NRigid = 0
        
    def __UpdateRigid(self):
        """Updates variables describing rigid bonds."""
        self.HasRigid = False
        self.NRigid = 0
        for Bond in self:
            if Bond.Rigid:
                self.HasRigid = True
                self.NRigid += 1        
        
    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("Bond list object already locked and loaded.")
        for Bond in self:
            Bond.Load()
        self.__UpdateRigid()
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Locked = False
        for Bond in self:
            Bond.Cleanup()
        del self[:]
       
    def append(self, Bond1):
        """Adds a bond to the list."""
        if not isinstance(Bond1, Bond):
            raise TypeError("Did not find a Bond type.")
        if self.Locked:
            raise ValueError("Cannot add bond because bond list is locked.")
        #check if it's already there
        if Bond1 in self:
            raise ValueError("Bond %s already added to bond list." % str(Bond1))
        list.append(self, Bond1)
        if Bond1.Rigid:
            self.HasRigid = True
            self.NRigid += 1

    def __add__(self, other):
        self.append(other)
        return self

    def extend(self, otherlist):
        for other in otherlist:
            self.append(other)

    def __delitem__(self, key):
        if self.Locked:
            raise ValueError("Cannot delete bond because bond list is locked.")
        if key in self:
            ind = self.index(key)
            list.__delitem__(self, ind)
        elif len(key) == 2:
            b = Bond(*key)
            list.__delitem__(self, self.index(b))
        else:
            list.__delitem__(self, key)
        self.__UpdateRigid()

    def __repr__(self):
        return list.__repr__(self)
        
    def __setattr__(self, name, val):
        if name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify atom type property %s because it is locked." % name)
        else:
            list.__setattr__(self, name, val)



class MolType(list):
    LockedVars = ["MID", "Rigid", "Bonds", "BondOrdData", "BondOrdStart", "BondOrdShift"]
    
    def __init__(self, Name = "", AtomTypes = [],
                 Rigid = False, COMPos = None):
        """Creates a new molecule type.
Arguments:
    Name: string name of molecule
    AtomTypes: list of Atom atom types
"""
        list.__init__(self)
        #variable to indicate if mol type is locked; set automatically
        self.Locked = False
        self.Name = Name
        self.extend(AtomTypes)
        #unique type identifier for this molecule
        self.MID = 0
        self.Rigid = Rigid
        self.COMPos = COMPos
        self.Bonds = BondList()
        self.__Update()
        #bond order data
        self.__ClearBondOrd()

    def __repr__(self):
        return "%s:%s" % (self.Name, list.__repr__(self))
        
    def __ClearBondOrd(self):
        """Clears internal bond order information."""
        self.BondOrdData = None
        self.BondOrdStart = None
        self.BondOrdShift = None  

    def __UpdateBondOrd(self):
        """UPdates internal bond order information."""
        BondOrdMatrix = self.GetBondOrdMatrix(ShowRigid = True)
        self.BondOrdData, self.BondOrdStart, self.BondOrdShift = ReduceBondOrdMatrix(BondOrdMatrix)
            
    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("Mol type object already locked and loaded.")
        self.__Update()
        for SType in self:
            SType.Load()
        self.Bonds.Load()
        self.__UpdateBondOrd()
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        self.Bonds.Cleanup()
        for SType in self:
            SType.Cleanup()
        self.Locked = False
        del self[:]        

    def __Update(self):
        """Updates internal numbering."""
        for (i, SType) in enumerate(self):
            SType.AInd = i

    def New(self):
        """Returns a new instance of self."""
        return Mol(self)

    def append(self, AType):
        """Adds an atom type to the list."""
        if not isinstance(AType, AtomType):
            raise TypeError("Did not find a AtomType type.")
        if self.Locked:
            raise ValueError("Cannot add atom type to mol type because it is locked.")
        list.append(self, SiteType(AType, self))
        self[-1].AInd = len(self) - 1
        self.__ClearBondOrd()

    def __add__(self, other):
        self.append(other)
        return self

    def extend(self, otherlist):
        for other in otherlist:
            self.append(other)

    def __delitem__(self, key):
        if self.Locked:
            raise ValueError("Cannot delete atom type from mol type because it is locked.")
        if key in self:
            ind = self.index(key)
            list.__delitem__(self, ind)
        else:
            list.__delitem__(self, key)
        self.__Update()
        self.__ClearBondOrd()
        
    def __getattr__(self, name):
        if name == "HasRigid":
            return self.Rigid or self.Bonds.HasRigid
        else:
            raise AttributeError(name)
        
    def __setattr__(self, name, val):
        if name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify mol type property %s because it is locked." % name)
        else:
            list.__setattr__(self, name, val)
            
    def Bond(self, SType1, SType2, RLength = None):
        """Adds a bond between SType1 and SType2."""
        if type(SType1) is int: SType1 = self[SType1]
        if type(SType2) is int: SType2 = self[SType2]
        self.Bonds.append(Bond(SType1, SType2, RLength))
        self.__ClearBondOrd()
            
    def GetBondOrdMatrix(self, ShowRigid = False):
        """Returns a matrix with the order between each pair of sites.
ShowRigid: If True, will use (-2) for rigid bonds and (+2) for nonrigid."""
        #make a matrix of bond orders between atoms
        n = len(self)
        m = np.identity(n, dtype=int)
        #fill with bonds
        for Bond in self.Bonds:
            ID1 = Bond.SType1.AInd
            ID2 = Bond.SType2.AInd
            if ID1 < 0 or ID1 >= n or ID2 < 0 or ID2 >= n:
                continue
            if ShowRigid and Bond.Rigid:
                m[ID1, ID2] = -2
                m[ID2, ID1] = -2
            else:
                m[ID1, ID2] = 2
                m[ID2, ID1] = 2
        #keep filling in new bond orders until we can't find any more
        CurOrd = 3
        Found = True
        while Found:
            Found = False
            #find elements in matrix that still have no bond order
            Ind1, Ind2 = (m==0).nonzero()
            #filter out only one half of matrix since symmetric
            Mask = Ind1 < Ind2
            Ind1, Ind2 = Ind1[Mask], Ind2[Mask]
            #loop through these elements
            Ind12 = np.concatenate((Ind1[:,np.newaxis], Ind2[:,np.newaxis]), axis=1)
            #make an absolute value matrix
            mabs = np.abs(m)
            for (ID1, ID2) in Ind12:
                #see if there are any bonds to a common third atom
                BondedMask = np.logical_and(mabs[ID1,:] > 1, mabs[ID2,:] > 1)
                #see if the sum of the bond orders equals the current
                OrderMask = np.logical_and(BondedMask, mabs[ID1,:] + mabs[ID2,:] - 1 == CurOrd)
                #check if there were any links at this bond order
                if np.any(OrderMask):
                    m[ID1, ID2] = CurOrd
                    m[ID2, ID1] = CurOrd
                    Found = True
            CurOrd += 1
        return m
        
    def BondOrder(self, SID1, SID2, ShowRigid = False):
        """Returns the order of a bond.  If ShowRigid==True, then will return 
-2 for rigid bonds of order 2; otherwise, will return positive numbers."""
        if type(SID1) is SiteType: SID1 = SID1.AInd
        if type(SID2) is SiteType: SID2 = SID2.AInd
        if self.BondOrdData is None:
            self.__UpdateBondOrd()
        StartInd = self.BondOrdStart[SID1]
        StopInd = self.BondOrdStart[SID1+1]
        ShiftInd = self.BondOrdShift[SID1]
        DataInd = SID2 - ShiftInd + StartInd
        if DataInd < StartInd or DataInd >= StopInd:
            return BondOrdLimit + 1
        else:
            bo = self.BondOrdData[DataInd]
            if ShowRigid:
                return bo
            else:
                return abs(bo)
        
    def Bonded(self, SID1, SID2, Rigid = None):
        """Returns true if SID1 and SID2 are bonded."""
        bo = self.BondOrder(SID1, SID2, ShowRigid = True)
        if Rigid is None:
            return abs(bo) == 2
        elif Rigid == True:
            return bo == -2
        else:
            return bo == 2

    def AtomBonds(self, SID, Rigid = None, BondOrd = 2):
        """Returns all bonds from given atom."""
        if self.BondOrdData is None:
            self.__UpdateBondOrd()   
        if type(SID) is SiteType: SID = SID.AInd
        if BondOrd > BondOrdLimit:
            raise ValueError("Specified bond order of %d is greater than tabulated limit %d." % (BondOrd, BondOrdLimit))
        StartInd = self.BondOrdStart[SID]
        StopInd = self.BondOrdStart[SID+1]
        ShiftInd = self.BondOrdShift[SID]
        ThisBondOrd = self.BondOrdData[StartInd:StopInd]
        if BondOrd > 2:
            Ind = np.where(ThisBondOrd == BondOrd)[0]
        elif Rigid is None:
            Ind = np.where(np.abs(ThisBondOrd) == 2)[0]
        elif Rigid == True:
            Ind = np.where(ThisBondOrd == -2)[0]
        else:
            Ind = np.where(ThisBondOrd == 2)[0]
        return [self[x + ShiftInd] for x in Ind]  

    def BondMap(self, Rigid = None, BondOrd = 2):
        """Returns a mapping of bonds from sites i to j."""
        return [self.AtomBonds(a.AInd, Rigid, BondOrd) for a in self]

    def WithinBondOrd(self, SID1, SID2, MinBondOrd = None, MaxBondOrd = None):
        """Returns true if two site types are bonded within MinBondOrd and MaxBondOrd."""
        #check to make sure not larger than tabulated
        if not MinBondOrd is None and MinBondOrd > BondOrdLimit + 1:
            raise ValueError("MinBondOrd of %d is larger than tabulated limit of %d in BondOrdLimit" % (MinBondOrd, BondOrdLimit))
        if not MaxBondOrd is None and MaxBondOrd > BondOrdLimit:
            raise ValueError("MaxBondOrd of %d is larger than tabulated limit of %d in BondOrdLimit" % (MaxBondOrd, BondOrdLimit))
        bo = self.BondOrder(SID1, SID2, ShowRigid = False)
        return (MinBondOrd is None or bo >= MinBondOrd) and (MaxBondOrd is None or bo <= MaxBondOrd) 

    
        
class Mol(list):
    MLinkVars = ["Name", "MID", "Rigid", "Bonds"]
    LockedVars = ["MType", "Parent"]
    
    def __init__(self, MType):
        """Creates a new molecule based on MType."""
        list.__init__(self)
        #variable to indicate if mol is locked; set automatically
        self.Locked = False
        self.MType = MType
        #index of molecule
        self.MInd = 0
        #make new atoms from site types
        self.extend([a.New(Parent = self) for a in MType])
        self.Pos = None
        self.Orient = None
        self.Parent = None
        #variable to indicate if mol is locked; set automatically
        self.Locked = True

    def __repr__(self):
        return "%s%d:%s" % (self.Name, self.MInd, list.__repr__(self.MType))

    def __delitem__(self, key):
        if self.Locked:
            raise ValueError("Cannot delete atom from molecule because it is locked.")
        if key in self:
            ind = self.index(key)
            list.__delitem__(self, ind)
        else:
            list.__delitem__(self, key)

    def __getattr__(self, name):
        if name in self.MLinkVars:
            return getattr(self.MType, name)
        elif name == "HasRigid":
            return self.MType.HasRigid
        else:
            raise AttributeError(name)

    def __setattr__(self, name, val):
        if name in self.MLinkVars:
            raise AttributeError("Cannot set %s from this object" % name)
        elif name in self.LockedVars and getattr(self, "Locked", False) and not getattr(self, name, None) is None:
            raise ValueError("Cannot modify mol property %s because it is locked." % name)
        else:
            list.__setattr__(self, name, val)

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        for Atom in self:
            Atom.Cleanup()
        self.Locked = False
        del self[:]
        self.MType = None
        self.Parent = None

    def Bonded(self, Atom1, Atom2, Rigid = None):
        """Returns true if Atom1 and Atom2 are bonded."""
        if type(Atom1) is int: Atom1 = self[Atom1]
        if type(Atom2) is int: Atom2 = self[Atom2]
        if not Atom1.Parent is Atom2.Parent:
            return False
        return self.MType.Bonds.Bonded(Atom1.SType.AInd, Atom2.SType.AInd, Rigid)

    def AtomBonds(self, Atom, Rigid = None, BondOrd = 2):
        """Returns all bonds from given atom."""
        if type(Atom) is int: Atom = self[Atom]
        STypes = self.MType.AtomBonds(Atom.SType.AInd, Rigid, BondOrd)
        return [self[SType.AInd] for SType in STypes]

    def BondMap(self, Rigid = None, BondOrd = 2):
        """Returns a mapping of bonds from atom sites i to j."""
        return [self.AtomBonds(a, Rigid, BondOrd) for a in self]

    def WithinBondOrd(self, Atom1, Atom2, MinBondOrd = None, MaxBondOrd = None):
        """Returns true if two site types are bonded within MaxBondOrd."""
        if type(Atom1) is int: Atom1 = self[Atom1]
        if type(Atom2) is int: Atom2 = self[Atom2]
        if not Atom1.Parent is Atom2.Parent:
            return False
        return self.MType.WithinBondOrd(Atom1.SType.AInd, Atom2.SType.AInd,
                                        MinBondOrd, MaxBondOrd)       



class World(list):
    LockedVars = ["Dim", "Units", "BondOrdData", "BondOrdStart", "BondOrdShift"]
    ModuleConsts = ["Dim", "BondOrdData", "BondOrdStart", "BondOrdShift",
                    "BondOrdLimit", "NAID", "NSID", "NMID"]
    
    def __init__(self, MolTypes, Dim, Units):
        """Creates a new list of molecular components."""
        list.__init__(self)
        #variable to indicate if world is locked; set automatically
        self.Locked = False
        self.Dim = Dim
        if not isinstance(Units, units.UnitsClass):
            raise TypeError("Units is not a UnitsClass instance.")
        self.Units = Units
        self.extend(MolTypes) 
        #variable to indicate if locked; set automatically
        self.Locked = False
        #bond order data
        self.BondOrdData = None
        self.BondOrdStart = None
        self.BondOrdShift = None
        self.BondOrdLimit = BondOrdLimit

    def __Update(self):
        """Updates numbering and variables."""
        self.SiteTypes = []
        self.AtomTypes = []
        #set all atom ids to zero
        for MType in self:
            for SType in MType:
                SType.AType.AID = None
        #count and renumber
        self.NAID = 0
        self.NSID = 0
        self.NMID = 0
        for MType in self:
            MType.MID = self.NMID
            self.NMID += 1
            for SType in MType:
                self.SiteTypes.append(SType)
                SType.SID = self.NSID
                self.NSID += 1
                AType = SType.AType
                if AType.AID is None:
                    self.AtomTypes.append(AType)
                    AType.AID = self.NAID
                    self.NAID += 1
                    
    def __UpdateBondOrd(self):
        """Updates internal bond order information."""
        BondOrdMatrix = self.GetBondOrdMatrix(ShowRigid = True)
        self.BondOrdData, self.BondOrdStart, self.BondOrdShift = ReduceBondOrdMatrix(BondOrdMatrix)

    def Load(self):
        """Prepares data for fast access and locks object."""
        if self.Locked:
            raise ValueError("World object already locked and loaded.")
        self.__Update()
        for MType in self:
            MType.Load()
        #make bond order data
        self.__UpdateBondOrd()
        self.Locked = True

    def Cleanup(self):
        """Removes links and prepares for garbage collection."""
        for MType in self:
            MType.Cleanup()
        for SType in self.SiteTypes:
            SType.Cleanup()
        for AType in self.AtomTypes:
            AType.Cleanup()
        self.Locked = False
        del self[:]
        del self.AtomTypes[:]
        del self.SiteTypes[:]              
        
    def append(self, MType):
        """Adds a molecule to the list."""
        if not isinstance(MType, MolType):
            raise TypeError("Did not find a MolType type.")
        if self.Locked:
            raise ValueError("Cannot add mol type to world because it is locked.")
        list.append(self, MType)
        self.__Update()

    def __add__(self, other):
        self.append(other)
        return self

    def extend(self, otherlist):
        for other in otherlist:
            self.append(other)

    def __delitem__(self, key):
        if self.Locked:
            raise ValueError("Cannot delete mol type from world because it is locked.")
        if key in self:
            ind = self.index(key)
            list.__delitem__(self, ind)
        else:
            list.__delitem__(self, key)
        self.__Update()
        
    def __setattr__(self, name, val):
        if hasattr(self, name):
            if name == "BondOrdLimit" or (name in self.LockedVars and getattr(self, "Locked", False)):
                raise ValueError("Cannot modify world property %s because it is locked." % name)
            else:
                list.__setattr__(self, name, val)
        else:
            list.__setattr__(self, name, val)

    def AtomBonds(self, Atom, Rigid = None, BondOrd = 2):
        """Returns a list of bonds from atom."""
        return Atom.Parent.AtomBonds(Atom, Rigid, BondOrd)

    def Bonds(self, Rigid = None, BondOrd = 2):
        """Returns a list of sites that are bonded."""
        ret = []
        for MType in self:
            ret.extend(MType.BondMap(Rigid, BondOrd))
        return ret

    def BondMap(self, Rigid = None, BondOrd = 2):
        """Returns a list of lists of bonded Sites."""
        ret = []
        for MType in self:
            ret.extend(MType.BondMap(Rigid, BondOrd))
        return ret
        
    def GetBondOrdMatrix(self, ShowRigid = False):
        """Returns a matrix with the order between each pair of site types.
ShowRigid: If True, will use (-2) for rigid bonds and (+2) for nonrigid."""
        n = len(self.SiteTypes)
        m = np.zeros((n,n), dtype=int)
        for MType in self:
            StartSID = MType[0].SID
            StopSID = StartSID + len(MType)
            mMol = MType.GetBondOrdMatrix(ShowRigid = True)
            m[StartSID:StopSID,StartSID:StopSID] = mMol
        if not ShowRigid:
            m = np.abs(m)
        return m
    
    def HasMolType(self, M):
        """Returns True if M is a molecule or molecule type in this world."""
        if isinstance(M, Mol):
            return M.MType in self
        elif isinstance(M, MolType):
            return M in self
        else:
            raise TypeError("I do not recognize M as Mol or MolType.")
            
    def GetPairTypes(self):
        """Returns NTypes, TypeIDs, TypeLabels.
NTypes is the number of unique atom pair types.
TypeIDs are tuples of the atom IDs.  
TypeLabels gives the atom name combinations for these pairs."""
        NTypes = 0
        TypeIDs = []
        TypeLabels = []
        nAID = len(self.AtomTypes)
        for i in range(nAID):
            ai = self.AtomTypes[i]
            for j in range(0, i+1):
                aj = self.AtomTypes[j]
                NTypes += 1
                TypeIDs += [(ai.AID, aj.AID)]
                TypeLabels += ["%s-%s" % (ai.Name, aj.Name)]
        return NTypes, TypeIDs, TypeLabels
             



        
