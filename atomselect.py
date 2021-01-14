#!usr/bin/env python

### Class definitions for atoms in SIM suite.
### coded by MSS

import copy

import numpy as np

import sim.chem as chem



def AtomIter(Atoms, N = 1, Bonded = None, Intra = None,
             MinBondOrd = None, MaxBondOrd = None,
             Filters = None, Prev = []):
    """Returns an iterator over either site or atom groups.
Atoms: list of atoms from which to make groups
N: number of atoms in each group
Bonded: whether successive atoms are bonded or not
Intra: whether all atoms are in the same molecule or not
MinBondOrd: minimum bond order to consider (inclusive)
MaxBondOrd: maximum bond order to consider (inclusive)
Filters: optional filters to apply
Prev: reserved variable for recursive search"""
    #make lists of current atoms
    ThisAtoms = Atoms
    #filter out already used atoms if these are not just types because
    #types could interact between the same molecule
    if Bonded == True or isinstance(Atoms[0], chem.Atom):
        ThisAtoms = [a for a in ThisAtoms if not a in Prev]
    #check user-supplied filter functions
    if not Filters is None:
        ThisAtoms = [a for a in ThisAtoms if Filters[0](a)]
    #apply filters
    if len(Prev):
        #See if this is an actual atom rather than a site type;
        #Don't do any bond order checking, etc, on site types because
        #these could be interactions between different molecules
        IsAtom = isinstance(Prev[0], chem.Atom)
        #check for bonded spec
        if Bonded == True:
            ThisAtoms = [a for a in ThisAtoms if a in Prev[-1].AtomBonds()]
        elif IsAtom and Bonded == False:
            ThisAtoms = [a for a in ThisAtoms if not a in Prev[-1].AtomBonds()]
        #check for bond order spec
        if IsAtom and not (MinBondOrd is None and MaxBondOrd is None):
            Parent = Prev[-1].Parent
            ThisAtoms = [a for a in ThisAtoms 
                         if Parent.WithinBondOrd(Prev[-1], a, MinBondOrd, MaxBondOrd)]
        #check for intra spec
        if Intra == True:
            ThisAtoms = [a for a in ThisAtoms if a in Prev[-1].Parent]
        elif IsAtom and Intra == False:
            ThisAtoms = [a for a in ThisAtoms if not a in Prev[-1].Parent]    
    #now select the atoms
    if N == 1:
        for a in ThisAtoms:
            yield [a]
        return
    elif N > 1:
        for a1 in ThisAtoms:
            for a2list in AtomIter(Atoms, N-1, Bonded, Intra,
                                   MinBondOrd, MaxBondOrd, 
                                   Filters = Filters[1:], Prev = Prev + [a1]):
                yield [a1] + a2list
        return
    else:
        raise IndexError("Need to specify one or more atoms/sites.")
   


class Filter:
    def __init__(self, Types = None, MolTypes = None,
                 Name = None, Element = None,
                 HasCharge = None, HasFCharge = None, Negate = False):
        """Creates a filter for atoms meeting specified criteria."""
        self.Types = Types
        if not self.Types is None and not hasattr(self.Types, "__getitem__"):
            self.Types = [self.Types]
        self.MolTypes = MolTypes 
        if isinstance(self.MolTypes, chem.MolType):
            self.MolTypes = [self.MolTypes]
        self.Name = Name
        self.Element = Element
        self.HasCharge = HasCharge
        self.HasFCharge = HasFCharge
        self.Negate = Negate

    def copy(self):
        return copy.copy(self)

    def HasAtomType(self):
        if self.Types is None: return False
        return any([isinstance(Type, chem.AtomType) for Type in self.Types])

    def HasSiteType(self):
        if self.Types is None: return False
        return any([isinstance(Type, chem.SiteType) for Type in self.Types])

    def HasMolType(self):
        return not self.MolTypes is None
      
    def Query(self, Atom, MType = None):
        """Returns True if applies to a given atom or site type."""
        if isinstance(Atom, chem.Atom):
            AType = Atom.AType
            SType = Atom.SType
            MType = Atom.MType
        elif isinstance(Atom, chem.SiteType):
            AType = Atom.AType
            SType = Atom
            MType = Atom.MType
        elif isinstance(Atom, chem.AtomType):
            AType = Atom
            SType = None
        else:
            raise TypeError("Type of Atom not recognized.")
        ret = True
        if not self.Types is None:
            ret = ret and (AType in self.Types or SType in self.Types)
        if not self.MolTypes is None and not MType is None:
            ret = ret and (MType in self.MolTypes)
        if not self.Name is None:
            ret = ret and (self.Name.strip() == AType.Name.strip())
        if not self.Element is None:
            ret = ret and (self.Element.strip() == AType.Element.strip())
        if not self.HasCharge is None:
            ret = ret and (self.HasCharge == (AType.Charge != 0))
        if not self.HasFCharge is None:
            ret = ret and (self.HasFCharge == (AType.FCharge != 0))
        return ret == (not self.Negate)

    def __call__(self, arg):
        return self.Query(arg)

    def Select(self, Atoms):
        """Returns all matching atoms or sites in arg."""
        for Atom in Atoms:
            if self.Query(Atom):
                yield Atom

    def __add__(self, other):
        if isinstance(other, Filter):
            return PolyFilter(Filters = [self, other])
        elif isinstance(other, PolyFilter):
            newother = other.copy()
            newother.append(self)
            return newother
        elif isinstance(other, FilterList):
            if other.N != self.N:
                raise TypeError("Added filter must be the same length as existing.")
            newother = other.copy()
            newother.append(self)
            return newother
        else:
            raise TypeError("Expected Filter or PolyFilter object.")

    def __repr__(self):
        l = []
        if not self.Types is None:
            l.append("Types:%s" % repr(self.Types))
        if not self.MolTypes is None:
            l.append("MolTypes:%s" % repr(self.MolTypes))
        if not self.Name is None:
            l.append("Name:%s" % self.Name)
        if not self.Element is None:
            l.append("Element:%s" % repr(self.Element))
        if not self.HasCharge is None:
            if self.HasCharge:
                l.append("Charge")
            else:
                l.append("NoCharge")
        if not self.HasFCharge is None:
            if self.HasFCharge:
                l.append("FCharge")
            else:
                l.append("NoFCharge")
        if self.Negate:
            l = ["NOT"] + l
        return "{" + " ".join(l) + "}"

    def __getattr__(self, name):
        if name == "N":
            return 1
        else:
            raise AttributeError(name)

    def __neg__(self):
        ret = self.copy()
        ret.Negate = (not ret.Negate)
        return ret


class PolyFilter(list):
    def __init__(self, Filters = [], N = None, 
                 Intra = None, Bonded = None,
                 MinBondOrd = None, MaxBondOrd = None, Ordered = False,
                 AIDPairs = None):
        """Creates a multiple-atom filter."""
        list.__init__(self)
        NewFilters = []
        for x in Filters:
            if isinstance(x, Filter):
                NewFilters.append(x)
            else:
                NewFilters.append(Filter(x))
        Filters = NewFilters
        if N is None:
            self.extend(Filters)
        else:
            self.extend([Filter() for i in range(N)])
        self.Intra = Intra
        self.Bonded = Bonded
        self.MinBondOrd = MinBondOrd
        self.MaxBondOrd = MaxBondOrd
        #check bond order specs
        if not self.MinBondOrd is None:
            #check to make sure not larger than tabulated
            if self.MinBondOrd > chem.BondOrdLimit + 1:
                raise ValueError("MinBondOrd of %d is larger than tabulated limit of %d in chem.BondOrdLimit" % (self.MinBondOrd, chem.BondOrdLimit))
        if not self.MaxBondOrd is None:
            self.Intra = True
            #check to make sure not larger than tabulated
            if self.MaxBondOrd > chem.BondOrdLimit:
                raise ValueError("MaxBondOrd of %d is larger than tabulated limit of %d in chem.BondOrdLimit" % (self.MaxBondOrd, chem.BondOrdLimit))               
        #check if we just want bonded
        if self.MaxBondOrd == 2:
            self.Bonded = True
            self.MaxBondOrd = None
            self.MinBondOrd = None
        #ordered
        self.Ordered = Ordered
        #check for a matrix of AType
        self.AIDPairs = None
        if not AIDPairs is None:
            #check for a pair
            if not len(self) == 2:
                raise ValueError("AIDPairs only valid for pair interactions.")
            if not type(AIDPairs) is np.ndarray:
                raise ValueError("AIDPairs must be a numpy array.")
            if not np.all(AIDPairs == AIDPairs.T):
                raise ValueError("AIDPairs is not symmetric.")
            self.AIDPairs = np.array(AIDPairs, dtype=bool)

    def copy(self):
        return copy.copy(self)
        
    def Reverse(self):
        c = self.copy()
        c.reverse()
        if not c.AIDPairs is None:
            c.AIDPairs = c.AIDPairs[::-1,::-1]
        return c

    def HasSiteType(self):
        return any([Filter.HasSiteType() for Filter in self])

    def HasAtomType(self):
        return any([Filter.HasAtomType() for Filter in self])

    def HasMolType(self):
        return any([Filter.HasMolType() for Filter in self])
        
    def HasIntra(self):
        return not (self.Intra is None)

    def Query(self, Atoms):
        """Returns true if applies to atoms or site types."""
        if len(Atoms) != self.N:
            raise IndexError("Expected %d atoms in Atoms" % self.N)

        AllAtoms = all([isinstance(Atom, chem.Atom) for Atom in Atoms])
        AllSites = all([isinstance(Atom, chem.SiteType) for Atom in Atoms])
        if not (AllAtoms or AllSites):
            raise TypeError("Do not recognize types in Atoms") 
        
        Mol = Atoms[0].Parent
        SameMol = all([Atom.Parent is Mol for Atom in Atoms])

        ret = all([Filter(Atom) for (Filter, Atom) in zip(self, Atoms)])
        if not self.Intra is None:
            ret = ret and (self.Intra == SameMol)
        if not self.Bonded is None:
            if SameMol:
                for (Atom1, Atom2) in zip(Atoms[:-1], Atoms[1:]):
                    ret = ret and Mol.Bonded(Atom1, Atom2) == self.Bonded
            else:
                ret = ret and (self.Bonded == False)
        if not self.MaxBondOrd is None or not self.MinBondOrd is None:
            if SameMol:
                ret = ret and Mol.WithinBondOrd(Atoms[0], Atoms[1],
                                                self.MinBondOrd,
                                                self.MaxBondOrd)
            else:
                ret = False  
      
        if not self.AIDPairs is None:
            print Atoms[0].AID, Atoms[1].AID, self.AIDPairs[Atoms[0].AID, Atoms[1].AID] 
            ret = ret and self.AIDPairs[Atoms[0].AID, Atoms[1].AID]       
            
        return ret


    def __call__(self, arg):
        return self.Query(arg)


    def Select(self, Atoms):
        """Returns all matching atoms in Atoms."""
        #now look at all groups of these atoms
        for Group in AtomIter(Atoms, self.N, self.Bonded, self.Intra,
                              self.MinBondOrd, self.MaxBondOrd, Filters = self):
            if not self.AIDPairs is None:
                if not self.AIDPairs[Group[0].AID, Group[1].AID]: continue                                      
            if self.Ordered:
                yield Group
            else:
                #find the reverse query
                QueryReverse = all([Filter(a) for (Filter, a) in zip(self, Group[::-1])])
                #if this is true only yield if QueryReverse has lower site id first
                if not QueryReverse:
                    yield Group
                elif Group[0].SID <= Group[-1].SID or self.N == 1:
                    yield Group      

    def __add__(self, other):
        if isinstance(other, Filter):
            newself = self.copy()
            newself.append(other)
            return newself
        elif isinstance(other, PolyFilter):
            return FilterList([self, other])
        elif isinstance(other, FilterList):
            if other.N:
                if other.N != self.N:
                    raise TypeError("Added filter must be the same length as existing.")
            newother = other.copy()
            newother.append(self)
            return newother
        else:
            raise TypeError("Expected another Filter, PolyFilter, or FilterList object.")      

    def __getattr__(self, name):
        if name == "N":
            return len(self)
        else:
            raise AttributeError(name)


class FilterList(list):
    def __init__(self, arg = []):
        """Creates a filter list, applied using Boolean OR."""
        list.__init__(self)
        self.extend(arg)

    def copy(self):
        return copy.copy(self)

    def HasSiteType(self):
        return any([Filter.HasSiteType() for Filter in self])

    def HasAtomType(self):
        return any([Filter.HasAtomType() for Filter in self])

    def HasMolType(self):
        return any([Filter.HasMolType() for Filter in self]) 
    
    def HasIntra(self):
        return any([Filter.HasIntra()] for Filter in self)

    def Query(self, Atoms):
        """Returns true if applies to atoms."""
        return any([Filter.Query(Atoms) for Filter in self])

    def __call__(self, arg):
        return self.Query(arg)

    def Select(self, Atoms):
        """Returns all matching atoms in Atoms."""
        if not len(self): return
        if not all([self[0].N == x.N for x in self]):
            raise TypeError("All filters must be the same length.")
        Found = []
        for Filter in self:
            for Group in Filter.Select(Atoms):
                if not Group in Found:
                    Found.append(Group)
                    yield Group
        del Found

    def __add__(self, other):
        if isinstance(other, PolyFilter) or isinstance(other, Filter):
            if self.N and other.N != self.N:
                raise TypeError("Added filter must be the same length as existing.")
            newself = self.copy()
            newself.append(other)
            return newself
        else:
            raise TypeError("Expected a PolyFilter object.")

    def __getattr__(self, name):
        if name == "N":
            if len(self):
                return self[0].N
            else:
                return None
        elif name == "Bonded":
            l = [Filter.Bonded == True for Filter in self]
            return all(l)
        elif name == "Intra":
            l = [Filter.Intra == True for Filter in self]
            return all(l)  
        else:
            raise AttributeError(name)

 
        

All = Filter()
Pairs = PolyFilter(Filters = [All, All])
Triples = PolyFilter(Filters = [All, All, All])
Quartets = PolyFilter(Filters = [All, All, All, All])
BondPairs = PolyFilter(Filters = [All, All], Bonded = True)
NonbondPairs = PolyFilter(Filters = [All, All], Bonded = False)
InterPairs = PolyFilter(Filters = [All, All], Intra = False)
IntraPairs = PolyFilter(Filters = [All, All], Intra = True)

BondTriples = PolyFilter(Filters = [All, All, All], Bonded = True)
BondQuartets = PolyFilter(Filters = [All, All, All, All], Bonded = True)

#all nonbond pairs except 1-2 and 1-3
NonbondPairs3 = PolyFilter(Filters = [All, All], MinBondOrd = 4)
NonbondPairs13 = NonbondPairs3
#all nonbond pairs except 1-2, 1-3, and 1-4
NonbondPairs4 = PolyFilter(Filters = [All, All], MinBondOrd = 5)
NonbondPairs14 = NonbondPairs4
