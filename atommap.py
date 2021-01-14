#/usr/bin/env python


### Class definitions for mapping atoms in SIM suite.
### coded by MSS

import numpy as np
import types

import sim.fortran as fortran
import sim.chem as chem


KindCOM = 0
KindCentroid = 1


#======== FORTRAN CODE ========
Lib = fortran.Module("maplib")
Lib.Code = """
subroutine Map(AtomMap1, Weight, StartInd1, NAtom2, Pos1, BoxL, Pos2, &
    & NAtomMap1, NStartInd1, Dim, NAtom1)
    implicit none
    integer, intent(in) :: NAtomMap1, NStartInd1, Dim, NAtom1, NAtom2
    integer, dimension(0:NAtomMap1-1), intent(in) :: AtomMap1
    real(8), dimension(0:NAtomMap1-1), intent(in) :: Weight
    integer, dimension(0:NStartInd1-1), intent(in) :: StartInd1
    real(8), dimension(0:NAtom1-1, 0:Dim-1), intent(in) :: Pos1
    real(8), dimension(0:Dim-1), intent(in) :: BoxL
    real(8), dimension(0:NAtom2-1, 0:Dim-1), intent(out) :: Pos2
    real(8), dimension(Dim) :: iBoxL, Pos1reimaged, Pos1ref
    integer :: i, j, ind1, ind2
    Pos2 = 0.d0
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL > 0.d0)
    do i = 0, NStartInd1 - 2
        ind1 = StartInd1(i)
        ind2 = StartInd1(i+1) - 1
        if (ind1 < 0 .or. ind1 > NAtomMap1 - 1 .or. &
          & ind2 < 0 .or. ind2 > NAtomMap1 - 1) cycle
        do j = ind1, ind2
            Pos1reimaged = Pos1(AtomMap1(j), :)
            if (j == ind1) then
                Pos1ref = Pos1reimaged
            else
                Pos1reimaged = Pos1reimaged - Pos1ref 
                Pos1reimaged = Pos1reimaged - BoxL * anint(Pos1reimaged * iBoxL)
                Pos1reimaged = Pos1reimaged + Pos1ref
            endif
            Pos2(i,:) = Pos2(i,:) + Weight(j) * Pos1reimaged
        enddo
        Pos2(i, :) = Pos2(i,:) - BoxL * anint(Pos2(i,:) * iBoxL)
    enddo
end subroutine
"""


#======== COMPILING ========
Lib.KeepSource = True
Lib.Load()


#======== PYTHON FUNCTIONS ========

class AtomMap(object):
    
    def __init__(self, Atoms1, Atom2, Mass1 = None, Kind = KindCOM, Atom2Name = None):
        """Returns a mapping object for center of masses of Atoms1 to Atom2.
Atoms1: list of Atoms or list of atom indices
Atom2: Atom or atom index
Mass1: list of masses of Atoms1
Kind: whether or not to do center of mass or centroid
Atom2Name: name of atom 2
"""
        if type(Atoms1) is list:
            if all([isinstance(a, chem.Atom) for a in Atoms1]):
                Sys1 = Atoms1[0].Parent.Parent
                if not all([a.Parent.Parent is Sys1 for a in Atoms1]):
                    raise TypeError("All atoms in Atoms1 must be from the same system.")
                self.Atoms1 = np.array([a.Ind for a in Atoms1], int)
                Mass1 = np.array([a.Mass for a in self.Atoms1], float)
            elif all([isinstance(a, int) or isinstance(a, np.integer) for a in Atoms1]):
                self.Atoms1 = np.array(Atoms1, int)
            else:
                raise TypeError("Invalid type for Atoms1.")   
        elif isinstance(Atoms1, int) or isinstance(Atoms1, np.integer):
            self.Atoms1 = np.array([Atoms1], int)
        elif isinstance(Atoms1, chem.Atom):
            self.Atoms1 = np.array([Atoms1.Ind], int)
            Mass1 = np.array([Atoms1.Mass], float)
        elif type(Atoms1) is types.FunctionType:
            #used for generators
            try:
                Atoms1 = [a for a in Atoms1]
            except TypeError:
                raise TypeError("Invalid type for Atoms1.")
            if all([isinstance(a, chem.Atom) for a in Atoms1]):
                if not all([a.Parent.Parent is Sys1 for a in Atoms1]):
                    raise TypeError("All atoms in Atoms1 must be from the same system.")
                self.Atoms1 = np.array([a.Ind for a in Atoms1], int)
                Mass1 = np.array([a.Mass for a in self.Atoms1], float)
            elif all([type(a) is int for a in self.Atoms1]):
                self.Atoms1 = np.array(Atoms1, int)
            else:
                raise TypeError("Invalid type for Atoms1.") 
        else:
            raise TypeError("Invalid type for Atoms1.")

        if isinstance(Atom2, int) or isinstance(Atom2, np.integer):
            self.Atom2 = Atom2
            if Atom2Name is None:
                self.Atom2Name = 'A'
            else:
                self.Atom2Name = Atom2Name
        elif isinstance(Atom2, chem.Atom):
            self.Atom2 = Atom2.Ind
            if Atom2Name is None:
                self.Atom2Name = Atom2.Name
            else:
                self.Atom2Name = Atom2Name
        else:
            raise TypeError("Invalid type for Atom2.")

        if Mass1 is None:
            Mass1 = np.ones(len(self.Atoms1), float)
        else:
            if len(Mass1) != len(self.Atoms1):
                raise ValueError("Length of Mass must be same as number of Atoms1")
        Mass1 = np.array(Mass1, float)
        
        #make weights
        if Kind == KindCOM:
            self.Weight = Mass1 / Mass1.sum()
        elif Kind == KindCentroid:
            self.Weight = np.ones(len(self.Atoms1), float) / float(len(self.Atoms1))
        else:
            raise ValueError("Kind of map not recognized.")



    def Update(self, Pos1, Pos2):
        """Updates Pos2 based on positions of atoms in Pos1."""
        Pos2[self.Atom2,:] = np.sum(Pos1[self.Atoms1,:] * self.Weight[:,np.newaxis], axis=0)
        



class PosMap(list):

    def __init__(self, Warn = True):
        """Returns a mapping list for translating Pos1 to Pos2."""
        list.__init__(self)
        self.Warn = Warn
        self.__AtomMap1 = None
        self.NAtom2 = None
        
    def Add(self, *args, **kwargs):
        am = AtomMap(*args, **kwargs)
        self.append(am)
        self.__AtomMap1 = None

    def Check(self):
        """Checks to see if definitions established for all atoms."""
        NotUpdated = []
        for am in self:
            n = len(NotUpdated)
            if (n < am.Atom2 + 1):
                NotUpdated = NotUpdated + [True] * (am.Atom2 + 1 - n)
            if NotUpdated[am.Atom2]:
                NotUpdated[am.Atom2] = False
            else:
                raise IndexError("Found multiple maps for Atom2 %d" % am.Atom2)
        if any(NotUpdated):
            s = "Warning: found no updates for Atom2 "
            s += ",".join([str(i) for (i,v) in enumerate(NotUpdated) if v])
            print s

    def Print(self):
        """Prints out mappings."""
        self.Prepare()
        for (i, Ind1) in enumerate(self.__StartInd1[:-1]):
            Ind2 = self.__StartInd1[i+1]
            print "CG atom %d contains all atoms:" % i
            for Ind in range(Ind1, Ind2):
                print "  %d (weight %f)" % (self.__AtomMap1[Ind], self.__Weight[Ind])
            

    def Prepare(self):
        """Updates internal variables for fast use."""
        self.sort(cmp = lambda x,y: cmp(x.Atom2, y.Atom2))
        AtomMap1 = []
        StartInd1 = []
        Weight = []
        Atom2 = -1
        for am in self:
            while Atom2 < am.Atom2:
                StartInd1.append(len(AtomMap1))
                Atom2 += 1
            Weight.extend(am.Weight)
            AtomMap1.extend(am.Atoms1)
        StartInd1.append(len(AtomMap1))
        self.__AtomMap1 = np.array(AtomMap1, int)
        self.__StartInd1 = np.array(StartInd1, int)
        self.__Weight = np.array(Weight, float)
        self.MaxAtom1 = np.max(self.__AtomMap1) - 1
        self.NAtom2 = Atom2 + 1
        if self.Warn: self.Check()


    def Map(self, Pos1, BoxL):
        """Returns Pos2 based on positions of atoms in Pos1."""
        if self.__AtomMap1 is None:
            self.Prepare()
        if len(Pos1) < self.MaxAtom1:
            raise ValueError("Pos1 has %d atoms but atom map expects %d" % (len(Pos1), self.MaxAtom1))
        return Lib.Module.map(self.__AtomMap1, self.__Weight, self.__StartInd1, self.NAtom2, Pos1, BoxL)
                                     
    



            
            
    