#/usr/bin/env python


### Routines for initial positions in SIM suite.
### coded by MSS

import numpy as np




def GetL(Sys, L = None):
    """Returns box lengths."""
    if L is None:
        L = Sys.BoxL
        if any(L == 0.):
            raise TypeError("System lacks periodicity in one or more directions; must specify L")
    if not hasattr(L, "__getitem__"):
        L = [L for x in range(Sys.Dim)]
    return L   
    
    
def CenterOnMol(Sys, Mol=0):
    """Centers the simulation box such that the centroid of molecule Mol is in the middle.
Sys: A System instance
Mol: Either a molecule index or object"""
    if type(Mol) is int:
        Mol = Sys.Mol[Mol]
    Center = np.zeros(Sys.Dim, dtype=float)
    LastPos = np.zeros(Sys.Dim, dtype=float)
    for Atom in Mol:
        ThisPos = LastPos + Sys.MinImage(Atom.Pos - LastPos)
        Center += ThisPos
        LastPos = ThisPos
    Center = Center / len(Mol)
    Sys.Pos = Sys.Pos - Center
            

def CubicLattice(Sys, L = None, Random = 0.):
    """Places all atom / molecule positions on a cubic lattice."""
    L = GetL(Sys, L)
    #compute integer grid # of locations for cubic lattice
    Units = [x for x in Sys.RigidUnitsIter()]
    N = len(Units)
    NLat = int(N**(1./Sys.Dim) + 1.)
    #make an array of lattice sites
    r = (np.arange(NLat, dtype=float) + 0.1) / NLat - 0.5
    #loop through positions in lattice until done
    #for every atom in the system
    ind = np.zeros(Sys.Dim, int)
    i = 0
    Orient = np.zeros(Sys.Dim, float)
    while i < N:
        Pos = L * r.take(ind) + Random * (np.random.rand(Sys.Dim) - 0.5)
        Sys.SetPos(Units[i], Pos, Orient)
        #increment
        ind[0] += 1
        for j in range (Sys.Dim-1):
            if ind[j] >= NLat:
                ind[j] = 0
                ind[j+1] = ind[j+1] + 1
        #if done placing atoms, return
        i += 1
        if i >= N:
            return 
            
def CubicLatticeFill(Sys, L = None, Random = 0.):
    """Places all atom / molecule positions on a cubic lattice;
number of lattice sites on each axis scaled by box length."""
    L = GetL(Sys, L)
    #compute integer grid # of locations for cubic lattice
    Units = [x for x in Sys.RigidUnitsIter()]
    N = len(Units)
    GridSize = (np.product(L) / N)**(1./Sys.Dim)
    NLat = (L / GridSize + 1.).astype(int)
    #make an array of lattice sites
    r = []
    for i in range(Sys.Dim):
        r.append((np.arange(NLat[i], dtype=float) + 0.1) / NLat[i] - 0.5)
    #loop through positions in lattice until done
    #for every atom in the system
    ind = np.zeros(Sys.Dim, int)
    i = 0
    Orient = np.zeros(Sys.Dim, float)
    while i < N:
        Pos = np.array([thisr[j] for (thisr, j) in zip(r, ind)], float)
        Pos = L * Pos + Random * (np.random.rand(Sys.Dim) - 0.5)
        Sys.SetPos(Units[i], Pos, Orient)
        #increment
        ind[0] += 1
        for j in range (Sys.Dim-1):
            if ind[j] >= NLat[j]:
                ind[j] = 0
                ind[j+1] = ind[j+1] + 1
        #if done placing atoms, return
        i += 1
        if i >= N:
            return 
