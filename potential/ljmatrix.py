#/usr/bin/env python


import numpy as np
import base

import sim.atomselect as atomselect
import sim.chem as chem


class LJMatrix(base.PotentialClass):
    """Lennard-Jones interactions with matrices for A and B."""

    Names = ["ljmatrix"]

    Source = """
>>> defs
float idist2
float idist6
float idist12

>>> ljmatrixenergydparam
idist2 = 1. / DIJSQ
idist6 = idist2 * idist2 * idist2
idist12 = idist6 * idist6
THISU = LJA(AIDIJ) * idist12 - LJB(AIDIJ) * idist6 + UShift0(AIDIJ)
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    THISW = -12. * LJA(AIDIJ) * idist12 + 6. * LJB(AIDIJ) * idist6
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM) then
    DU_LJA(AIDIJ) = idist12 + UShift1(AIDIJ)
    DU_LJB(AIDIJ) = -idist6 + UShift2(AIDIJ)
endif
if (CALCDWPARAM) then
    DW_LJA(AIDIJ) = -12. * idist12
    DW_LJB(AIDIJ) = 6. * idist6
endif

>>> mainlooppair
!pair interactions for potential %(PName)s
[ljmatrixenergydparam]
if (CALCFORCE) then
    FORCEI = (RIJ * THISW / DIJSQ) * SCALE
    FORCE(I,:) = FORCE(I,:) + FORCEI
    FORCE(J,:) = FORCE(J,:) - FORCEI
endif

>>> argmainlooppair
ARGVAL = DIJ
ARGTYPE = AIDIJ
[ARGGET]

>>> argevalconst
AIDIJ = ArgType

>>> argeval
DIJ = ARGVAL
DIJSQ = DIJ * DIJ
[ljmatrixenergydparam]
"""

    TestArgs = {"Cut":2.5, "Shift":True,
                "LJA":np.array([0.8, 0.9, 1.0, 1.1, 1.2, 1.3]),
                "LJB":np.array([0.9, 1.2, 1.0, 0.8, 1.1, 1.3])}
    Type = base.ptypes.PairPotential
    UsesATypes = True

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Shift = False, LJA = None, LJB = None):
        """Initializes a Lennard-Jones with interaction matrices for A and B."""
        NType, TypeAIDs, TypeLabels = Sys.World.GetPairTypes()
        base.PotentialClass.__init__(self, Sys, Label, Cut, Filter,
                                     NType = NType, TypeLabels = TypeLabels,
                                     Source = LJMatrix.Source)
        self.ModuleVars += ["UShift0", "UShift1", "UShift2"]

        self.Shift = Shift
        
        NAID = Sys.World.NAID
        NType = self.Arg.NType
             
        
        def CheckMatrix(M, Exp, Name):
            ret = np.ones(NType, dtype = float)
            if M is None:
                return Sys.Units.LScale**(-Exp) * ret
            elif type(M) is float:
                return M * ret
            elif type(M) is np.ndarray:
                if M.ndim == 1:
                    if len(M) == NType:
                        return M.copy()
                    else:
                        raise ValueError("%s is wrong length for 1D version." % Name)
                elif M.ndim == 2:
                    if not len(M) == NAID:
                        raise ValueError("%s is wrong length for 2D version." % Name)
                    if not np.all(M.T == M):
                        raise ValueError("%s is not a symmetric matrix." % Name)
                    for i in range(NType):
                        for j in range(NType):
                            ret[chem.PairID(i, j)] = M[i,j]
                    return ret
            else:
                raise ValueError("%s type is not recognized." % Name)        
            
        LJA = CheckMatrix(LJA, 12, "LJA")
        LJB = CheckMatrix(LJB, 6, "LJB")
        
        self.Param.Add("LJA", NType, Value = LJA, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale**12)
        self.Param.Add("LJB", NType, Value = LJB, Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale**6)
        self.UShift0 = np.zeros(NType, dtype=float)
        self.UShift1 = np.zeros(NType, dtype=float)
        self.UShift2 = np.zeros(NType, dtype=float)
        self.Update()   
        self.SetTypeInd()

    def SetParam(self, LJA = None, LJB = None):
        """Sets parameters for a specific atom filter."""
        if not LJA is None:
            self.LJA = LJA
        if not LJB is None:
            self.LJB = LJB
        self.Update()
        
    def GetTypeInd(self, Type1, Type2):
        """Returns the type number for this potential."""
        return chem.PairID(Type1.AID, Type2.AID)
        
    def SetTypeInd(self, TypeInd = 0):
        """Sets the type to be used in Val and DVal, or None for generic."""
        self.CurrentLJA = self.LJA[TypeInd]
        self.CurrentLJB = self.LJB[TypeInd] 
        self.CurrentShift = self.UShift0[TypeInd]

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        else:
            return self.CurrentLJA * x**(-12) - self.CurrentLJB * x**(-6) + self.CurrentShift
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        else:
            return -12. * self.CurrentLJA * x**(-13) + self.CurrentLJB * x**(-7)

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        M = zip(self.LJA, self.LJB)
        for (i, (A, B)) in enumerate(M):
            A.SoftMin = 0.001 * self.Sys.Units.EScale * (0.8 * A.ArgMin)**12 
            if not MaxFracChange is None:
                A.MaxChange = self.Sys.Units.EScale * (self.Sys.Units.LScale**12) * MaxFracChange
                B.MaxChange = self.Sys.Units.EScale * (self.Sys.Units.LScale**6) * MaxFracChange

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            rc = self.Cut
            self.UShift0[:] = -self.LJA * rc**(-12) + self.LJB * rc**(-6)
            self.UShift1[:] = -rc**(-12)
            self.UShift2[:] = rc**(-6)
        else:
            self.UShift0.fill(0.)
            self.UShift1.fill(0.)
            self.UShift2.fill(0.)

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        raise ValueError("This function not yet implemented.")
        

 