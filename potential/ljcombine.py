#/usr/bin/env python


import numpy as np
import base

import sim.atomselect as atomselect
import sim.chem as chem


class LJCombine(base.PotentialClass):
    """Lennard-Jones interactions with LB combining rules."""

    Names = ["ljcombine"]

    Source = """
>>> defs
float idist2
float idist6
float idist12
float Sigma
float Eps
float HSigmai
float SqrtEpsi
float HSigmaj
float SqrtEpsj
float rval1
float rval2

>>> ljcombineenergydparam
Sigma = HSigmai + HSigmaj
Eps = SqrtEpsi * SqrtEpsj
idist2 = Sigma * Sigma / DIJSQ
idist6 = idist2 * idist2 * idist2
idist12 = idist6 * idist6
THISU = Eps * (4. * (idist12 - idist6) + UShift0(AIDIJ))
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    THISW = Eps * (24. * idist6 - 48. * idist12)
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM) then
    rval1 = 1. / Sigma
    if (AIDI == AIDJ) then
        rval2 = (48. * idist12 - 24. * idist6) * rval1 + UShift1(AIDIJ)
        DU_HSigma(AIDI) = 2. * rval2 * Eps
        DDU_HSigma_SqrtEps(AIDI,AIDI) = 4. * rval2 * SqrtEpsi
        rval2 = Eps * ((528.*idist12 - 120.*idist6) * rval1*rval1 + UShift2(AIDIJ))
        DDU_HSigma_HSigma(AIDI,AIDI) = 4. * rval2
        rval2 = 4. * (idist12 - idist6) + UShift0(AIDIJ)
        DU_SqrtEps(AIDI) = 2. * rval2 * SqrtEpsi
        DDU_SqrtEps_SqrtEps(AIDI,AIDI) = 2. * rval2
    else
        rval2 = (48. * idist12 - 24. * idist6) * rval1 + UShift1(AIDIJ)
        DU_HSigma(AIDI) = rval2 * Eps
        DU_HSigma(AIDJ) = rval2 * Eps
        DDU_HSigma_SqrtEps(AIDI,AIDI) = rval2 * SqrtEpsj
        DDU_HSigma_SqrtEps(AIDI,AIDJ) = rval2 * SqrtEpsi
        DDU_HSigma_SqrtEps(AIDJ,AIDI) = rval2 * SqrtEpsj
        DDU_HSigma_SqrtEps(AIDJ,AIDJ) = rval2 * SqrtEpsi
        rval2 = Eps * ((528.*idist12 - 120.*idist6) * rval1*rval1 + UShift2(AIDIJ))
        DDU_HSigma_HSigma(AIDI,AIDJ) = rval2
        DDU_HSigma_HSigma(AIDI,AIDI) = rval2
        DDU_HSigma_HSigma(AIDJ,AIDJ) = rval2
        rval2 = 4. * (idist12 - idist6) + UShift0(AIDIJ)
        DU_SqrtEps(AIDI) = rval2 * SqrtEpsj
        DU_SqrtEps(AIDJ) = rval2 * SqrtEpsi
        DDU_SqrtEps_SqrtEps(AIDI,AIDJ) = rval2
    endif
endif
if (CALCDWPARAM) then
    rval1 = 1. / Sigma
    idist12 = -12. * idist12
    idist6 = -6. * idist6
    if (AIDI == AIDJ) then
        rval2 = (48. * idist12 - 24. * idist6) * rval1
        DW_HSigma(AIDI) = 2. * rval2 * Eps
        DDW_HSigma_SqrtEps(AIDI,AIDI) = 4. * rval2 * SqrtEpsi
        rval2 = Eps * (528.*idist12 - 120.*idist6) * rval1*rval1 
        DDW_HSigma_HSigma(AIDI,AIDI) = 4. * rval2
        rval2 = 4. * (idist12 - idist6) 
        DW_SqrtEps(AIDI) = 2. * rval2 * SqrtEpsi
        DDW_SqrtEps_SqrtEps(AIDI,AIDI) = 2. * rval2
    else
        rval2 = (48. * idist12 - 24. * idist6) * rval1 
        DW_HSigma(AIDI) = rval2 * Eps
        DW_HSigma(AIDJ) = rval2 * Eps
        DDW_HSigma_SqrtEps(AIDI,AIDI) = rval2 * SqrtEpsj
        DDW_HSigma_SqrtEps(AIDI,AIDJ) = rval2 * SqrtEpsi
        DDW_HSigma_SqrtEps(AIDJ,AIDI) = rval2 * SqrtEpsj
        DDW_HSigma_SqrtEps(AIDJ,AIDJ) = rval2 * SqrtEpsi
        rval2 = Eps * (528.*idist12 - 120.*idist6) * rval1*rval1 
        DDW_HSigma_HSigma(AIDI,AIDJ) = rval2
        DDW_HSigma_HSigma(AIDI,AIDI) = rval2
        DDW_HSigma_HSigma(AIDJ,AIDJ) = rval2
        rval2 = 4. * (idist12 - idist6) 
        DW_SqrtEps(AIDI) = rval2 * SqrtEpsj
        DW_SqrtEps(AIDJ) = rval2 * SqrtEpsi
        DDW_SqrtEps_SqrtEps(AIDI,AIDJ) = rval2
    endif
endif

>>> mainloopbeforepair
HSigmai = HSigma(AIDI)
SqrtEpsi = SqrtEps(AIDI)

>>> mainlooppair
!pair interactions for potential %(PName)s
HSigmaj = HSigma(AIDJ)
SqrtEpsj = SqrtEps(AIDJ)
[ljcombineenergydparam]
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
call GetijFromPairInd(AIDIJ, AIDi, AIDj)
HSigmai = HSigma(AIDI)
SqrtEpsi = SqrtEps(AIDI)
HSigmaj = HSigma(AIDJ)
SqrtEpsj = SqrtEps(AIDJ)

>>> argeval
DIJ = ARGVAL
DIJSQ = DIJ * DIJ
[ljcombineenergydparam]
"""

    TestArgs = {"Cut":2.5, "Shift":True,
                "Sigma":np.array([0.7,0.9,0.8]),
                "Epsilon":np.array([2.3,1.7,1.1])}
    Type = base.ptypes.PairPotential
    UsesATypes = True

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Shift = False, Sigma = None, Epsilon = None):
        """Initializes a Lennard-Jones w/ LB combining potential."""
        NType, TypeAIDs, TypeLabels = Sys.World.GetPairTypes()
        base.PotentialClass.__init__(self, Sys, Label, Cut, Filter,
                                     NType = NType, TypeLabels = TypeLabels,
                                     Source = LJCombine.Source)
        self.ModuleVars += ["UShift0", "UShift1", "UShift2"]

        self.Shift = Shift
        if Sigma is None:
            Sigma = Sys.Units.LScale
        elif type(Sigma) is list:
            Sigma = np.array(Sigma, dtype=float)
        if Epsilon is None:
            Epsilon = Sys.Units.EScale
        elif type(Epsilon) is list:
            Epsilon = np.array(Epsilon, dtype=float)            
        self.Param.Add("SqrtEps", Sys.World.NAID, Value = np.sqrt(Epsilon),
                       Min = 0., Fixed = Fixed,
                       Scale = np.sqrt(self.Sys.Units.EScale))
        self.Param.Add("HSigma", Sys.World.NAID, Value = Sigma / 2.,
                       Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.UShift0 = np.zeros(self.Arg.NType, dtype=float)
        self.UShift1 = np.zeros(self.Arg.NType, dtype=float)
        self.UShift2 = np.zeros(self.Arg.NType, dtype=float)

        self.Update()   
        self.SetTypeInd()

    def SetParam(self, HSigma = None, SqrtEps = None):
        """Sets parameters for a specific atom filter."""
        if not HSigma is None:
            self.HSigma = HSigma
        if not SqrtEps is None:
            self.SqrtEps = SqrtEps
        self.Update()
        
    def GetTypeInd(self, Type1, Type2):
        """Returns the type number for this potential."""
        return chem.PairID(Type1.AID, Type2.AID)
        
    def SetTypeInd(self, TypeInd = None):
        """Sets the type to be used in Val and DVal, or None for generic."""
        if TypeInd is None:
            AID1, AID2 = 0, 0
        else:
            AID1, AID2 = chem.GetPair(TypeInd)
        self.CurrentSig = self.HSigma[AID1] + self.HSigma[AID2]
        self.CurrentEps = self.SqrtEps[AID1] * self.SqrtEps[AID2]
        self.CurrentShift = self.UShift0[chem.PairID(AID1,AID2)]

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        else:
            r = x / self.CurrentSig 
            return self.CurrentEps * (4. * r**(-12) - 4. * r**(-6) + self.CurrentShift)
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        else:
            r = x / self.CurrentSig 
            return 4. * self.CurrentEps * self.CurrentSig * (-12. * r**(-13) + 6. * r**(-7))

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        M = zip(self.SqrtEps, self.HSigma)
        for (i, (SqrtEpsi, HSigmai)) in enumerate(M):
            HSigmai.SoftMin = (0.8 * HSigmai.ArgMin) / 2.
            SqrtEpsi.SoftMin = np.sqrt(0.001 * self.Sys.Units.EScale)
            if not MaxFracChange is None:
                HSigmai.MaxChange = self.Sys.Units.LScale * MaxFracChange
                SqrtEpsi.MaxChange = np.sqrt(self.Sys.Units.EScale) * MaxFracChange

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            M = zip(self.SqrtEps, self.HSigma)
            for (i, (SqrtEpsi, HSigmai)) in enumerate(M):
                for (j, (SqrtEpsj, HSigmaj)) in enumerate(M):                  
                    Eps = SqrtEpsi * SqrtEpsj
                    Sigma = HSigmai + HSigmaj                   
                    k = chem.PairID(i,j)
                    if Sigma <= 0. or Eps <= 0.:
                        self.UShift0[k] = 0.
                        self.UShift1[k] = 0.
                        self.UShift2[k] = 0.
                    else:
                        rc = self.Cut / Sigma
                        self.UShift0[k] = -4. * (rc**(-12) - rc**(-6))
                        self.UShift1[k] = -(48. * rc**(-12) - 24. * rc**(-6)) / Sigma
                        self.UShift2[k] = -(528. * rc**(-12) - 120.*rc**(-6)) / Sigma**2
        else:
            self.UShift0.fill(0.)
            self.UShift1.fill(0.)
            self.UShift2.fill(0.)

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        raise ValueError("This function not yet implemented.")
        

 