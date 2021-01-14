#/usr/bin/env python

import numpy as np
import base
import sim.chem as chem


class ScreenedCoulomb(base.PotentialClass):
    """Screened Coulombic interactions."""

    Names = ["screenedcoulomb"]

    Source = """
>>> defs
float Chargei
float Chargej
float idist
float val1
float val2
float val3
float val4
float val5


>>> screenedcoulombenergydparam
idist = 1.d0 / DIJ
val1 = Chargei * Chargej
val2 = exp(-Kappa(0) * DIJ) * idist
val3 = Coef(0) * val1
THISU = val3 * (val2 + UShift(0))
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    THISW = -val3 * val2 * (1.d0 + Kappa(0) * DIJ)
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM) then
    DU_Coef(0) = val1 * (val2 + UShift(0))
    DU_Kappa(0) = -val3 * (val2 * DIJ + UShift(1))
    DDU_Kappa_Kappa(0,0) = val3 * (val2 * DIJSQ + UShift(2))
    DDU_Kappa_Coef(0,0) = -val1 * (val2 * DIJ + UShift(1))
    if (AIDI==AIDJ) then
        DU_Charge(AIDI) = 2.d0 * Coef(0) * Chargej * (val2 + UShift(0))
        DDU_Charge_Kappa(AIDI,0) = -2.d0 * Coef(0) * Chargej * (val2 * DIJ + UShift(1))
        DDU_Charge_Coef(AIDI,0) = 2.d0 * Chargej * (val2 + UShift(0))
        DDU_Charge_Charge(AIDI,AIDJ) = 2.d0 * Coef(0) * (val2 + UShift(0))
    else
        DU_Charge(AIDI) = Coef(0) * Chargej * (val2 + UShift(0))
        DU_Charge(AIDJ) = Coef(0) * Chargei * (val2 + UShift(0))
        DDU_Charge_Kappa(AIDI,0) = -Coef(0) * Chargej * (val2 * DIJ + UShift(1))
        DDU_Charge_Kappa(AIDJ,0) = -Coef(0) * Chargei * (val2 * DIJ + UShift(1))
        DDU_Charge_Coef(AIDI,0) = Chargej * (val2 + UShift(0))
        DDU_Charge_Coef(AIDJ,0) = Chargei * (val2 + UShift(0))
        DDU_Charge_Charge(AIDI,AIDJ) = Coef(0) * (val2 + UShift(0))
    endif
endif
if (CALCDWPARAM) then
    val5 = 1.d0 + Kappa(0) * DIJ
    DW_Coef(0) = -val1 * val2 * val5
    val4 = val3 * val2 * DIJSQ
    DW_Kappa(0) = val4 * Kappa(0)
    DDW_Kappa_Kappa(0,0) = val4 * (1 - DIJ * Kappa(0))
    DDW_Kappa_Coef(0,0) = val1 * val2 * DIJSQ * Kappa(0)
    if (AIDI==AIDJ) then
        DW_Charge(AIDI) = -2.d0 * Coef(0) * Chargej * val2 * val5
        DDW_Charge_Kappa(AIDI,0) = 2.d0 * Coef(0) * Chargej * val2 * DIJSQ * Kappa(0)
        DDW_Charge_Coef(AIDI,0) = -2.d0 * Chargej * val2 * val5
        DDW_Charge_Charge(AIDI,AIDJ) = -2.d0 * Coef(0) * val2 * val5
    else
        DW_Charge(AIDI) = -Coef(0) * Chargej * val2 * val5
        DW_Charge(AIDJ) = -Coef(0) * Chargei * val2 * val5
        val4 = Coef(0) * val2 * DIJSQ * Kappa(0)
        DDW_Charge_Kappa(AIDI,0) = val4 * Chargej 
        DDW_Charge_Kappa(AIDJ,0) = val4 * Chargei 
        DDW_Charge_Coef(AIDI,0) = -Chargej * val2  * val5
        DDW_Charge_Coef(AIDJ,0) = -Chargei * val2 * val5
        DDW_Charge_Charge(AIDI,AIDJ) = -Coef(0) * val2 * val5
    endif
endif

>>> mainloopbeforepair
Chargei = Charge(AIDI)

>>> mainlooppair
!pair interactions for potential %(PName)s
Chargej = Charge(AIDJ)
[screenedcoulombenergydparam]
if (CALCFORCE) then
    FORCEI = RIJ * THISW / DIJSQ
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
Chargei = Charge(AIDI)
Chargej = Charge(AIDJ)

>>> argeval
DIJ = ARGVAL
DIJSQ = DIJ * DIJ
[screenedcoulombenergydparam]
"""

    TestArgs = {"Cut":2.5, "Kappa":2., "Coef":1.1, "Shift":True}
    Type = base.ptypes.PairPotential
    UsesCharge = True
    UsesATypes = True

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False,
                 FixedCoef = True, FixedKappa = False,
                 Kappa = None, Coef = None, Shift = False):
        """Initializes a screened Coulomb potential."""
        NType, TypeAIDs, TypeLabels = Sys.World.GetPairTypes()
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Filter,
                                     NType = NType, TypeLabels = TypeLabels,
                                     Source = ScreenedCoulomb.Source)
        self.ModuleVars += ["UShift"]

        self.Shift = Shift
        #must fix coef if charges vary
        self.FixedCoef = FixedCoef
        self.FixedKappa = FixedKappa

        if Kappa is None:
            Kappa = 1. / Sys.Units.LScale
        if Coef is None:
            Coef = 1. / Sys.Units.FPE
        self.Param.Add("Coef", 1, Value = Coef, Min = 0., Fixed = Fixed or FixedCoef,
                       Scale = 1. / self.Sys.Units.FPE)
        self.Param.Add("Kappa", 1, Value = Kappa, Min = 0., Fixed = Fixed or FixedKappa,
                       Scale = 1. / self.Sys.Units.LScale)
        self.UShift = np.zeros(3, dtype=float)
        
        self.Update()
        self.SetTypeInd()
        

    def SetParam(self, Kappa = None, Coef = None):
        """Sets parameters."""
        if not Kappa is None:
            if self.FixedKappa:
                raise ValueError("Kappa is Fixed.")
            self.Kappa = Kappa
        if not Coef is None:
            if self.FixedCoef:
                raise ValueError("Coefficient is Fixed.")
            self.Coef = Coef
        self.Update()
        
        
    def GetTypeInd(self, Type1, Type2):
        """Returns the type number for this potential."""
        return chem.PairID(Type1.AID, Type2.AID)
        
    def SetTypeInd(self, TypeInd = None):
        """Sets the type to be used in Val and DVal, or None for generic."""
        if TypeInd is None:
            self.CurrentChargeCoef = self.Coef[0] * self.Sys.Units.ECharge**2
        else:
            AID1, AID2 = chem.GetPair(TypeInd)
            Charge1 = self.Sys.World.AtomTypes[AID1].Charge 
            Charge2 = self.Sys.World.AtomTypes[AID2].Charge 
            self.CurrentChargeCoef = self.Coef[0] * Charge1 * Charge2      

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        return self.CurrentChargeCoef * np.exp(-self.Kappa[0] * x) / x
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        return -self.CurrentChargeCoef * np.exp(-self.Kappa[0] * x) * (1+self.Kappa[0]*x) / x**2            

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Coef.SoftMin = 0.0001 / self.Sys.Units.FPE
        self.Kappa.SoftMin = 0.00001 / self.Sys.Units.LScale
        if not MaxFracChange is None:
            self.Coef.MaxChange = MaxFracChange / self.Sys.Units.FPE
            self.Kappa.MaxChange = MaxFracChange / self.Sys.Units.LScale
        

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            self.UShift[0] = -np.exp(-self.Kappa[0] * self.Cut) / self.Cut
            self.UShift[1] = -np.exp(-self.Kappa[0] * self.Cut)
            self.UShift[2] = -np.exp(-self.Kappa[0] * self.Cut) * self.Cut
        else:
            self.UShift = 0.

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Kappa = None
        Coef = None
        if not self.FixedCoef: Coef = 1. / self.Sys.Units.FPE
        if not self.FixedKappa: Kappa = 1. / self.Kappa.ArgAvg
        self.SetParam(Kappa, Coef)
        


 