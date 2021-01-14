#/usr/bin/env python
#Added 2019.05.15 by Kevin Shen
#To simulate regularized smeared Coulomb interaction
##`a` parameter is defined s.t. reproduces Born energy lb/2a
#Interaction is erf(sqrt(pi)r/2a)/r
#Right now assumes that ewald choice of alpha and rcut is chosen to appropriately reproduce 1/r
#Further, default is to use a shift, such that the rcut and shift errors incurred by Ewald are immaterial (i.e. if know the Ewald shift, those shifts will just be substracted off again by the shift of this correction potential)
#
#TODO: allow for species-specific ion sizes, will need to retrieve appropriate a-size and add necessary chain derivatives

import numpy as np
from scipy.special import erf
import base
import sim.chem as chem


class SmearedCoulombEwCorr(base.PotentialClass):
    """Smeared Coulombic interaction correction for Ewald"""

    Names = ["smearedcoulombEwCorr"]

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
float val6
float fac
float fac2
float iA
float tmp

>>> smearedcoulombEwCorrenergydparam
idist = 1.d0 / DIJ
fac  = sqrt(4.d0*atan(1.0_8))/2.d0/BornA(0)
fac2 = fac*fac
val1 = Chargei * Chargej
external :: erf
call erf(fac*DIJ, tmp)
val2 = tmp * idist
val3 = Coef(0) * val1
THISU = val3 * (val2 - idist + UShift(0) + CoulShift(0))
[UPDATEENERGY]
iA = 1.d0 / BornA(0)
val5 = exp(-fac2*DIJSQ)*iA*iA
val6 = exp(-fac2*DIJSQ)/BornA(0)
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    THISW = val3 * ( val6 - val2 + idist)   
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM) then
    DU_Coef(0) = val1 * (val2 - idist + UShift(0) + CoulShift(0))
    DU_BornA(0) =  val3*(-val5 + UShift(1)  ) 
    DDU_BornA_BornA(0,0) = val3 * ( val6*2.d0/BornA(0)**2.d0 - val6*2.d0/BornA(0)**2.d0*fac2*DIJSQ  + UShift(2)  )   
    DDU_BornA_Coef(0,0) = val1*(-val5 + UShift(1) ) 
    if (AIDI==AIDJ) then
        DU_Charge(AIDI) = 2.d0 * Coef(0) * Chargej * (val2 + UShift(0) - idist + CoulShift(0))
        DDU_Charge_BornA(AIDI,0) = 2.d0 * Coef(0) * Chargej * (-val5 + UShift(1))
        DDU_Charge_Coef(AIDI,0) = 2.d0 * Chargej * (val2 + UShift(0) - idist + CoulShift(0))
        DDU_Charge_Charge(AIDI,AIDJ) = 2.d0 * Coef(0) * (val2 + UShift(0) - idist + CoulShift(0))
    else
        DU_Charge(AIDI) = Coef(0) * Chargej * (val2 + UShift(0) -idist + CoulShift(0))
        DU_Charge(AIDJ) = Coef(0) * Chargei * (val2 + UShift(0) -idist + CoulShift(0))
        DDU_Charge_BornA(AIDI,0) = Coef(0) * Chargej * (-val5 + UShift(1))
        DDU_Charge_BornA(AIDJ,0) = Coef(0) * Chargei * (-val5 + UShift(1))
        DDU_Charge_Coef(AIDI,0) = Chargej * (val2 + UShift(0) - idist + CoulShift(0))
        DDU_Charge_Coef(AIDJ,0) = Chargei * (val2 + UShift(0) - idist + CoulShift(0))
        DDU_Charge_Charge(AIDI,AIDJ) = Coef(0) * (val2 + UShift(0) - idist + CoulShift(0))
    endif
endif
if (CALCDWPARAM) then

    DW_Coef(0) = val1 * (val6 - val2 + idist)
    val4 = 2.d0 * val5 * fac2 * DIJSQ
    DW_BornA(0) = val3*val4

    DDW_BornA_BornA(0,0) = 2.d0*iA * (-2.d0 + fac2*DIJSQ) * val3*val4 
    DDW_BornA_Coef(0,0) = val1 * val4

    if (AIDI==AIDJ) then
        DW_Charge(AIDI) = 2.d0 * Chargej * Coef(0) * (val6 - val2 + idist)
        DDW_Charge_BornA(AIDI,0) = 2.d0 * Chargej * Coef(0) * val4
        DDW_Charge_Coef(AIDI,0) = 2.d0 * Chargej * (val6 - val2 + idist)
        DDW_Charge_Charge(AIDI,AIDJ) = 2.d0 * Coef(0) * (val6 - val2 + idist)
    else
        DW_Charge(AIDI) = Coef(0) * Chargej * (val6 - val2 + idist)
        DW_Charge(AIDJ) = Coef(0) * Chargei * (val6 - val2 + idist)
        
        DDW_Charge_BornA(AIDI,0) = Coef(0) * val4 * Chargej 
        DDW_Charge_BornA(AIDJ,0) = Coef(0) * val4 * Chargei 
        DDW_Charge_Coef(AIDI,0) = Chargej * (val6 - val2 + idist)
        DDW_Charge_Coef(AIDJ,0) = Chargei * (val6 - val2 + idist)
        DDW_Charge_Charge(AIDI,AIDJ) = Coef(0) * (val6 - val2 + idist)
    endif
endif

>>> mainloopbeforepair
Chargei = Charge(AIDI)

>>> mainlooppair
!pair interactions for potential %(PName)s
Chargej = Charge(AIDJ)
[smearedcoulombEwCorrenergydparam]
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
[smearedcoulombEwCorrenergydparam]
"""

    TestArgs = {"Cut":2.5, "BornA":1., "Coef":7.0, "Shift":True}
    Type = base.ptypes.PairPotential
    UsesCharge = True
    UsesATypes = True

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False,
                 FixedCoef = True, FixedBornA = False,
                 BornA = None, Coef = None, Shift = True):
        """Initializes a smeared Coulomb potential correction for Ewald."""
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Filter,
                                     NType = chem.NPairID(Sys.World.NAID),
                                     Source = SmearedCoulombEwCorr.Source)
        self.ModuleVars += ["UShift","CoulShift"]
        self.Shift = Shift
        #must fix coef if charges vary
        self.FixedCoef = FixedCoef
        self.FixedBornA = FixedBornA

        if BornA is None:
            BornA = Sys.Units.LScale
        if Coef is None:
            Coef = 1. / Sys.Units.FPE
        self.Param.Add("Coef", 1, Value = Coef, Min = 0., Fixed = Fixed or FixedCoef,
                       Scale = 1. / self.Sys.Units.FPE)
        self.Param.Add("BornA", 1, Value = BornA, Min = 0.0, Fixed = Fixed or FixedBornA,
                       Scale = 1. / self.Sys.Units.LScale)
        #self.Const.Add("UShift", 3, Value = 0.)
        self.UShift = np.zeros(3, dtype=float)
        self.CoulShift = np.zeros(3, dtype=float)

        self.Update()
        self.SetTypeInd()
        

    def SetParam(self, BornA = None, Coef = None):
        """Sets parameters."""
        if not BornA is None:
            if self.FixedBornA:
                raise ValueError("BornA is Fixed.")
            self.BornA = BornA
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
        if x == 0.:
            #print('x==0, val is {}'.format(self.currentChargeCoef /self.BornA[0]))
            return self.CurrentChargeCoef /self.BornA[0]
        #print('ChargeCoef {}, x {}'.format(self.CurrentChargeCoef,x))
        v = self.CurrentChargeCoef * (erf(np.sqrt(np.pi)*x/2/self.BornA[0]) -1)/ x + self.CurrentChargeCoef*(self.UShift[0] + self.CoulShift[0])
        #print('val is {},{}'.format(self.CurrentChargeCoef/x, v))
        return v
        #return self.CurrentChargeCoef * np.exp(-self.Kappa[0] * x) / x
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        if x == 0.:
            return 0.
        fac = np.sqrt(np.pi)/2/self.BornA[0]
        return self.CurrentChargeCoef * ( np.exp(-fac**2 * x**2) / (x*self.BornA[0]) - erf(fac*x)/x**2 + 1/x**2)
        #return -self.CurrentChargeCoef * np.exp(-self.Kappa[0] * x) * (1+x*self.Kappa[0]) / x**2  #FIXED KS 2019.03.12         

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Coef.SoftMin = 0.0001 / self.Sys.Units.FPE
        self.BornA.SoftMin = 0.001 / self.Sys.Units.LScale
        if not MaxFracChange is None:
            self.Coef.MaxChange = MaxFracChange / self.Sys.Units.FPE
            self.BornA.MaxChange = MaxFracChange / self.Sys.Units.LScale
        

    def Update(self):
        """Updates the potential."""
        if self.BornA is None:
            self.BornA = 1.0
        if self.Shift:
            fac = np.sqrt(np.pi)/2.0/self.BornA
            self.CoulShift[0] = 1/self.Cut
            self.CoulShift[1] = 1
            self.CoulShift[2] = self.Cut

            self.UShift[0] = -(erf(fac*self.Cut))/ self.Cut #-np.exp(-self.Kappa[0] * self.Cut) / self.Cut
            self.UShift[1] = np.exp(-fac*fac*self.Cut*self.Cut)/self.BornA[0]/self.BornA[0]
            self.UShift[2] = -np.exp(-fac*fac*self.Cut*self.Cut)*2.0/self.BornA[0]**3*(1.0 - fac**2 * self.Cut**2)
        else:
            self.UShift = 0.

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        BornA = None
        Coef = None
        if not self.FixedCoef: Coef = 1. / self.Sys.Units.FPE
        if not self.FixedBornA: BornA = self.BornA.ArgAvg
        self.SetParam(BornA, Coef)
        


 
