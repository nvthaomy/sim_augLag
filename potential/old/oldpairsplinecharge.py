#/usr/bin/env python

import numpy as np
import base

import sim.chem as chem



class OldPairSplineCharge(base.SplinePotentialClass):
    """Pairwise spline interactions multiplied by atom charges."""

    Names = ["oldpairsplinecharge"]    

    Source = """
>>> defs
float SPt
int SPInd
float SPx
float SPdm1
float SPd0
float SPdp1
float SPdp2
float val0
float val1
float val2
float val3
float val4
float val5
float val6
float Chargei
float Chargej
float ChargeCoef


>>> pairsplinechargeenergydparam
SPx = DIJ * SPiDist
SPind = max(min(int(SPx), [[NKnot]] - 1), 0)
SPt = SPx - real(SPind)
ChargeCoef = Chargei * Chargej * Coef
val4 = (SPC0(SPInd) + SPt * (SPC1(SPInd) + SPt * (SPC2(SPInd) + SPt * SPC3(SPInd))))
val5 = (SPC1(SPInd) + SPt * (2.*SPC2(SPInd) + SPt * 3.*SPC3(SPInd))) * SPx
THISU = val4 * ChargeCoef
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    THISW = val5 * ChargeCoef
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM) then
    SPdm1 = 0.166666666666667d0 + SPt * (-0.5d0 + SPt * (0.5d0 - SPt * 0.166666666666667d0))
    SPd0 = 0.666666666666667d0 + SPt * SPt * (-1.0d0 + SPt * 0.5d0)
    SPdp1 = 0.166666666666667d0 + SPt * (0.5d0 + SPt * (0.5d0 - SPt * 0.5d0))
    SPdp2 = SPt * SPt * SPt * 0.166666666666667d0
    if (SPind == 0) then
        !note: y(-1) = 2 y(0) - y(1) for linearity condition
        val0 = 0.d0
        val1 = SPd0 + 2.0d0 * SPdm1
        val2 = SPdp1 - SPdm1
        val3 = SPdp2
    elseif (SPind == [[NKnot]] - 1) then
        !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
        val0 = SPdm1
        val1 = SPd0 - 0.5d0 * SPdp1 + SPdp2
        val2 = 0.d0
        val3 = 0.d0
    elseif (SPind == [[NKnot]] - 2) then
        !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
        val0 = SPdm1
        val1 = SPd0
        val2 = SPdp1 - 0.5d0 * SPdp2
        val3 = 0.d0
    else
        val0 = SPdm1
        val1 = SPd0
        val2 = SPdp1
        val3 = SPdp2
    endif
    DU_Knots(SPInd - 1) = val0 * ChargeCoef
    DU_Knots(SPInd    ) = val1 * ChargeCoef
    DU_Knots(SPInd + 1) = val2 * ChargeCoef
    DU_Knots(SPInd + 2) = val3 * ChargeCoef
    DDU_Charge_Knots(AIDI, SPInd - 1) = val0 * Chargej 
    DDU_Charge_Knots(AIDI, SPInd    ) = val1 * Chargej
    DDU_Charge_Knots(AIDI, SPInd + 1) = val2 * Chargej 
    DDU_Charge_Knots(AIDI, SPInd + 2) = val3 * Chargej 
    DDU_Charge_Knots(AIDJ, SPInd - 1) = val0 * Chargei 
    DDU_Charge_Knots(AIDJ, SPInd    ) = val1 * Chargei 
    DDU_Charge_Knots(AIDJ, SPInd + 1) = val2 * Chargei 
    DDU_Charge_Knots(AIDJ, SPInd + 2) = val3 * Chargei 
    DU_Charge(AIDI) = val4 * Chargej
    DU_Charge(AIDJ) = val4 * Chargei
    if (AIDI == AIDJ) then
        DDU_Charge_Charge(AIDI, AIDJ) = val4 * 2.d0
    else
        DDU_Charge_Charge(AIDI, AIDJ) = val4
    endif
endif
if (CALCDWPARAM) then
    SPdm1 = (-0.5d0 + SPt * (1.d0 - 0.5d0 * SPt)) * SPx
    SPd0 = (SPt * (-2.0d0 + SPt * 1.5d0)) * SPx
    SPdp1 = (0.5d0 + SPt * (1.d0 - 1.5d0 * SPt)) * SPx
    SPdp2 = (SPt * SPt * 0.5d0) * SPx
    if (SPind == 0) then
        !note: y(-1) = 2 y(0) - y(1) for linearity condition
        val0 = 0.d0
        val1 = SPd0 + 2.0d0 * SPdm1
        val2 = SPdp1 - SPdm1
        val3 = SPdp2
    elseif (SPind == [[NKnot]] - 1) then
        !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
        val0 = SPdm1
        val1 = SPd0 - 0.5d0 * SPdp1 + SPdp2
        val2 = 0.d0
        val3 = 0.d0
    elseif (SPind == [[NKnot]] - 2) then
        !note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
        val0 = SPdm1
        val1 = SPd0
        val2 = SPdp1 - 0.5d0 * SPdp2
        val3 = 0.d0
    else
        val0 = SPdm1
        val1 = SPd0
        val2 = SPdp1
        val3 = SPdp2
    endif
    DW_Knots(SPInd - 1) = val0 * ChargeCoef
    DW_Knots(SPInd    ) = val1 * ChargeCoef
    DW_Knots(SPInd + 1) = val2 * ChargeCoef
    DW_Knots(SPInd + 2) = val3 * ChargeCoef
    DDW_Charge_Knots(AIDI, SPInd - 1) = val0 * Chargej
    DDW_Charge_Knots(AIDI, SPInd    ) = val1 * Chargej
    DDW_Charge_Knots(AIDI, SPInd + 1) = val2 * Chargej
    DDW_Charge_Knots(AIDI, SPInd + 2) = val3 * Chargej
    DDW_Charge_Knots(AIDJ, SPInd - 1) = val0 * Chargei
    DDW_Charge_Knots(AIDJ, SPInd    ) = val1 * Chargei
    DDW_Charge_Knots(AIDJ, SPInd + 1) = val2 * Chargei
    DDW_Charge_Knots(AIDJ, SPInd + 2) = val3 * Chargei
    DW_Charge(AIDI) = val5 * Chargej
    DW_Charge(AIDJ) = val5 * Chargei
    if (AIDI == AIDJ) then
        DDW_Charge_Charge(AIDI, AIDJ) = val5 * 2.d0
    else
        DDW_Charge_Charge(AIDI, AIDJ) = val5
    endif
endif

>>> mainloopbeforepair
Chargei = Charge(AIDI)

>>> mainlooppair
Chargej = Charge(AIDJ)
[pairsplinechargeenergydparam]
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
Chargei = Charge(AIDI)
Chargej = Charge(AIDJ)

>>> argeval
DIJ = ARGVAL
DIJSQ = DIJ * DIJ
[pairsplinechargeenergydparam]
"""

    Type = base.ptypes.PairPotential
    KnotsShiftable = False
    KnotsCyclic = False
    UsesCharge = True
    UsesATypes = True

    TestArgs = {"Cut":2.5, "Knots":[5., 3., 2., 1., 0.]}    
    

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Knots = None, NKnot = 20, Coef = None):
        """Initializes a B-spline pair potential multiplied by charges."""
        NType, TypeAIDs, TypeLabels = Sys.World.GetPairTypes()

        base.SplinePotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                           Filter = Filter, Fixed = Fixed,
                                           Knots = Knots, NKnot = NKnot,
                                           InnerBC = 2, OuterBC = 3,
                                           NType = NType, TypeLabels = TypeLabels,
                                           Source = OldPairSplineCharge.Source)   
        self.ModuleVars += ["Coef"]
                             
        #coefficient behavior
        if Coef is None:
            Coef = 1. / Sys.Units.FPE
        self.Coef = Coef
        #default spline treatment
        self.ConstrainSpline = False
        #update the potential
        self.Update()   
        self.SetTypeInd()
            

    def SetParam(self, Knots):
        """Sets parameters for this potential."""
        try:
            self.Knots = Knots
        except ValueError:
            raise ValueError("Improper value for Knots. Must be single float or length-%d" % self.NKnot)
        self.Update()

    def EmulateCoulomb(self, Shift = True):
        """Sets knot values to roughly emulate a Coulomb potential."""
        import sim.mathutil
        #This fitting doesn't really work -- results in an oscillatory potential
        if Shift:
            ShiftVal = -1./self.Cut
        else:
            ShiftVal = 0.
        tVals = np.arange(0.1, self.NKnot - 0.1, 0.25) #HERE 
        DistVals = tVals * self.SPDist
        iVals = tVals.astype(int)
        tVals = tVals - iVals
        CoulVals = 1. / (DistVals + 1.e-300) + ShiftVal
        #get rid of steep core and just do linear ramp
        for i in range(len(CoulVals) - 3, -1, -1):
            if CoulVals[i] > 3*CoulVals[i+1]:
                CoulVals[i] = 2*CoulVals[i+1] - CoulVals[i+2]
        #least squares determination
        self.FitSpline(DistVals, CoulVals)            
       
    def GetTypeInd(self, Type1, Type2):
        """Returns the type number for this potential."""
        return chem.PairID(Type1.AID, Type2.AID)      
        
    def SetTypeInd(self, TypeInd = None):
        """Sets the type to be used in Val and DVal, or None for generic."""
        if TypeInd is None:
            self.CurrentChargeCoef = self.Coef * self.Sys.Units.ECharge**2
        else:
            AID1, AID2 = chem.GetPair(TypeInd)
            Charge1 = self.Sys.World.AtomTypes[AID1].Charge 
            Charge2 = self.Sys.World.AtomTypes[AID2].Charge 
            self.CurrentChargeCoef = self.Coef * Charge1 * Charge2           

    def KnotInd(self, x):
        """Returns the knot index for this potential."""
        if x >= self.Cut:
            return None
        return min(int(x * self.SPiDist), self.NKnot - 1)    
        
    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        return self.CurrentChargeCoef * (self.SPC0[i] + t * (self.SPC1[i] + t * (self.SPC2[i] + t * self.SPC3[i])))
        
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        return self.CurrentChargeCoef * (self.SPC1[i] + t * (2.*self.SPC2[i] + t * 3.*self.SPC3[i])) * self.SPiDist

    def DUParam(self, x):
        """Returns the values of the derivative with respect to spline knots."""
        if x >= self.Cut:
            return 0.
        x = x * self.SPiDist
        i = min(int(x), self.NKnot - 1)
        t = x - float(i)
        d = np.zeros(self.NKnot, float)
        SPdm1 = 0.166666666666667 + t * (-0.5 + t * (0.5 - t * 0.166666666666667))
        SPd0 = 0.666666666666667 + t * t * (-1.0 + t * 0.5)
        SPdp1 = 0.166666666666667 + t * (0.5 + t * (0.5 - t * 0.5))
        SPdp2 = t * t * t * 0.166666666666667
        if i == 0:
            #note: y(-1) = 2 y(0) - y(1) for linearity condition
            d[i    ] = SPd0 + 2.0 * SPdm1
            d[i + 1] = SPdp1 - SPdm1
            d[i + 2] = SPdp2
        elif i == self.NKnot - 1:
            #note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
            d[i - 1] = SPdm1
            d[i    ] = SPd0 - 0.5 * SPdp1 + SPdp2
        elif i == self.NKnot - 2:
            #note: y(N) = -0.5 y(N-1) and y(N+1) = y(N-1) for zero slope and value condition
            d[i - 1] = SPdm1
            d[i    ] = SPd0
            d[i + 1] = SPdp1 - 0.5 * SPdp2
        else:
            d[i - 1] = SPdm1
            d[i    ] = SPd0
            d[i + 1] = SPdp1
            d[i + 2] = SPdp2
        return d

    def Update(self):
        """Updates the potential."""
        #check the cutoff
        if self.Cut is None:
            raise TypeError("The cutoff distance for PairSpline %s not defined." % self.Label)
        #compute the distances
        self.SPDist = self.Cut / float(self.NKnot)
        self.SPiDist = 1. / self.SPDist        
        #make a shortcut array
        Y = np.zeros(self.NKnot+3, float)
        #fill in the y-values
        Y[1:-2] = self.Knots
        #make Y(0) for zero zecond deriv
        Y[0] = 2.*Y[1] - Y[2]
        #make the last two knots for zero val and first deriv
        Y[-2] = -0.5 * Y[-3]
        Y[-1] = Y[-3]
        #make the coefficients
        C = np.zeros((4, self.NKnot), float)
        for i in range(1, self.NKnot + 1):
            C[:,i-1] = np.dot(self.BFunc, Y[i-1:i+3])
        self.SPC0 = C[0]
        self.SPC1 = C[1]
        self.SPC2 = C[2]
        self.SPC3 = C[3]

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = 2. * MaxFracChange * self.Sys.Units.EScale
            #make the first knot have a much bigger change because it will have a high value
            self.Knots.MaxChange[0] = 50. * MaxFracChange * self.Sys.Units.EScale
            if self.Filter.Bonded:
                #make the last knot have a much bigger change because it will have a high value
                self.Knots.MaxChange[-1] = 50. * MaxFracChange * self.Sys.Units.EScale


         