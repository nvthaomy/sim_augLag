#/usr/bin/env python

import numpy as np
import base

import sim.chem as chem



class PairSplineCharge(base.SplinePotentialClass):
    """Pairwise spline interactions multiplied by atom charges."""

    Names = ["pairsplinecharge"]    

    Source = """
>>> defs
float val4
float val5
float val6
float Chargei
float Chargej
float ChargeiCoef
float ChargejCoef


>>> pairsplinechargeenergydparam
Thisx = DIJ
SPCoef = Chargei * Chargej * Coef
[spindex]
[spenergy]
[UPDATEENERGY]
if (CALCVIRIAL .or. CALCFORCE .or. CALCDWPARAM) then
    [spvirial]
    [UPDATEVIRIAL]
endif
if (CALCDUPARAM .or. CALCDWPARAM) then
    ChargeiCoef = Chargei * Coef
    ChargejCoef = Chargej * Coef
endif
if (CALCDUPARAM) then
    [spduparam]
    if (SPdm1 .ne. 0.) DDU_Charge_Knots(AIDI, SPInd - 1) = SPdm1 * ChargejCoef 
    if (SPd0  .ne. 0.) DDU_Charge_Knots(AIDI, SPInd    ) = SPd0  * ChargejCoef
    if (SPdp1 .ne. 0.) DDU_Charge_Knots(AIDI, SPInd + 1) = SPdp1 * ChargejCoef
    if (SPdp2 .ne. 0.) DDU_Charge_Knots(AIDI, SPInd + 2) = SPdp2 * ChargejCoef
    if (SPdm1 .ne. 0.) DDU_Charge_Knots(AIDJ, SPInd - 1) = SPdm1 * ChargeiCoef
    if (SPd0  .ne. 0.) DDU_Charge_Knots(AIDJ, SPInd    ) = SPd0  * ChargeiCoef
    if (SPdp1 .ne. 0.) DDU_Charge_Knots(AIDJ, SPInd + 1) = SPdp1 * ChargeiCoef
    if (SPdp2 .ne. 0.) DDU_Charge_Knots(AIDJ, SPInd + 2) = SPdp2 * ChargeiCoef
    DU_Charge(AIDI) = SPThisU * ChargejCoef
    DU_Charge(AIDJ) = SPThisU * ChargeiCoef
    if (AIDI == AIDJ) then
        DDU_Charge_Charge(AIDI, AIDJ) = SPThisU * 2.d0 * Coef
    else
        DDU_Charge_Charge(AIDI, AIDJ) = SPThisU * Coef
    endif
endif
if (CALCDWPARAM) then
    [sppairdwparam]
    if (SPdm1 .ne. 0.) DDW_Charge_Knots(AIDI, SPInd - 1) = SPdm1 * ChargejCoef 
    if (SPd0  .ne. 0.) DDW_Charge_Knots(AIDI, SPInd    ) = SPd0  * ChargejCoef
    if (SPdp1 .ne. 0.) DDW_Charge_Knots(AIDI, SPInd + 1) = SPdp1 * ChargejCoef
    if (SPdp2 .ne. 0.) DDW_Charge_Knots(AIDI, SPInd + 2) = SPdp2 * ChargejCoef
    if (SPdm1 .ne. 0.) DDW_Charge_Knots(AIDJ, SPInd - 1) = SPdm1 * ChargeiCoef
    if (SPd0  .ne. 0.) DDW_Charge_Knots(AIDJ, SPInd    ) = SPd0  * ChargeiCoef
    if (SPdp1 .ne. 0.) DDW_Charge_Knots(AIDJ, SPInd + 1) = SPdp1 * ChargeiCoef
    if (SPdp2 .ne. 0.) DDW_Charge_Knots(AIDJ, SPInd + 2) = SPdp2 * ChargeiCoef
    DW_Charge(AIDI) = SPThisW * ChargejCoef
    DW_Charge(AIDJ) = SPThisW * ChargeiCoef
    if (AIDI == AIDJ) then
        DDW_Charge_Charge(AIDI, AIDJ) = SPThisW * 2.d0 * Coef
    else
        DDW_Charge_Charge(AIDI, AIDJ) = SPThisW * Coef
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
                                           Source = PairSplineCharge.Source)         
        self.ModuleVars += ["Coef"]
                             
        #coefficient behavior
        if Coef is None:
            Coef = 1. / Sys.Units.FPE
        self.Coef = Coef

        #update the potential
        self.Update()   
        self.SetTypeInd()
            

    def EmulateCoulomb(self, Shift = True):
        """Sets knot values to roughly emulate a Coulomb potential."""
        #This fitting doesn't really work -- results in an oscillatory potential
        if Shift:
            ShiftVal = -1./self.Cut
        else:
            ShiftVal = 0.
        tVals = np.arange(0.1, self.NKnot - 0.1, 0.25)  
        DistVals = tVals * self.SPDist[0]
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
            self.SPCoef = self.CurrentChargeCoef


    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = 2. * MaxFracChange * self.Sys.Units.EScale
            #make the first knot have a much bigger change because it will have a high value
            self.Knots.MaxChange[0] = 50. * MaxFracChange * self.Sys.Units.EScale
            if self.Filter.Bonded:
                #make the last knot have a much bigger change because it will have a high value
                self.Knots.MaxChange[-1] = 50. * MaxFracChange * self.Sys.Units.EScale


         