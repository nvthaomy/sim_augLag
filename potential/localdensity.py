#/usr/bin/env python

import numpy as np
import base


class LocalDensity(base.SplinePotentialClass):
    """Local density interaction."""

    Names = ["LocalDensity"]

    Source = """
>>> defs
float rho
float val1
float WTerm
float DPhiDRDivR

>>> init
if (CALCFORCE .or. CALCVIRIAL .or. CALCDWPARAM) DUDRho = 0.d0
if (MODE==0) then
    LocalDensity = 0.d0
else
    !subtract the current local density energy off because we will recalculate the whole thing
    PEnergy = PEnergy - THISPENERGY
    !zero this energy term
    THISPENERGY = 0.
endif

>>> mainlooppair
if (dijsq < InnerCut(1)) then
    rho = 1.d0
else
    rho = LDCoef(0) + dijsq * (LDCoef(1) + dijsq * (LDCoef(2) + LDCoef(3) * dijsq))
endif
if (Applyi) LocalDensity(I) = LocalDensity(I) + rho * SCALE
if (Applyj) LocalDensity(J) = LocalDensity(J) + rho * SCALE

>>> localdensityenergyduparam
SPCoef = 1.
Thisx = rho
[spindex]
[spenergy]
[UPDATEENERGY]
if (CalcDUParam) then
    [spduparam]
endif

>>> localdensityloop1
rho = LocalDensity(I)
[localdensityenergyduparam]
if (CALCFORCE .or. CALCVIRIAL .or. CALCDWPARAM) then
    [spdenergy]
    DUDRho(I) = ThisDU
endif

>>> localdensityloop2
if (dijsq >= InnerCut(1)) then
    DPhiDRDivR = LDCoef(4) + dijsq * (LDCoef(5) + dijsq * LDCoef(6))
    ThisDU = 0.d0
    if (Applyi) ThisDU = ThisDU + DUDRho(I)
    if (Applyj) ThisDU = ThisDU + DUDRho(J)
    if (CALCVIRIAL .or. CALCDWPARAM) then
        THISW = ThisDU * DPhiDRDivR * dijsq
        [UPDATEVIRIAL]
    endif
    if (CALCFORCE) then
        ForceI = ThisDU * DPhiDRDivR * rij 
        Force(I,:) = Force(I,:) + ForceI
        Force(J,:) = Force(J,:) - ForceI
    endif  
    if (CALCDWPARAM) then
        do K = 1, 2
            if (K==1 .and. Applyi) then
                rho = LocalDensity(I)
            elseif (K==2 .and. Applyj) then
                rho = LocalDensity(J)
            endif
            Thisx = rho 
            SPCoef = DPhiDRDivR * DIJ
            [spindex]
            [sppairdwparam]   
        enddo
    endif
endif

>>> saveenergystate
OldLocalDensity = LocalDensity

>>> revertenergystate
LocalDensity = OldLocalDensity

>>> arginit
LocalDensity = 0.d0

>>> argmainlooppair
[mainlooppair]

>>> arglocaldensityloop1
ARGVAL = LocalDensity(I)
ARGTYPE = 0
[ARGGET]

>>> argeval
rho = ARGVAL
[localdensityenergyduparam]
if (CalcDWParam) then
    print *, "Local density potential not set up to do dW/dParam using histograms."
    stop
endif
"""

    TestArgs = {"Cut":1.5, "InnerCut":0.5, "Knots":[0., 18., 4., 20., 5.], "RhoMin":10., "RhoMax":10.5}
    
    Type = base.ptypes.LocalDensityPotential
    KnotsShiftable = False
    KnotsCyclic = False
    NoEnergyUpdateScaling = True


    def __init__(self, Sys, Label = "", Cut = None, Filter = None, Fixed = False,
                 Knots = None, NKnot = 20, RhoMin = 0., RhoMax = 12., InnerCut = 0.):
        """Initializes a local density potential."""
        base.SplinePotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                           Filter = Filter, Fixed = Fixed,
                                           Knots = Knots, NKnot = NKnot,
                                           InnerBC = 1, OuterBC = 1,
                                           Min = RhoMin, Max = RhoMax,
                                           Source = LocalDensity.Source)
        
        #add module vars
        self.ModuleVars += ["LocalDensity", "OldLocalDensity", "DUDRho", "RhoMin", "RhoMax",
                            "InnerCut", "LDCoef"]
        
        ##set the default outer slope and inner value
        self.EneSlopeOuter = 0.
        self.EneSlopeInner = 0.
    
        #check to see if range of densities is specified
        if RhoMax == None or RhoMin == None:
            raise ValueError("Must specify RhoMin and RhoMax for LocalDensityPotential.")
        self.RhoMin = RhoMin
        self.RhoMax = RhoMax
        
        #set potential reporting options
        self.SetReport(ArgMin = self.RhoMin, ArgMax = self.RhoMax)
        
        #set inner cutoff
        self.InnerCut = np.ones(2, dtype=float) * InnerCut
        
        #set the coefficients for the local density indicator function
        self.LDCoef = np.zeros(7, dtype=float)

        self.Update()


    def PreLoad(self):
        #make the local density arrays
        NAtom = self.Sys.NAtom
        self.LocalDensity = np.zeros((NAtom,), dtype=float)
        self.OldLocalDensity = np.zeros((NAtom,), dtype=float)
        self.DUDRho = np.zeros((NAtom,), dtype=float) 
        base.SplinePotentialClass.PreLoad(self)      
        
        
    def RhoVal(self, r):
        """Returns the value of the density function."""
        if r >= self.Cut: 
            return 0.
        elif r <= self.InnerCut[0]:
            return 1.
        LDC = self.LDCoef
        rsq = r*r
        return LDC[0] + rsq * (LDC[1] + rsq * (LDC[2] + LDC[3] * rsq))
        

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = MaxFracChange * self.Sys.Units.EScale
            

    def SetCut(self, Cut = None, InnerCut = None):
        """Sets the cutoffs for the LD potential."""
        if not Cut is None:
            self.Cut = Cut
        if not InnerCut is None:
            self.InnerCut[0] = InnerCut
        self.Update()


    def Update(self):
        """Updates the potential."""
        #check the cutoff
        if self.Cut is None:
            raise TypeError("The cutoff distance for LocalDensityPotential %s not defined." % self.Label)
        if self.InnerCut[0] >= self.Cut:
            raise ValueError("The inner cutoff must be less than the cutoff distance.")
        #compute inner cut sq
        self.InnerCut[1] = self.InnerCut[0]**2
        #compute the coefficients
        Cut1Sq = self.InnerCut[0]**2
        Cut2Sq = self.Cut**2
        d = (Cut2Sq - Cut1Sq)**3
        self.LDCoef[0] = (Cut2Sq**3 - 3.*Cut1Sq * Cut2Sq**2) / d
        self.LDCoef[1] = 6.*Cut1Sq*Cut2Sq / d
        self.LDCoef[2] = -3.*(Cut1Sq + Cut2Sq) / d
        self.LDCoef[3] = 2. / d
        self.LDCoef[4] = 2. * self.LDCoef[1]
        self.LDCoef[5] = 4. * self.LDCoef[2]
        self.LDCoef[6] = 6. * self.LDCoef[3]
        #compute the spline parameters
        base.SplinePotentialClass.Update(self)
 
