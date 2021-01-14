#/usr/bin/env python


import numpy as np
import base


class PairSpline(base.SplinePotentialClass):
    """Pairwise spline interactions."""

    Names = ["pairspline", "spline", "bspline", "pairbspline"]    

    Source = base.basefortran.SourcePairTemplate + """
>>> pairenergy
SPCoef = 1.
Thisx = DIJ
[spindex]
[spenergy]

>>> pairvirial
[spvirial]

>>> pairduparam
[spduparam]

>>> pairdwparam
[sppairdwparam]
"""

    Type = base.ptypes.PairPotential
    KnotsShiftable = False
    KnotsCyclic = False


    TestArgs = {"Cut":2.5, "NKnot":6}    
    

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Knots = None, NKnot = 20,
                 NonbondEneSlopeInit = "50kTperA", BondEneSlope = "200kTperA"):
        """Initializes a B-spline pair potential."""
        base.SplinePotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                           Filter = Filter, Fixed = Fixed,
                                           Knots = Knots, NKnot = NKnot,
                                           InnerBC = 2, OuterBC = 3, 
                                           Source = PairSpline.Source)                                  
        #default spline treatments
        if self.Filter.Bonded: 
            self.EneSlopeInner = BondEneSlope 
            self.ForceInner = True
            self.EneSlopeOuter = BondEneSlope   
            self.ForceOuter = True
        else:
            #constrain the slope of the inner core
            self.EneSlopeInnerInit = NonbondEneSlopeInit
            self.EneSlopeInnerMin = 0.
            
        #make the minimum histogram requirement a little less stringent
        self.MinHistFrac = 0.001            
                       
        #update the potential
        self.Update()        


    def Estimate(self, MinKnot = None, MaxKnot = None, NonbondEneSlope = "50kTperA"):
        """Estimates spline potential parameters based on arg counts and histogram."""
        if self.Filter.Bonded:
            base.SplinePotentialClass.Estimate(MinKnot=MinKnot, MaxKnot=MaxKnot)
        else:
            KnotHist = self.KnotHist()
            if np.any(KnotHist):
                ind = np.nonzero(KnotHist)[0][0]
                r = self.Cut * (ind + 2.0) / self.NKnot
                DeltaEne = self.Sys.Units.ConvertEne(NonbondEneSlope, self.Sys.TempSet) * r / (ind + 1.0)
                for i in range(ind+1):
                    self.Knots[i] = (ind-i) * DeltaEne
                self.Knots[ind+1:] = 0.
                self.Update()                            


    def EmulateLJ(self, Sigma = 1., Epsilon = 1., Shift = False, MaxKnot = None):
        """Sets knot values to roughly emulate a Lennard-Jones potential."""
        if Shift:
            ShiftVal = -4. * ((self.Cut/Sigma)**(-12) - (self.Cut/Sigma)**(-6))
        else:
            ShiftVal = 0.
        tVals = np.arange(0, self.NKnot - 0.1, 0.25) 
        DistVals = tVals * self.SPDist
        DistVals[0] = self.SPDist*0.5
        iVals = tVals.astype(int)
        tVals = tVals - iVals
        LJVals = (4. * Epsilon * ((DistVals/Sigma)**(-12) -
                                  (DistVals/Sigma)**(-6)) + ShiftVal)
        #get rid of strongly repulsive core and just do linear ramp
        for i in range(len(LJVals) - 1, -1, -1):
            if LJVals[i] > 10*Epsilon and i < self.NKnot - 2:
                LJVals[i] = 2*LJVals[i+1] - LJVals[i+2]
        #filter for maximum energy
        if not MaxKnot is None:
            for (i,v) in enumerate(LJVals):
                if v > MaxKnot:
                    LJVals[i] = MaxKnot
        #iterate to match
        Weights = (LJVals / Epsilon + 2)
        Weights = 1. / (Weights + 1)**3
        Weights = Weights / Weights.max()
        #least squares determination
        self.FitSpline(DistVals, LJVals, Weights, MaxKnot=MaxKnot)           


    def SetupTest(self):
        """Sets parameters for a test run."""
        self.EmulateLJ(Sigma = 0.9, Epsilon = 1.1)
        

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = 2. * MaxFracChange * self.Sys.Units.EScale
            #make the first knot have a much bigger change because it will have a high value
            self.Knots.MaxChange[0] = 50. * MaxFracChange * self.Sys.Units.EScale
            if self.Filter.Bonded:
                #make the last knot have a much bigger change because it will have a high value
                self.Knots.MaxChange[-1] = 50. * MaxFracChange * self.Sys.Units.EScale


        

        