#/usr/bin/env python


import numpy as np
import base


class AngleSpline(base.SplinePotentialClass):
    """Spline angle interactions."""

    Names = ["anglespline"]    

    Source = """
>>> angleenergydparam
Thisx = Theta
SPCoef = 1.
[spindex]
[spenergy]
[UPDATEENERGY]
if (CalcDUParam) then
    [spduparam]
endif

>>> argangleloop
ARGVAL = Theta
ARGTYPE = 0
[ARGGET]

>>> argeval
Theta = ARGVAL
[angleenergydparam]

>>> angleloop
[angleenergydparam]
if (CalcForce .and. sintheta /= 0.d0) then
    [spdenergy]
    Forcei =  (ThisdU / (dij * sintheta)) * (costheta * rij / dij + rjk / djk)
    Forcek = -(ThisdU / (djk * sintheta)) * (costheta * rjk / djk + rij / dij)
    Force(i,:) = Force(i,:) + Forcei
    Force(j,:) = Force(j,:) - Forcei - Forcek
    Force(k,:) = Force(k,:) + Forcek
endif
"""
    
    TestArgs = {"Knots":[0., 10., 4., 20., 5.]} 

    Type = base.ptypes.AnglePotential
    KnotsShiftable = True
    KnotsCyclic = False
    

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False, Knots = None, NKnot = 20,
                 AngleEneSlopeInit = "50kT"):
        """Initializes a spline angle interaction."""
        base.SplinePotentialClass.__init__(self, Sys, Label, Filter = Filter,
                                           Fixed = Fixed, Knots = Knots, NKnot = NKnot,
                                           Min = 0., Max = np.pi, InnerBC = 1, OuterBC = 1,
                                           Source = AngleSpline.Source)           
        #default spline treatments
        self.EneSlopeInnerInit = AngleEneSlopeInit
        self.EneSlopeOuterInit = AngleEneSlopeInit
        self.EneSlopeInnerMin = 0.
        self.EneSlopeOuterMin = 0.        
        self.ForceInner = True
        self.ForceOuter = True
        #update
        self.Update()        


    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = MaxFracChange * self.Sys.Units.EScale


    def SetupTest(self):
        """Sets parameters for a test run."""
        KnotVals = [1.8*np.sin(0.1+5*np.pi*i/float(self.NKnot))
                    for i in range(self.NKnot)]
        self.SetParam(KnotVals)
        
        


                
            
        
        

        