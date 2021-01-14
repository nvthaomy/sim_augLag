#/usr/bin/env python


import numpy as np
import base



class TorsionSpline(base.SplinePotentialClass):
    """Spline torsion angle interactions."""

    Names = ["torsionspline"]    

    Source = """
>>> torsionenergydparam
Thisx = phi
SPCoef = 1.
[spindex]
[spenergy]
[UPDATEENERGY]
if (CalcDUParam) then
    [spduparam]
endif

>>> argtorsionloop
ARGVAL = Phi
ARGTYPE = 0
[ARGGET]

>>> argeval
Phi = ARGVAL
[torsionenergydparam]

>>> torsionloop
if (TorDefined) then
    [torsionenergydparam]
    if (CalcForce) then
        [spdenergy]
        Forcei =  (djk * ThisdU) * crossijk / crosssqijk
        Forcel = -(djk * ThisdU) * crossjkl / crosssqjkl
        Forcej = (dot_product(rkl,rjk) / djk**2) * Forcel &
             & - (dot_product(rij,rjk) / djk**2 + 1.d0) * Forcei
        Forcek = (dot_product(rij,rjk) / djk**2) * Forcei &
             & - (dot_product(rkl,rjk) / djk**2 + 1.d0) * Forcel	
        Force(i,:) = Force(i,:) + Forcei
        Force(j,:) = Force(j,:) + Forcej
        Force(k,:) = Force(k,:) + Forcek
        Force(l,:) = Force(l,:) + Forcel
    endif
endif
"""
    
    TestArgs = {"Knots" : [-5, -3, 0, 1, -2]}
    
    Type = base.ptypes.TorsionPotential
    KnotsShiftable = True
    KnotsCyclic = True
    

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, Knots = None, NKnot = 20):
        """Initializes a spline torsion angle interaction."""
        base.SplinePotentialClass.__init__(self, Sys, Label = Label, 
                                           Filter = Filter, Fixed = Fixed,
                                           Knots = Knots, NKnot = NKnot,
                                           Min = -np.pi, Max = np.pi,
                                           Source = TorsionSpline.Source)            
        self.Update()    


    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.Knots.MaxChange = MaxFracChange * self.Sys.Units.EScale
        

  
   
    
    