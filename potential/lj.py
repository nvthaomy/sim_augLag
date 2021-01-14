#/usr/bin/env python


import base
import numpy as np


class LJ(base.PotentialClass):
    """Lennard-Jones interactions."""

    Names = ["lj"]    

    Source = base.basefortran.SourcePairTemplate + """
>>> defs
float idist2
float idist6
float idist12
float val1
float val2
float val3
float Sig
float Eps

>>> pairenergy
Sig = Sigma(0)
Eps = Epsilon(0)
idist2 = Sig**2 / DIJSQ
idist6 = idist2 * idist2 * idist2
idist12 = idist6 * idist6
val1 = 4.d0 * (idist12 - idist6)
val2 = 24.d0 * idist6 - 48.d0 * idist12
THISU = Eps * val1 + UShift(0)

>>> pairvirial
THISW = Eps * val2

>>> pairduparam
DU_Epsilon(0) = val1 + UShift(1)
val1 = 1.d0 / Sig
DU_Sigma(0) = -Eps * Val2 * val1 + UShift(2)           
DDU_Sigma_Sigma(0,0) = Eps * (528.d0*idist12 - 120.d0*idist6) * val1*val1 + UShift(3)
DDU_Sigma_Epsilon(0,0) = -val2 * val1 + UShift(4)

>>> pairdwparam
DW_Epsilon(0) = val2
val1 = 1.d0 / Sig
val3 = (144.d0 * idist6 - 576.d0 * idist12 ) * val1    
DW_Sigma(0) = Eps *  val3    
DDW_Sigma_Sigma(0,0) = Eps * (720.d0*idist6 - 6336.d0*idist12) * val1*val1 
DDW_Sigma_Epsilon(0,0) = val3
""" 

    TestArgs = {"Sigma":0.8, "Epsilon":1.2, "Cut":2.5, "Shift":True}
    Type = base.ptypes.PairPotential

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False,
                 Sigma = None, Epsilon = None,
                 Shift = False):
        """Initializes a Lennard-Jones potential."""
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Filter, Source = LJ.Source)
        self.ModuleVars += ["UShift"]

        self.Shift = Shift
        if Sigma is None:
            Sigma = Sys.Units.LScale
        if Epsilon is None:
            Epsilon = Sys.Units.EScale
        self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)
        self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.UShift = np.zeros(5, dtype=float)
 
        self.Update()      
        

    def SetParam(self, Sigma = None, Epsilon = None):
        """Sets parameters for this potential."""
        if not Epsilon is None:  self.Epsilon = Epsilon
        if not Sigma is None: self.Sigma = Sigma
        self.Update()

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0.:
            x = x / self.Sigma[0]
            return 4.*self.Epsilon[0]*(x**(-12) - x**(-6)) + self.UShift[0]
        else:
            return 1.e300
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0.:
            x = x / self.Sigma[0]
            return self.Epsilon[0]*(-48. * x**(-13) + 24. * x**(-7)) / self.Sigma[0]
        else:
            return -1.e300            

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            Cut = self.Cut / self.Sigma
            self.UShift[0] = -4.*self.Epsilon * (Cut**(-12)-Cut**(-6))
            self.UShift[1] = -4.*(Cut**(-12)-Cut**(-6))
            self.UShift[2] = -self.Epsilon * (48.*Cut**(-12) - 24.*Cut**(-6)) / self.Sigma
            self.UShift[3] = -self.Epsilon * (528.*Cut**(-12) - 120.*Cut**(-6)) / self.Sigma**2
            self.UShift[4] = -(48.*Cut**(-12) - 24.*Cut**(-6)) / self.Sigma
        else:
            self.UShift[:] = 0.

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Sigma.SoftMin = 0.8 * self.Arg.Min
        self.Epsilon.SoftMin = 0.01 * self.Sys.Units.EScale
        if not MaxFracChange is None:
            self.Sigma.MaxChange = self.Sys.Units.LScale * MaxFracChange
            self.Epsilon.MaxChange = self.Sys.Units.EScale * MaxFracChange

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Sigma = self.Arg.Min * 2**(1./6.)
        Epsilon = self.Epsilon
        self.SetParam(Sigma, Epsilon)

        

                
            
        
        

        