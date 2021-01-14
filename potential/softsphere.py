#/usr/bin/env python


import numpy as np
import base


class SoftSphere(base.PotentialClass):
    """Soft sphere interactions."""

    Names = ["softsphere", "ss"]    

    SourceEven = base.basefortran.SourcePairTemplate + """
>>> defs
float idist2
float val1
float val2

>>> pairenergy
idist2 = Sigma(0)*Sigma(0) / DIJSQ
val1 = Epsilon(0) * idist2**[[Exponent/2]]
THISU = val1 + UShift(0)

>>> pairvirial
THISW = -[[Exponent]].0d0 * val1

>>> pairduparam
val2 = 1.d0 / Sigma(0)
DU_Sigma(0) = val1 * [[Exponent]].0d0 * val2 + UShift(1)
DDU_Sigma_Sigma(0,0) = val1 * [[DDCoef]].0d0 * val2*val2 + UShift(2)

>>> pairdwparam
val2 = 1.d0 / Sigma(0)
DW_Sigma(0) = -val1 * val2 * ([[Exponent]].0d0)**2
DDW_Sigma_Sigma(0,0) = -val1 * [[DDCoef]].0d0 * [[Exponent]].0d0 * val2*val2 
"""

    SourceOdd = base.basefortran.SourcePairTemplate + """
>>> defs
float idist
float val1
float val2

>>> pairenergy
idist = Sigma(0) / DIJ
val1 = Epsilon(0) * idist**[[Exponent]]
THISU = val1 + UShift(0)

>>> pairvirial
THISW = -[[Exponent]].0d0 * val1

>>> pairduparam
val2 = 1.d0 / Sigma(0)
DU_Sigma(0) = val1 * [[Exponent]].0d0 * val2 + UShift(1)
DDU_Sigma_Sigma(0,0) = val1 * [[DDCoef]].0d0 * val2*val2 + UShift(2)

>>> pairdwparam
val2 = 1.d0 / Sigma(0)
DW_Sigma(0) = -val1 * val2 * ([[Exponent]].0d0)**2
DDW_Sigma_Sigma(0,0) = -val1 * [[DDCoef]].0d0 * [[Exponent]].0d0 * val2*val2 
"""

    SourceVar = base.basefortran.SourcePairTemplate + """
>>> defs
float idist
float val1
float val2
float val3
float Exp1

>>> pairenergy
Exp1 = Exponent(0)
idist = Sigma(0) / DIJ
val1 = Epsilon(0) * idist**Exp1
THISU = val1 + UShift(0)

>>> pairvirial
THISW = -Exp1 * val1

>>> pairduparam
val2 = 1.d0 / Sigma(0)
val3 = log(idist)
DU_Sigma(0) = val1 * Exp1 * val2 + UShift(1)
DDU_Sigma_Sigma(0,0) = val1 * Exp1*(Exp1 - 1.d0) * val2*val2 + UShift(2)
DU_Exponent(0) = val1 * val3 + UShift(3)
DDU_Exponent_Exponent(0,0) = val1 * val3*val3 + UShift(4)
DDU_Exponent_Sigma(0,0) = val1 * val2 * (1.d0 + val3*Exp1) + UShift(5)

>>> pairdwparam
val2 = 1.d0 / Sigma(0)
val3 = log(idist)
DW_Sigma(0) = -val1 * val2 * Exp1**2
DDW_Sigma_Sigma(0,0) = -val1 * val2*val2 * (Exp1 - 1.d0) * Exp1 * Exp1
DW_Exponent(0) = -val1 * (Exp1 * val3 + 1.d0)
DDW_Exponent_Exponent(0,0) = -val1 * val3 * (2.d0 + Exp1 * val3)
DDW_Exponent_Sigma(0,0) = -val1 * val2 * Exp1 * (2.d0 + Exp1 * val3)
"""

    TestArgs = {"Exponent":9, "Sigma":1.1, "Epsilon":1.2, "Cut":2.5, "Shift":True,
                "FixedExp" : False}
    Type = base.ptypes.PairPotential

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False, FixedExp = True,
                 Exponent = 9, Sigma = 1., Epsilon = 1.,
                 Shift = False):
        """Initializes a soft-sphere potential."""       
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut, 
                                     Filter = Filter)
        self.ModuleVars += ["UShift"]
        
        self.FixedExp = FixedExp
        self.Shift = Shift
        
        if self.FixedExp:
            if Exponent % 2 == 0:
                Source = SoftSphere.SourceEven
            else:
                Source = SoftSphere.SourceOdd
            Source = Source.replace("[[Exponent/2]]", "%d" % (Exponent/2))
            Source = Source.replace("[[Exponent]]", "%d" % Exponent)
            Source = Source.replace("[[DDCoef]]", "%d" % (Exponent*(Exponent-1)))
        else:
            Source = SoftSphere.SourceVar
        self.SetSource(Source)
        
        if self.FixedExp:
            self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = True,
                           Scale = self.Sys.Units.EScale)
            self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                           Scale = self.Sys.Units.LScale)
            self.UShift = np.zeros(3, dtype=float)
        else:
            self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = True,
                           Scale = self.Sys.Units.EScale)
            self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                           Scale = self.Sys.Units.LScale)
            self.UShift = np.zeros(7, dtype=float)
            
        self.Param.Add("Exponent", 1, Value = Exponent, Min = 1., SoftMin = 2.,
                       Max = 100., SoftMax = 90., Fixed = Fixed or FixedExp)
                  
        self.Update()       
        

    def SetParam(self, Sigma = None, Epsilon = None, Exponent = None):
        """Sets parameters for this potential."""
        if not Epsilon is None:  self.Epsilon = Epsilon
        if not Sigma is None: self.Sigma = Sigma
        if not Exponent is None:
            if self.FixedExp:
                raise ValueError("The exponent is fixed in this potential.")
            else:
                self.Exponent = Exponent
        self.Update()

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0:
            x = x / self.Sigma[0]
            return self.Epsilon[0] * x**(-self.Exponent[0]) + self.UShift[0]
        else:
            return 1.e300
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0:
            x = x / self.Sigma[0]
            return -self.Epsilon[0] * x**(-self.Exponent[0] - 1) * self.Exponent[0]
        else:
            return -1.e300            

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            Cut = self.Cut / self.Sigma[0]
            self.UShift[0] = -self.Epsilon[0] * Cut**(-self.Exponent[0])
            self.UShift[1] = self.UShift[0] * self.Exponent[0] / self.Sigma[0]
            self.UShift[2] = self.UShift[1] * (self.Exponent[0] - 1.) / self.Sigma[0]
            if not self.FixedExp:
                LogTerm = -np.log(Cut)
                self.UShift[3] = self.UShift[0] * LogTerm
                self.UShift[4] = self.UShift[0] * LogTerm**2
                self.UShift[5] = self.UShift[0] * (1. + self.Exponent[0] * LogTerm) / self.Sigma[0]
        else:
            self.UShift = 0.

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Sigma.SoftMin = 0.8 * self.Arg.Min
        if not MaxFracChange is None:
            self.Sigma.MaxChange = self.Sys.Units.LScale * MaxFracChange
            if not self.FixedExp:
                self.Exponent.MaxChange = 2. * MaxFracChange

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Sigma = self.Arg.Avg
        Epsilon = self.Epsilon
        Epsilon = self.Sys.Units.EScale
        Exponent = None
        if not self.FixedExp: Exponent = 9
        self.SetParam(Sigma, Epsilon, Exponent)        

                
            
        
        

        