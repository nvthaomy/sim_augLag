#/usr/bin/env python


import base
import numpy as np


class LJExponents(base.PotentialClass):
    """Lennard-Jones interactions with variable exponents."""

    Names = ["ljexponents"]    

    Source = base.basefortran.SourcePairTemplate + """
>>> defs
float idist
float idist1
float idist2
float val1
float val2
float val3
float val4
float Sig
float Eps

>>> pairenergy
Sig = Sigma(0)
Eps = Epsilon(0)
idist = Sig / sqrt(DIJSQ)
idist1 = idist**(Exp1)
idist2 = idist**(Exp2)
val1 = 4.d0 * (idist1 - idist2)
val2 = 4.d0 * (Exp2*idist2 - Exp1*idist1)
THISU = Eps * val1 + UShift(0)

>>> pairvirial
THISW = Eps * val2

>>> pairduparam
DU_Epsilon(0) = val1 + UShift(1)
val1 = 1.d0 / Sig
DU_Sigma(0) = -Eps * val2 * val1 + UShift(2)   
val3 = log(idist)
DU_Exp1(0) = 4.d0 * Eps * idist1 * val3 + UShift(3)
DU_Exp2(0) = -4.d0 * Eps * idist2 * val3 + UShift(4)
DDU_Sigma_Sigma(0,0) = Epsilon(0) * 4.d0 * (Exp1*(Exp1-1.d0)*idist1 - Exp2*(Exp2-1.d0)*idist2) * val1*val1 + UShift(5)
DDU_Sigma_Epsilon(0,0) = -val2 * val1 + UShift(6)
DDU_Sigma_Exp1(0,0) = 4.d0 * Eps * idist1 * val1 * (1.d0 + val3*Exp1) + UShift(7)
DDU_Sigma_Exp2(0,0) = -4.d0 * Eps * idist2 * val1 * (1.d0 + val3*Exp2) + UShift(8)
DDU_Epsilon_Exp1(0,0) = 4.d0 * idist1 * val3 + UShift(9)
DDU_Epsilon_Exp2(0,0) = -4.d0 * idist2 * val3 + UShift(10)
DDU_Exp1_Exp1(0,0) = 4.d0 * Eps * idist1 * val3 * val3 + UShift(11)
DDU_Exp2_Exp2(0,0) = -4.d0 * Eps * idist2 * val3 * val3 + UShift(12)

>>> pairdwparam
DW_Epsilon(0) = val2
val1 = 1.d0 / Sig
val4 = 4.d0 * (Exp2*Exp2 * idist2 - Exp1*Exp1 * idist1 ) * val1 
DW_Sigma(0) = Eps *  val4    
DDW_Sigma_Sigma(0,0) = Eps * 4.d0 * (Exp2*Exp2*(Exp2-1)*idist2 - Exp1*Exp1*(Exp1-1)*idist1) * val1*val1 
DDW_Sigma_Epsilon(0,0) = val4
val3 = log(idist)
val4 = -4.d0 * idist1 * (1.d0 + Exp1 * val3)
DW_Exp1(0) = Eps * val4
DDW_Epsilon_Exp1(0,0) = val4
val4 = 4.d0 * idist2 * (1.d0 + Exp2 * val3)
DW_Exp2(0) = Eps * val4
DDW_Epsilon_Exp2(0,0) = val4
val4 = -4.d0 * Eps * idist1 * (2.d0 + val3*Exp1)
DDW_Sigma_Exp1(0,0) =  Exp1 * val4 * val1
DDW_Exp1_Exp1(0,0) = val4 * val3 
val4 = 4.d0 * Eps * idist2 * (2.d0 + val3*Exp2)
DDW_Sigma_Exp2(0,0) = Exp2 * val4 * val1
DDW_Exp2_Exp2(0,0) = val4 * val3
""" 

    TestArgs = {"Sigma":1.1, "Epsilon":1.2, "Cut":2.5, "Shift":True, "Exp1":18, "Exp2":8}
    Type = base.ptypes.PairPotential

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False,
                 Sigma = None, Epsilon = None, Exp1 = 12, Exp2 = 6,
                 Shift = False):
        """Initializes a Lennard-Jones potential."""
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Filter, Source = LJExponents.Source)
        self.ModuleVars = ["UShift"]

        self.Shift = Shift
        if Sigma is None:
            Sigma = Sys.Units.LScale
        if Epsilon is None:
            Epsilon = Sys.Units.EScale
        self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)
        self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.Param.Add("Exp1", 1, Value = Exp1, Min = 3., SoftMin = 4.,
                       Max = 100., SoftMax = 90., Fixed = Fixed)
        self.Param.Add("Exp2", 1, Value = Exp2, Min = 3., SoftMin = 4.,
                       Max = 100., SoftMax = 90., Fixed = Fixed)                       
        self.UShift = np.zeros(13, dtype=float)
 
        self.Update()      
        

    def SetParam(self, Sigma = None, Epsilon = None, Exp1 = None, Exp2 = None):
        """Sets parameters for this potential."""
        if not Epsilon is None:  self.Epsilon = Epsilon
        if not Sigma is None: self.Sigma = Sigma
        if not Exp1 is None: self.Exp1 = Exp1
        if not Exp2 is None: self.Exp2 = Exp2        
        self.Update()

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0.:
            x = x / self.Sigma[0]
            return 4.*self.Epsilon[0]*(x**(-self.Exp1[0]) - x**(-self.Exp2[0])) + self.UShift[0]
        else:
            return 1.e300
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0.:
            x = x / self.Sigma[0]
            return -4.*self.Epsilon[0] * (self.Exp1[0] * x**(-self.Exp1[0]-1) 
                                        - self.Exp2[0] * x**(-self.Exp2[0]-1)) / self.Sigma[0]
        else:
            return -1.e300

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            Exp1 = self.Exp1[0]
            Exp2 = self.Exp2[0]
            Cut = self.Cut / self.Sigma
            val3 = -np.log(Cut)
            self.UShift[0] = -4.*self.Epsilon * (Cut**(-Exp1)-Cut**(-Exp2))
            self.UShift[1] = -4.*(Cut**(-Exp1)-Cut**(-Exp2))
            self.UShift[2] = -self.Epsilon * 4. * (Exp1*Cut**(-Exp1) - Exp2*Cut**(-Exp2)) / self.Sigma
            self.UShift[3] = -4.*self.Epsilon * Cut**(-Exp1) * val3
            self.UShift[4] = 4.*self.Epsilon * Cut**(-Exp2) * val3
            self.UShift[5] = -self.Epsilon * 4. * (Exp1*(Exp1-1.)*Cut**(-Exp1) - Exp2*(Exp2-1.)*Cut**(-Exp2)) / self.Sigma**2
            self.UShift[6] = -4.*(Exp1*Cut**(-Exp1) - Exp2*Cut**(-Exp2)) / self.Sigma
            self.UShift[7] = -4.*self.Epsilon * Cut**(-Exp1) * (1. + val3*Exp1) / self.Sigma
            self.UShift[8] = 4.*self.Epsilon * Cut**(-Exp2) * (1. + val3*Exp2) / self.Sigma
            self.UShift[9] = -4.* Cut**(-Exp1) * val3
            self.UShift[10] = 4.* Cut**(-Exp2) * val3
            self.UShift[11] = -4.*self.Epsilon * Cut**(-Exp1) * val3 * val3
            self.UShift[12] = 4.*self.Epsilon * Cut**(-Exp2) * val3 * val3
        else:
            self.UShift = 0.

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Sigma.SoftMin = 0.8 * self.Arg.Min
        self.Epsilon.SoftMin = 0.01 * self.Sys.Units.EScale
        if not MaxFracChange is None:
            self.Sigma.MaxChange = self.Sys.Units.LScale * MaxFracChange
            self.Epsilon.MaxChange = self.Sys.Units.EScale * MaxFracChange
            self.Exp1.MaxChange = 2. * MaxFracChange
            self.Exp2.MaxChange = 2. * MaxFracChange

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Sigma = self.Arg.Min * 2**(1./6.)
        Epsilon = self.Epsilon
        self.SetParam(Sigma, Epsilon)

        

                
            
        
        

        