#/usr/bin/env python


import numpy as np
import base


class LJGaussian(base.PotentialClass):
    """Lennard-Jones plus Gaussian interactions."""

    Names = ["ljg"]    

    Source = base.basefortran.SourcePairTemplate + """
>>> defs
float idist2
float idist6
float idist12
float val1
float val2
float Sig
float Eps
float val3
float val4
float val5
float val6
float val7

>>> pairenergy
Sig = Sigma(0)
Eps = Epsilon(0)
idist2 = Sig**2 / DIJSQ
idist6 = idist2 * idist2 * idist2
idist12 = idist6 * idist6
val1 = 4.d0 * (idist12 - idist6)
val2 = 24.d0 * idist6 - 48.d0 * idist12
val3 = (DIJ - Dist0(0))
val4 = val3*val3
val5 = exp(-Kappa(0) * val4)
THISU = Eps * val1 + B(0) * val5 + UShift(0)

>>> pairvirial
THISW = Eps * val2 - 2.d0 * B(0) * Kappa(0) * DIJ * val3 * val5

>>> pairduparam
DU_Epsilon(0) = val1 + UShift(1)
val6 = 1.d0 / Sig
DU_Sigma(0) = -Eps * val2 * val6 + UShift(2)           
DDU_Sigma_Sigma(0,0) = Epsilon(0) * (528.d0*idist12 - 120.d0*idist6) * val6*val6 + UShift(3)
DDU_Sigma_Epsilon(0,0) = -val2 * val6 + UShift(4)
DU_B(0) = val5 + UShift(5)
val6 = -val4 * val5
DU_Kappa(0) = B(0) * val6 + UShift(6)
DDU_B_Kappa(0,0) = val6 + UShift(7)
val6 = 2.d0 * Kappa(0) * val3 * val5
DU_Dist0(0) = B(0) * val6 + UShift(8)
DDU_B_Dist0(0,0) = val6 + UShift(9)
DDU_Kappa_Kappa(0,0) = B(0) * val4*val4 * val5 + UShift(10)
val6 = 2.d0*B(0) * val5 
DDU_Kappa_Dist0(0,0) = val6 * val3 * (1.d0 - val4*Kappa(0)) + UShift(11)
DDU_Dist0_Dist0(0,0) = val6 * Kappa(0) * (2.d0 * Kappa(0) * val4 - 1.d0) + UShift(12)

>>> pairdwparam
DW_Epsilon(0) = val2
val6 = 1.d0 / Sig
val7 = (144.d0 * idist6 - 576.d0 * idist12 ) * val6   
DW_Sigma(0) = Eps *  val7    
DDW_Sigma_Sigma(0,0) = Eps * (720.d0*idist6 - 6336.d0*idist12) * val6*val6 
DDW_Sigma_Epsilon(0,0) = val7
val7 = 2.d0 * DIJ * val3 * val5
DW_B(0) = -val7 * Kappa(0)
val6 = val4 * Kappa(0) 
DW_Kappa(0) = B(0) * val7 * (val6 - 1.d0)
DDW_Kappa_Kappa(0,0) = -B(0) * val7 * val4 * (val6 - 2.d0)
DDW_B_Kappa(0,0) = val7 * (val6 - 1.d0)
val7 = 2.d0 * DIJ * val5
DW_Dist0(0) = -val7 * B(0) * Kappa(0) * (2.d0 * Kappa(0) * val4 - 1.d0)
DDW_Dist0_Kappa(0,0) = B(0) * val7 * (1.d0 - 5.d0*val6 + 2.d0*val6*val6)
DDW_Dist0_Dist0(0,0) = -2.d0 * val7 * B(0) * Kappa(0)*Kappa(0) * val3 * (2.d0 * val6 - 3.d0)
DDW_Dist0_B(0,0) = -val7 * Kappa(0) * (2.d0 * val6 - 1.d0)
""" 

    TestArgs = {"Sigma":0.9, "Epsilon":1.2, "B":0.5,
                "Kappa":0.2, "Dist0":1.0, "Cut":1.235, "Shift":False} 
    Type = base.ptypes.PairPotential

    def __init__(self, Sys, Label = "", Cut = None, Filter = None,
                 Fixed = False,
                 Sigma = None, Epsilon = None, 
                 B = None, Kappa = None, Dist0 = None,
                 Shift = False):
        """Initializes a Lennard-Jones plus Gaussian potential."""
        base.PotentialClass.__init__(self, Sys, Label = Label, Cut = Cut,
                                     Filter = Filter, Source = LJGaussian.Source)
        self.ModuleVars += ["UShift"]
                                     
        self.Shift = Shift
        self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)
        self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.Param.Add("B", 1, Value = B, Min = -100., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)
        self.Param.Add("Kappa", 1, Value = Kappa, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale**2)
        self.Param.Add("Dist0", 1, Value = Dist0, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.UShift = np.zeros(13, dtype=float)
 
        self.Update()      
        

    def SetParam(self, Sigma = None, Epsilon = None,
                 B = None, Kappa = None, Dist0 = None,):
        """Sets parameters for this potential."""
        if not Epsilon is None:  self.Epsilon = Epsilon
        if not Sigma is None: self.Sigma = Sigma
        if not B is None: self.B = B
        if not Kappa is None: self.Kappa = Kappa
        if not Dist0 is None: self.Dist0 = Dist0
        self.Update()

    def Val(self, x):
        """Returns the value of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0:
            y = x / self.Sigma[0]
            ELJ = 4.*self.Epsilon[0]*(y**(-12) - y**(-6)) 
            y = x - self.Dist0[0]
            EGauss = self.B[0] * np.exp(-self.Kappa[0] * y**2)
            return ELJ + EGauss + self.UShift[0]
        else:
            return 1.e300
            
    def DVal(self, x):
        """Returns the derivative of the potential."""
        if x >= self.Cut:
            return 0.
        elif x > 0:
            y = x / self.Sigma[0]
            DELJ = self.Epsilon[0] * (-48.*y**(-13) + 24.*y**(-7)) / self.Sigma[0] 
            y = x - self.Dist0[0]
            DEGauss = -2. * self.B[0] * self.Kappa[0] * y *  np.exp(-self.Kappa[0] * y**2)
            return DELJ + DEGauss
        else:
            return -1.e300            

    def Update(self):
        """Updates the potential."""
        if self.Shift:
            idist2 = (self.Sigma[0] / self.Cut)**2
            idist6 = idist2 * idist2 * idist2
            idist12 = idist6 * idist6
            val1 = 4. * (idist12 - idist6)
            val2 = 24. * idist6 - 48. * idist12
            val3 = self.Cut - self.Dist0[0]
            val4 = val3*val3
            val5 = np.exp(-self.Kappa[0] * val4)
            self.UShift[0] = self.Epsilon[0] * val1 + self.B[0] * val5
            self.UShift[1] = val1 
            val6 = 1. / self.Sigma[0]
            self.UShift[2] = -self.Epsilon[0] * val2 * val6         
            self.UShift[3] = self.Epsilon[0] * (528.*idist12 - 120.*idist6) * val6*val6
            self.UShift[4] = -val2 * val6
            self.UShift[5] = val5
            val6 = -val4 * val5
            self.UShift[6] = self.B[0] * val6
            self.UShift[7] = val6
            val6 = 2. * self.Kappa[0] * val3 * val5
            self.UShift[8] = self.B[0] * val6
            self.UShift[9] = val6
            self.UShift[10] = self.B[0] * val4*val4 * val5
            val6 = 2.*self.B[0] * val5 
            self.UShift[11] = val6 * val3 * (1. - val4*self.Kappa[0])
            self.UShift[12] = val6 * self.Kappa[0] * (2. * self.Kappa[0] * val4 - 1.)
            self.UShift = -self.UShift
        else:
            self.UShift[:] = 0.

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Sigma.SoftMin = 0.8 * self.Arg.Min
        self.Epsilon.SoftMin = 0.01 * self.Sys.Units.EScale
        self.B.SoftMin = 0.01 * self.Sys.Units.EScale
        self.Kappa.SoftMin = (0.01 * self.Sys.Units.LScale)**2
        self.Dist0.SoftMin = 0.01 * self.Sys.Units.LScale
        if not MaxFracChange is None:
            self.Sigma.MaxChange = self.Sys.Units.LScale * MaxFracChange
            self.Epsilon.MaxChange = self.Sys.Units.EScale * MaxFracChange

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Sigma = self.Arg.Min * 2**(1./6.)
        Epsilon = self.Epsilon
        self.SetParam(Sigma, Epsilon)

        

                
            
        
        

        
