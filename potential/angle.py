#/usr/bin/env python

import numpy as np
import base


class Angle(base.PotentialClass):
    """Harmonic angle interactions."""

    Names = ["angle"]    

    Source = """
>>> defs
float dTheta

>>> angleenergydparam
dTheta = THETA - Theta0(0) 
ThisU = FConst(0) * dTheta*dTheta
[UPDATEENERGY]
if (CalcDUParam) then
    DU_FConst(0)  = dTheta * dTheta
    DU_Theta0(0)  = -2. * FConst(0) * dTheta
    DDU_Theta0_Theta0(0,0) = 2. * FConst(0)
    DDU_FConst_Theta0(0,0) = -2. * dTheta
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
    ThisdU = 2. * FConst(0) * dTheta 
    Forcei =  (ThisdU / (dij * sintheta)) * (costheta * rij / dij + rjk / djk)
    Forcek = -(ThisdU / (djk * sintheta)) * (costheta * rjk / djk + rij / dij)
    Force(i,:) = Force(i,:) + Forcei
    Force(j,:) = Force(j,:) - Forcei - Forcek
    Force(k,:) = Force(k,:) + Forcek
endif
"""

    TestArgs = {"Theta0":1.5, "FConst":3.1}
    Type = base.ptypes.AnglePotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False,
                 Theta0 = None, FConst = None):
        """Initializes a harmonic angle interaction."""
        base.PotentialClass.__init__(self, Sys, Label, Filter = Filter, 
                                     Source = Angle.Source)
        if Theta0 is None:
            Theta0 = Sys.Units.LScale
        if FConst is None:
            FConst = Sys.Units.EScale / Sys.Units.LScale**2
        self.Param.Add("Theta0", 1, Value = Theta0, Fixed = Fixed, Scale = 1., 
                       Min = 0., SoftMin = 0.001,
                       Max = np.pi, SoftMax = np.pi * 0.99999)
        self.Param.Add("FConst", 1, Value = FConst, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.EScale)

        self.Update()        

    def SetParam(self, Theta0, FConst):
        """Sets parameters for this potential."""
        self.Theta0 = Theta0
        self.FConst = FConst
        self.Update()

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.FConst.SoftMin = 0.0001 * self.Sys.Units.EScale / self.Sys.Units.LScale**2
        if not MaxFracChange is None:
            self.FConst.MaxChange = MaxFracChange * self.Sys.Units.EScale / self.Sys.Units.LScale**2
            self.Theta0.MaxChange = MaxFracChange * np.pi / 4. 
        
    def Val(self, x):
        """Returns the value of the potential."""
        dTheta = x - self.Theta0[0] 
        return self.FConst[0] * dTheta**2
        
    def DVal(self, x):
        """Returns the derivative of the potential."""
        dTheta = x - self.Theta0[0]  
        return 2. * self.FConst[0] * dTheta      

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        Theta0 = self.Arg.Avg
        FConst = self.Sys.Units.EScale
        self.SetParam(Theta0, FConst)
        


                
            
        
        

        