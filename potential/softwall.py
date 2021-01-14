#/usr/bin/env python


import base
import numpy as np


class SoftWall(base.PotentialClass):
    """Wall with inverse power repulsion."""

    Names = ["softwall"]    

    Source = """
>>> defs
float PlaneDist
float val1
float val2

>>> softwallenergydparam
val1 = abs(PlaneDist / Sigma)
THISU = Epsilon * (val1**(-Exponent))
[UPDATEENERGY]
if (CalcDUParam) then
    val2 = log(val1)
    DU_Sigma(0) = THISU * Exponent / Sigma
    DDU_Sigma_Sigma(0,0) = THISU * Exponent * (Exponent-1) / (Sigma * Sigma)
    DU_Exponent(0) = -THISU * val2
    DDU_Exponent_Exponent(0,0) = THISU * val2*val2
    DDU_Exponent_Sigma(0,0) = THISU * (1.d0 - val2 * Exponent) / Sigma   
endif

>>> argmainloopbeforepair
PlaneDist = POSI([[PlaneAxis]]) - PlaneLoc
if (BOXL([[PlaneAxis]]) > 0.d0) PlaneDist = PlaneDist - BOXL([[PlaneAxis]]) * dnint(PlaneDist * IBOXL([[PlaneAxis]]))
ARGVAL = PlaneDist
ARGTYPE = 0
[ARGGET]

>>> argeval
PlaneDist = ARGVAL
[softwallenergydparam]

>>> mainloopbeforepair
!interactions for field potential %(PName)s
PlaneDist = POSI([[PlaneAxis]]) - PlaneLoc
if (BOXL([[PlaneAxis]]) > 0.d0) PlaneDist = PlaneDist - BOXL([[PlaneAxis]]) * dnint(PlaneDist * IBOXL([[PlaneAxis]]))
[softwallenergydparam]
if (CalcForce) then
    THISFORCE = 0.d0
    THISFORCE([[PlaneAxis]]) = Exponent * (THISU / PlaneDist)
    Force(i,:) = Force(i,:) + ThisForce
endif
""" 

    TestArgs = {"Sigma":0.52, "Epsilon":2.33, "Exponent":2, "FixedExp":False,
                "PlaneAxis":0, "PlaneLoc":0.1}
    Type = base.ptypes.FieldPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False, FixedExp = True, Shift = True,
                 Sigma = 1., Epsilon = 1., Exponent = 9.,
                 PlaneAxis = 0, PlaneLoc = 0):
        """Initializes a soft wall field potential."""
        if PlaneAxis >= Sys.World.Dim:
            raise ValueError("PlaneAxis must be in range 0 to %d for this system." % (Sys.World.Dim-1))
        Source = SoftWall.Source.replace("[[PlaneAxis]]", "%d" % PlaneAxis)
        base.PotentialClass.__init__(self, Sys, Label = Label,
                                     Filter = Filter, 
                                     Source = Source)
        self.ModuleVars += ["PlaneLoc"]
               
        self.Param.Add("Epsilon", 1, Value = Epsilon, Min = 0., Fixed = True,
                       Scale = self.Sys.Units.EScale)
        self.Param.Add("Sigma", 1, Value = Sigma, Min = 0., Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.Param.Add("Exponent", 1, Value = Exponent, Min = 1., SoftMin = 3.,
                       Max = 100., SoftMax = 90., Fixed = FixedExp)
        self.PlaneLoc = PlaneLoc
 
        self.Update()      
        

    def SetParam(self, Sigma = None, Epsilon = None, Exponent = None):
        """Sets parameters for this potential."""
        if not Epsilon is None:  self.Epsilon = Epsilon
        if not Sigma is None: self.Sigma = Sigma
        if not Exponent is None: self.Exponent = Exponent
        self.Update()

    def Val(self, x):
        """Returns the value of the potential.
Here, x is equal to PlaneDist"""
        x = x / self.Sigma[0]
        if x==0:
            return 1.e300
        else:
            return self.Epsilon[0] * x**(-self.Exponent[0])
            
    def DVal(self, x):
        """Returns the derivative of the potential.
Here, x is equal to PlaneDist"""
        x = x / self.Sigma[0]
        if x==0:
            return 1.e300
        else:
            return -self.Exponent[0] * self.Epsilon[0] * x**(-self.Exponent[0] - 1)            

    def Update(self):
        """Updates the potential."""
        pass

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        self.Sigma.SoftMin = 0.5 * self.Arg.Min
        if not MaxFracChange is None:
            self.Sigma.MaxChange = self.Sys.Units.LScale * MaxFracChange
            self.Exponent.MaxChange = 2. * MaxFracChange

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        pass

        

                
            
        
        

        