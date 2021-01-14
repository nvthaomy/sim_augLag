#/usr/bin/env python


import base
import numpy as np


class GravityField(base.PotentialClass):
    """Gravitational field interactions."""

    Names = ["gravityfield"]    

    Source = """
>>> defs
float PlaneDistMass
float ThisU

>>> argmainloopbeforepair
PlaneDistMass = (POSIMINIMAGE([[PlaneAxis]]) - PlaneLoc) * MASSI
ARGVAL = PlaneDistMass
ARGTYPE = 0
[ARGGET]

>>> argeval
PlaneDistMass = ARGVAL
THISU = GravConst(0) * PlaneDistMass
if (CalcDUParam) then
    DU_GravConst(0) = PlaneDistMass
endif
[UPDATEENERGY]

>>> mainloopbeforepair
PlaneDistMass = (POSIMINIMAGE([[PlaneAxis]]) - PlaneLoc) * MASSI
THISU = GravConst(0) * PlaneDistMass
[UPDATEENERGY]
if (CalcDUParam) then
    DU_GravConst(0) = PlaneDistMass
endif
if (CalcForce) then
    FORCEI = 0.d0
    FORCEI([[PlaneAxis]]) = -MASSI * GravConst(0)
    Force(i,:) = Force(i,:) + Forcei
endif
""" 

    TestArgs = {"GravConst":1.1, "PlaneAxis":1, "PlaneLoc":-1000.0}
    Type = base.ptypes.FieldPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False,
                 GravConst = None, PlaneAxis = 0, PlaneLoc = -1000):
        """Initializes a gravitational field potential."""
        if PlaneAxis >= Sys.World.Dim:
            raise ValueError("PlaneAxis must be in range 0 to %d for this system." % (Sys.World.Dim-1))
        Source = GravityField.Source.replace("[[PlaneAxis]]", "%d" % PlaneAxis)
        base.PotentialClass.__init__(self, Sys, Label = Label, 
                                     Filter = Filter, 
                                     Source = Source)
        self.ModuleVars += ["PlaneLoc"]

        if GravConst is None:
            raise ValueError("Must specify GravConst.")
        
        self.Param.Add("GravConst", 1, Value = GravConst, Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale)
        self.PlaneLoc = PlaneLoc
 
        self.Update()      
        

    def SetParam(self, GravConst = None):
        """Sets parameters for this potential."""
        self.GravConst = GravConst
        self.Update()

    def Val(self, x, Types = None):
        """Returns the value of the potential."""
        if Types is None:
            Mass = 1.0
        else:
            Mass = Types[0].Mass
        return self.GravConst[0] * Mass * x
        
    def DVal(self, x, Types = None):
        """Returns the derivative of the potential."""
        if Types is None:
            Mass = 1.0
        else:
            Mass = Types[0].Mass
        return self.GravConst[0] * Mass      

    def Update(self):
        """Updates the potential."""
        pass

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.GravConst.MaxChange = MaxFracChange * self.Sys.Units.EScale / self.Sys.Units.LScale

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        pass

        

                
            
        
        

        