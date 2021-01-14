
import base
import numpy as np



class RgUmbrella(base.PotentialClass):
    """Umbrella potential in Rg."""

    Names = ["rgumbrella"]    

    Source = """
>>> defs
float RgPosSum__(Dim)
float RgPosSumSq__(Dim)
float RgPosAvg__(Dim)
float RgN__
float Rg__
float RgLastPos__(Dim)
float ThisPos(Dim)
float dPos(Dim)

>>> rgumbrellaenergydparam
ThisU = 0.5 * FConst(0) * (Rg__ - RgSet(0))**2 
[UPDATEENERGY]
if (CalcDUParam) then
    DU_FConst(0) = 0.5 * (Rg__ - RgSet)**2 
    DU_RgSet(0) = FConst * (RgSet - Rg__)
    DDU_RgSet_RgSet(0,0) = FConst
    DDU_RgSet_FConst(0,0) = RgSet - Rg__
endif

>>> arginit
[init]

>>> argpreloopbeforepair
[preloopbeforepair]

>>> argbeforemainloop
RgPosAvg__ = RgPosSum__ / RgN__
Rg__ = sqrt(sum(RgPosSumSq__) / RgN__ - sum(RgPosAvg__*RgPosAvg__))
ARGVAL = Rg__
ARGTYPE = 0
[ARGGET]

>>> argeval
Rg__ = ARGVAL
[rgumbrellaenergydparam]

>>> init
!initialize Rg averages
RgPosSum__ = 0.d0 
RgPosSumSq__ = 0.d0  
RgLastPos__ = 0.d0  
RgN__ = 0. 

>>> preloopbeforepair
!minimum image; compute distance from last particle
dPos = POSI - RgLastPos__
!minimum-image
dPos = dPos - BOXL * dnint(dPos * IBOXL)
!compute minimum-imaged position
ThisPos = RgLastPos__ + dPos
!add to running sums
RgPosSum__ = RgPosSum__ + ThisPos
RgPosSumSq__ = RgPosSumSq__ + ThisPos*ThisPos
!update last position
RgLastPos__ = ThisPos
RgN__ = RgN__ + 1.0

>>> beforemainloop
RgPosAvg__ = RgPosSum__ / RgN__
Rg__ = sqrt(sum(RgPosSumSq__) / RgN__ - sum(RgPosAvg__*RgPosAvg__))
[rgumbrellaenergydparam]
RgLastPos__ = 0.d0

>>> mainloopbeforepair
if (CalcForce) then
    !minimum image; compute distance from last particle
    dPos = POSI - RgLastPos__
    !minimum-image
    dPos = dPos - BOXL * dnint(dPos * IBOXL)
    !compute minimum-imaged position
    ThisPos = RgLastPos__ + dPos
    !do force calculation
    ThisForce = -(FConst(0) / RgN__) * (1.0 - RgSet(0)/Rg__) * (ThisPos - RgPosAvg__) 
    Force(i,:) = Force(i,:) + ThisForce
    !update last position
    RgLastPos__ = ThisPos
endif
""" 

    TestArgs = {"RgSet":5, "FConst":5}
    Type = base.ptypes.GlobalPotential

    def __init__(self, Sys, Label = "", Filter = None,
                 Fixed = False, RgSet = None, FConst = 0):
        """Initializes a harmonic umbrella potential in Rg."""
        Source = RgUmbrella.Source
        base.PotentialClass.__init__(self, Sys, Label = Label, 
                                     Filter = Filter, 
                                     Source = Source)
        if RgSet is None:
            raise ValueError("Must specify RgSet.")
        if FConst is None:
            raise ValueError("Must specify FConst.")
            
        self.Param.Add("RgSet", 1, Value = RgSet, Fixed = Fixed,
                       Scale = self.Sys.Units.LScale)
        self.Param.Add("FConst", 1, Value = FConst, Fixed = Fixed,
                       Scale = self.Sys.Units.EScale / self.Sys.Units.LScale**2)                       
 
        self.Update()      
        
        #set potential reporting options
        self.SetReport(ArgMin =0.0, ArgMax = self.Sys.Units.LScale * 10.)
        

    def SetParam(self, RgSet, FConst):
        """Sets parameters for this potential."""
        self.RgSet = RgSet
        self.FConst = FConst
        self.Update()

    def Val(self, x):
        """Returns the value of the potential.
Here, x is equal to Rg"""
        return 0.5 * self.FConst[0] * (x - self.RgSet[0])**2
        
    def DVal(self, x):
        """Returns the derivative of the potential.
Here, x is equal to Rg"""
        return -self.FConst[0] * (x - self.RgSet[0])   

    def Update(self):
        """Updates the potential."""
        pass

    def SetBounds(self, MaxFracChange = None):
        """Sets bounds on parameters based on argument ranges and scales."""
        if not MaxFracChange is None:
            self.RgSet.MaxChange = MaxFracChange * self.Sys.Units.LScale
            self.FConst.MaxChange = MaxFracChange * self.Sys.Units.EScale / self.Sys.Units.LScale**2

    def Estimate(self):
        """Estimates parameters based on argument ranges and scales."""
        pass
